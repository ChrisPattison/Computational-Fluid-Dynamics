program main_prg
	use kinds_mod
	use plplotlib_mod
	use case_mod
	use util_mod
	use types_mod
	implicit none
	
	integer,parameter::psize = 400
	type(coeff5),allocatable,dimension(:,:)::tcoeff,centercoeff,upwindcoeff
	type(coeff9),allocatable,dimension(:,:)::quickcoeff
	real(wp),allocatable,dimension(:,:)::T,Told,Tgraph
	real(wp),allocatable,dimension(:,:)::lowres,highres
	integer::i,j,k,m
	logical::quick = .false.
	real(wp)::vol
	real(wp)::residual
	real(wp)::beta = 0.95_wp
	
	allocate(tcoeff(psize,psize))
	allocate(centercoeff(psize,psize))
	allocate(upwindcoeff(psize,psize))
	allocate(quickcoeff(psize,psize))
	allocate(T(psize,psize))
	allocate(Tgraph(psize,psize))
	allocate(Told(psize,psize))
	allocate(highres(psize,psize))
	allocate(lowres(psize,psize))
	
	call zero(upwindcoeff)
	call upwind_convection(upwindcoeff)
	call boundary(upwindcoeff)

	call zero(centercoeff)
	call center_convection(centercoeff)
	call boundary(centercoeff)
	
	tcoeff = upwindcoeff

	call extzero(quickcoeff)
	call quick_convection(quickcoeff)
	call extboundary(quickcoeff)

	T = 0.0_wp
	Told = 0.0_wp
	do m = 1, 100000
		if(quick) then
			call extresidvec(quickcoeff,Told,highres)
		else
			call residvec(centercoeff,Told,highres)
		end if
		call residvec(upwindcoeff,Told,lowres)
		
		tcoeff%b = upwindcoeff%b - beta*(highres - lowres)
		tcoeff(psize,:)%b = upwindcoeff(psize,:)%b ! Don't correct the boundary
		tcoeff(:,psize)%b = upwindcoeff(:,psize)%b
				
		call solve(tcoeff, T)
				
		residual = sum((Told-T)**2)
		if(residual<1e-6_wp) then
			exit
		end if
		Told = T
		if(modulo(m,100)==0) then
			write(*,*) m,residual
		end if
	end do
	
	call solve(upwindcoeff,told)
	
	Tgraph = (T - Told)**2
	Tgraph = Tgraph/(sum(Tgraph)/size(Tgraph))
	
	call setup(device='svg',fileName='fig-%n.svg',figSize=[1200,900])
	call figure()
	call subplot(1,1,1)
	call xylim([0.0_wp,psize*L],[0.0_wp,psize*L])
	call contourf(linspace(0.0_wp,L*psize,psize),linspace(0.0_wp,L*psize,psize),Tgraph,10)
	call colorbar(Tgraph,10)
	call ticks()
	call labels('x [m]','y [m]','')

	call figure()
	call subplot(1,1,1)
	call xylim([0.0_wp,psize*L],[0.0_wp,psize*L])
	call contourf(linspace(0.0_wp,L*psize,psize),linspace(0.0_wp,L*psize,psize),T,10)
	call colorbar(T,10)
	call ticks()
	call labels('x [m]','y [m]','')
	
	call figure()
	
	call subplot(2,1,1)
	call xylim([0.0_wp,psize*L],[0.0_wp,psize*L])
	call contourf(linspace(0.0_wp,L*psize,psize),linspace(0.0_wp,L*psize,psize),T,10)
	call colorbar(T,10)
	call ticks()
	call labels('x [m]','y [m]','') 
   	
   	call subplot(2,1,2)
   	call xylim([0.0_wp,psize*L],mixval(T(200,:)))
   	call plot(linspace(0.0_wp,L*psize,psize),T(200,:), lineWidth=2.0_wp)
   	call ticks()
   	call labels('y [m]','T [K]','')
   	call show()
contains

subroutine GaussSeidel(tcoeff,T)
	type(coeff5),dimension(psize,psize),intent(in)::tcoeff
	real(wp),dimension(psize,psize),intent(inout)::T
	integer::i,j
	real(wp)::ndrho
	
	do i = 1, psize
		do j = 1, psize
			ndrho = tvalue(T,i,j+1)*tcoeff(i,j)%N + tvalue(T,i,j-1)*tcoeff(i,j)%S + &
				tvalue(T,i-1,j)*tcoeff(i,j)%W + tvalue(T,i+1,j) *tcoeff(i,j)%E
			T(i,j) = (tcoeff(i,j)%b-ndrho)/tcoeff(i,j)%P
		end do
	end do
end subroutine

subroutine ExtGaussSeidel(tcoeff,T)
	type(coeff9),dimension(psize,psize),intent(in)::tcoeff
	real(wp),dimension(psize,psize),intent(inout)::T
	integer::i,j
	real(wp)::ndrho
	
	do i = 1, psize
		do j = 1, psize
			ndrho = tvalue(T,i,j+1)*tcoeff(i,j)%N + tvalue(T,i,j-1)*tcoeff(i,j)%S + &
				tvalue(T,i-1,j)*tcoeff(i,j)%W + tvalue(T,i+1,j)*tcoeff(i,j)%E + &
				tvalue(T,i,j+2)*tcoeff(i,j)%NN + tvalue(T,i,j-2)*tcoeff(i,j)%SS + &
				tvalue(T,i-2,j)*tcoeff(i,j)%WW + tvalue(T,i+2,j)*tcoeff(i,j)%EE
			T(i,j) = (tcoeff(i,j)%b-ndrho)/tcoeff(i,j)%P
		end do
	end do
end subroutine

subroutine parsebc(bc, tcoeff)
	type(boundc),intent(in)::bc
	type(coeff5),intent(inout)::tcoeff
	select case(bc%bctype)
		case(0)
!           b = -kval*bc%bcvalue
			tcoeff%b = bc%bcvalue
			tcoeff%P = 1.0_wp
			tcoeff%W = 0.0_wp
			tcoeff%E = 0.0_wp
			tcoeff%S = 0.0_wp
			tcoeff%N = 0.0_wp
		case(1)
			tcoeff%b = tcoeff%b + bc%bcvalue
	end select
end subroutine

subroutine extparsebc(bc, tcoeff)
	type(boundc),intent(in)::bc
	type(coeff9),intent(inout)::tcoeff
	select case(bc%bctype)
		case(0)
!           b = -kval*bc%bcvalue
			tcoeff%b = bc%bcvalue
			tcoeff%P = 1.0_wp
			tcoeff%W = 0.0_wp
			tcoeff%E = 0.0_wp
			tcoeff%S = 0.0_wp
			tcoeff%NN = 0.0_wp
			tcoeff%WW = 0.0_wp
			tcoeff%EE = 0.0_wp
			tcoeff%SS = 0.0_wp
			tcoeff%NN = 0.0_wp
		case(1)
			tcoeff%b = tcoeff%b + bc%bcvalue
	end select
end subroutine

function resid(tcoeff,T) result(v)
	type(coeff5),dimension(psize,psize),intent(in)::tcoeff
	real(wp),dimension(psize,psize),intent(in)::T
	real(wp)::rho
	real(wp)::v
	integer::i,j
	v = 0.0_wp

	do i = 1, psize
		do j = 1, psize
			rho = tvalue(T,i,j+1)*tcoeff(i,j)%N + tvalue(T,i,j-1)*tcoeff(i,j)%S + &
				tvalue(T,i-1,j)*tcoeff(i,j)%W + tvalue(T,i+1,j) *tcoeff(i,j)%E
			v = v + (rho + T(i,j)*tcoeff(i,j)%P - tcoeff(i,j)%b)**2
		end do
	end do
end function

function extresid(tcoeff,T) result(v)
	type(coeff9),dimension(psize,psize),intent(in)::tcoeff
	real(wp),dimension(psize,psize),intent(in)::T
	real(wp)::ndrho
	real(wp)::v
	integer::i,j
	v = 0.0_wp

	do i = 1, psize
		do j = 1, psize
			ndrho = tvalue(T,i,j+1)*tcoeff(i,j)%N + tvalue(T,i,j-1)*tcoeff(i,j)%S + &
				tvalue(T,i-1,j)*tcoeff(i,j)%W + tvalue(T,i+1,j)*tcoeff(i,j)%E + &
				tvalue(T,i,j+2)*tcoeff(i,j)%NN + tvalue(T,i,j-2)*tcoeff(i,j)%SS + &
				tvalue(T,i-2,j)*tcoeff(i,j)%WW + tvalue(T,i+2,j)*tcoeff(i,j)%EE
			v = v + (ndrho + T(i,j)*tcoeff(i,j)%P - tcoeff(i,j)%b)**2
		end do
	end do
end function

subroutine residvec(tcoeff,T,resid)
	type(coeff5),dimension(psize,psize),intent(in)::tcoeff
	real(wp),dimension(psize,psize),intent(in)::T
	real(wp)::rho
	real(wp),dimension(psize,psize),intent(out)::resid
	integer::i,j

	do i = 1, psize
		do j = 1, psize
			rho = tvalue(T,i,j+1)*tcoeff(i,j)%N + tvalue(T,i,j-1)*tcoeff(i,j)%S + &
				tvalue(T,i-1,j)*tcoeff(i,j)%W + tvalue(T,i+1,j) *tcoeff(i,j)%E
			resid(i,j) = rho + T(i,j)*tcoeff(i,j)%P - tcoeff(i,j)%b
		end do
	end do
end subroutine

subroutine extresidvec(tcoeff,T,resid)
	type(coeff9),dimension(psize,psize),intent(in)::tcoeff
	real(wp),dimension(psize,psize),intent(in)::T
	real(wp)::ndrho
	real(wp),dimension(psize,psize),intent(out)::resid
	integer::i,j

	do i = 1, psize
		do j = 1, psize
			ndrho = tvalue(T,i,j+1)*tcoeff(i,j)%N + tvalue(T,i,j-1)*tcoeff(i,j)%S + &
				tvalue(T,i-1,j)*tcoeff(i,j)%W + tvalue(T,i+1,j)*tcoeff(i,j)%E + &
				tvalue(T,i,j+2)*tcoeff(i,j)%NN + tvalue(T,i,j-2)*tcoeff(i,j)%SS + &
				tvalue(T,i-2,j)*tcoeff(i,j)%WW + tvalue(T,i+2,j)*tcoeff(i,j)%EE
			resid(i,j) =  ndrho + T(i,j)*tcoeff(i,j)%P - tcoeff(i,j)%b
		end do
	end do
end subroutine

subroutine solve(tcoeff, T)
	type(coeff5),dimension(psize,psize)::tcoeff
	real(wp),dimension(psize,psize)::T
	integer::k
	
	do k = 1, 1000
		call GaussSeidel(tcoeff,T)
			
		residual = resid(tcoeff,T)
		if(residual < 1e-6_wp) then
			exit
		end if
	end do
end subroutine

function tvalue(T,i,j) result(v)
	real(wp),dimension(:,:),intent(in)::T
	integer,intent(in)::i,j
	real(wp)::v

	if(i <= ubound(T,1).and.i >= lbound(T,1) .and. j <= ubound(T,2) .and. j >= lbound(T,2)) then
		v = T(i,j)
	else
		v = 0.0_wp
	end if
end function

subroutine boundary(tcoeff)
	type(coeff5),dimension(:,:),intent(inout)::tcoeff
	integer::i,j
	do i = 1, psize
		do j = 1, psize
			if(j==psize) then
				vol = vol/2
				call parsebc(NBound, tcoeff(i,j))
			else if(j==1) then
				vol = vol/2
				call parsebc(SBound, tcoeff(i,j))
			end if

			if(i==psize) then
				vol = vol/2
				call parsebc(EBound, tcoeff(i,j))
			else if(i==1) then
				vol = vol/2
				call parsebc(WBound, tcoeff(i,j))
			end if
		end do
	end do
end subroutine

subroutine extboundary(tcoeff)
	type(coeff9),dimension(:,:),intent(inout)::tcoeff
	integer::i,j
	do i = 1, psize
		do j = 1, psize
			if(j==psize) then
				vol = vol/2
				call extparsebc(NBound, tcoeff(i,j))
			else if(j==1) then
				vol = vol/2
				call extparsebc(SBound, tcoeff(i,j))
			end if

			if(i==psize) then
				vol = vol/2
				call extparsebc(EBound, tcoeff(i,j))
			else if(i==1) then
				vol = vol/2
				call extparsebc(WBound, tcoeff(i,j))
			end if
		end do
	end do
end subroutine

subroutine zero(tcoeff)
	type(coeff5),dimension(:,:),intent(inout)::tcoeff
	integer::i,j
	
	do i = 1, psize
		do j = 1, psize
			tcoeff(i,j)%N = 0.0_wp
			tcoeff(i,j)%S = 0.0_wp
			tcoeff(i,j)%E = 0.0_wp
			tcoeff(i,j)%W = 0.0_wp
			tcoeff(i,j)%P = 0.0_wp
			tcoeff(i,j)%b = 0.0_wp
		end do
	end do          
end subroutine

subroutine extzero(tcoeff)
	type(coeff9),dimension(:,:),intent(inout)::tcoeff
	integer::i,j
	
	do i = 1, psize
		do j = 1, psize
			tcoeff(i,j)%N = 0.0_wp
			tcoeff(i,j)%S = 0.0_wp
			tcoeff(i,j)%E = 0.0_wp
			tcoeff(i,j)%W = 0.0_wp
			tcoeff(i,j)%NN = 0.0_wp
			tcoeff(i,j)%SS = 0.0_wp
			tcoeff(i,j)%EE = 0.0_wp
			tcoeff(i,j)%WW = 0.0_wp
			tcoeff(i,j)%P = 0.0_wp
			tcoeff(i,j)%b = 0.0_wp
		end do
	end do          
end subroutine

subroutine diffusion(tcoeff)
	type(coeff5),dimension(:,:),intent(inout)::tcoeff
	integer::i,j

	do i = 1, psize
		do j = 1, psize
			tcoeff(i,j)%N = -kval 
			tcoeff(i,j)%S = -kval
			tcoeff(i,j)%E = -kval
			tcoeff(i,j)%W = -kval
			tcoeff(i,j)%P = kval*4
			tcoeff(i,j)%b = qval*cellvol(tcoeff,i,j)
		end do
	end do
end subroutine

subroutine upwind_convection(tcoeff)
	type(coeff5),dimension(:,:),intent(inout)::tcoeff
	integer:: i,j

	do i = 1, psize
		do j = 1, psize
			tcoeff(i,j)%N = tcoeff(i,j)%N + max(0.0_wp,-rhoUCpVal(2))*L
			tcoeff(i,j)%S = tcoeff(i,j)%S + max(0.0_wp,rhoUCpVal(2))*L
			tcoeff(i,j)%E = tcoeff(i,j)%E + max(0.0_wp,-rhoUCpVal(1))*L
			tcoeff(i,j)%W = tcoeff(i,j)%W + max(0.0_wp,rhoUCpVal(1))*L
			tcoeff(i,j)%P = tcoeff(i,j)%P - (max(0.0_wp,-rhoUCpVal(2))*L &
				+ max(0.0_wp,rhoUCpVal(2))*L + max(0.0_wp,-rhoUCpVal(1))*L & 
				+ max(0.0_wp,rhoUCpVal(1))*L)
		end do
	end do
end subroutine

subroutine center_convection(tcoeff)
	type(coeff5),dimension(:,:),intent(inout)::tcoeff
	integer::i,j

	do i = 1, psize
		do j = 1, psize
			tcoeff(i,j)%N = tcoeff(i,j)%N - rhoUCpVal(2)*L/2.0_wp
			tcoeff(i,j)%S = tcoeff(i,j)%S + rhoUCpVal(2)*L/2.0_wp
			tcoeff(i,j)%E = tcoeff(i,j)%E - rhoUCpVal(1)*L/2.0_wp
			tcoeff(i,j)%W = tcoeff(i,j)%W + rhoUCpVal(1)*L/2.0_wp
		end do
	end do
end subroutine

subroutine quick_convection(tcoeff)
	type(coeff9),dimension(:,:),intent(inout)::tcoeff
	real(wp),dimension(3,3)::interp
	real(wp),dimension(3)::linterp,rinterp
	integer::i,j
	
!	interp = transpose(reshape([1.0_wp, 0.0_wp, 0.0_wp, -1.5_wp, 2.0_wp, -0.5_wp, 0.5_wp, -1.0_wp, 0.5_wp],[3,3]))
!	rinterp = matmul([1.0_wp, 0.5_wp, 0.25_wp], interp)
!	linterp = matmul([1.0_wp, 1.5_wp, 2.25_wp], interp)
	
	do i = 1, psize
		do j = 1, psize
			if (rhoUCpVal(1) > 0.0_wp) then 
				tcoeff(i,j)%WW= tcoeff(i,j)%WW- 0.125_wp*rhoUCpVal(1)*L
				tcoeff(i,j)%W = tcoeff(i,j)%W + 0.875_wp*rhoUCpVal(1)*L 
				tcoeff(i,j)%P = tcoeff(i,j)%P - 0.375_wp*rhoUCpVal(1)*L
				tcoeff(i,j)%E = tcoeff(i,j)%E - 0.375_wp*rhoUCpVal(1)*L

			else
				!TODO: implement
			end if
			
			if (rhoUCpVal(2) > 0.0_wp) then
				tcoeff(i,j)%SS= tcoeff(i,j)%SS- 0.125_wp*rhoUCpVal(2)*L
				tcoeff(i,j)%S = tcoeff(i,j)%S + 0.875_wp*rhoUCpVal(2)*L 
				tcoeff(i,j)%P = tcoeff(i,j)%P - 0.375_wp*rhoUCpVal(2)*L
				tcoeff(i,j)%N = tcoeff(i,j)%N - 0.375_wp*rhoUCpVal(2)*L

			else
				!TODO: implement
			end if
		end do
	end do
end subroutine

function cellvol(tcoeff,i,j) result(v)
	type(coeff5),dimension(:,:),intent(inout)::tcoeff
	integer,intent(in)::i,j
	real(wp)::v

	v = L*L
	if(i==ubound(tcoeff,1).or.i==lbound(tcoeff,1)) then
		v = v/2
	end if
	if(j==ubound(tcoeff,2).or.j==lbound(tcoeff,2)) then
		v = v/2
	end if
end function
end program main_prg
