program main_prg
	use kinds_mod
	use plplotlib_mod
	use case_mod
	use util_mod
	use types_mod
	implicit none
	
	integer,parameter::psize = 100
	type(coeff5),dimension(psize,psize)::tcoeff
	real(wp),dimension(psize,psize)::T
	real(wp),dimension(psize,psize)::b
	real(wp),dimension(psize,4)::A
	integer::i,j,k,m
	real(wp)::ndrho
	real(wp)::vol
	real(wp)::residual
	real(wp)::timestart,timestop

	do i = 1, psize
		do j = 1, psize
			b(i,j) = 0.0_wp
			vol = L*L
			tcoeff(i,j)%N = -kval
			tcoeff(i,j)%S = -kval
			tcoeff(i,j)%E = -kval
			tcoeff(i,j)%W = -kval
			tcoeff(i,j)%P = 0.0_wp
			if(j==psize) then
				tcoeff(i,j)%N = 0.0_wp
				vol = vol/2
				call parsebc(NBound, tcoeff(i,j), b(i,j))
			end if
			if(j==1) then
				tcoeff(i,j)%S = 0.0_wp
				vol = vol/2
				call parsebc(SBound, tcoeff(i,j), b(i,j))
			end if
			if(i==psize) then
				tcoeff(i,j)%E = 0.0_wp
				vol = vol/2
				call parsebc(EBound, tcoeff(i,j), b(i,j))
			end if
			if(i==1) then
				tcoeff(i,j)%W = 0.0_wp
				vol = vol/2
				call parsebc(WBound, tcoeff(i,j), b(i,j))
			end if
			
			if(abs(tcoeff(i,j)%P) < 1.0e-7_wp) then
				tcoeff(i,j)%P = -(tcoeff(i,j)%N + tcoeff(i,j)%S + tcoeff(i,j)%E + tcoeff(i,j)%W)
			end if
			b(i,j) = merge(b(i,j), qval*vol, b(i,j)/=0.0_wp)
		end do
	end do

!===========================================
!=----------Gauss-Seidel-------------------=
!===========================================

	T = 0.0_wp
	call cpu_time(timestart)
	do k = 1, 100000
		do i = 1, psize
			do j = 1, psize
				ndrho = 0.0_wp
				if (abs(tcoeff(i,j)%N) > 1.0e-7_wp) then
					ndrho = ndrho + T(i,j+1) * tcoeff(i,j)%N
				endif
				if (abs(tcoeff(i,j)%S) > 1.0e-7_wp) then
					ndrho = ndrho + T(i,j-1) * tcoeff(i,j)%S
				endif
				if (abs(tcoeff(i,j)%W) > 1.0e-7_wp) then
					ndrho = ndrho + T(i-1,j) * tcoeff(i,j)%W
				endif
				if (abs(tcoeff(i,j)%E) > 1.0e-7_wp) then
					ndrho = ndrho + T(i+1,j) * tcoeff(i,j)%E
				endif
				T(i,j) = (b(i,j)-ndrho)/tcoeff(i,j)%P
			end do
		end do

		residual = resid(tcoeff,T,b)
		if(residual < 1e-6_wp) then
			exit
		end if
	end do
	call cpu_time(timestop)
	write(*,*) "GS Finished in ",k-1," iterations(",timestop-timestart,"s) with residual: ", residual

!===========================================
!=----------------TDMA---------------------=
!===========================================

	T = 0.0_wp
	call cpu_time(timestart)
	do k = 1, 100000
		do j = 1, psize
			do i = 1, psize
				A(i,4) = b(i,j)
				if (abs(tcoeff(i,j)%N) > 1.0e-7_wp) then
					A(i,4) = A(i,4) - T(i,j+1) * tcoeff(i,j)%N
				endif
				if (abs(tcoeff(i,j)%S) > 1.0e-7_wp) then
					A(i,4) = A(i,4) - T(i,j-1) * tcoeff(i,j)%S
				endif
				A(i,1) = tcoeff(i,j)%W
				A(i,3) = tcoeff(i,j)%E
				A(i,2) = tcoeff(i,j)%P
			end do
			call TDMA(A(:,1:3),T(:,j),A(:,4))
		end do

		do i = 1, psize
			do j = 1, psize
				A(j,4) = b(i,j)
				if (abs(tcoeff(i,j)%E) > 1.0e-7_wp) then
					A(j,4) = A(j,4) - T(i+1,j) * tcoeff(i,j)%E
				endif
				if (abs(tcoeff(i,j)%W) > 1.0e-7_wp) then
					A(j,4) = A(j,4) - T(i-1,j) * tcoeff(i,j)%W
				endif
				A(j,1) = tcoeff(i,j)%S
				A(j,3) = tcoeff(i,j)%N
				A(j,2) = tcoeff(i,j)%P
			end do
			call TDMA(A(:,1:3),T(i,:),A(:,4))
		end do

		residual = resid(tcoeff,T,b)
		if(residual < 1e-6_wp) then
			exit
		end if
	end do
	call cpu_time(timestop)
	write(*,*) "TDMA Finished in ",k-1," iterations(",timestop-timestart,"s) with residual: ", residual

!===========================================
!=--------------Hybrid-1-------------------=
!===========================================


	T = 0.0_wp
	call cpu_time(timestart)
	do k = 1, 100000
		do j = 1, psize
			do i = 1, psize
				A(i,4) = b(i,j)
				if (abs(tcoeff(i,j)%N) > 1.0e-7_wp) then
					A(i,4) = A(i,4) - T(i,j+1) * tcoeff(i,j)%N
				endif
				if (abs(tcoeff(i,j)%S) > 1.0e-7_wp) then
					A(i,4) = A(i,4) - T(i,j-1) * tcoeff(i,j)%S
				endif
				A(i,1) = tcoeff(i,j)%W
				A(i,3) = tcoeff(i,j)%E
				A(i,2) = tcoeff(i,j)%P
			end do
			call TDMA(A(:,1:3),T(:,j),A(:,4))
		end do

		do i = 1, psize
			do j = 1, psize
				A(j,4) = b(i,j)
				if (abs(tcoeff(i,j)%E) > 1.0e-7_wp) then
					A(j,4) = A(j,4) - T(i+1,j) * tcoeff(i,j)%E
				endif
				if (abs(tcoeff(i,j)%W) > 1.0e-7_wp) then
					A(j,4) = A(j,4) - T(i-1,j) * tcoeff(i,j)%W
				endif
				A(j,1) = tcoeff(i,j)%S
				A(j,3) = tcoeff(i,j)%N
				A(j,2) = tcoeff(i,j)%P
			end do
			call TDMA(A(:,1:3),T(i,:),A(:,4))
		end do

		do m = 1, 16
			do i = 1, psize
				do j = 1, psize
					ndrho = 0.0_wp
					if (abs(tcoeff(i,j)%N) > 1.0e-7_wp) then
						ndrho = ndrho + T(i,j+1) * tcoeff(i,j)%N
					endif
					if (abs(tcoeff(i,j)%S) > 1.0e-7_wp) then
						ndrho = ndrho + T(i,j-1) * tcoeff(i,j)%S
					endif
					if (abs(tcoeff(i,j)%W) > 1.0e-7_wp) then
						ndrho = ndrho + T(i-1,j) * tcoeff(i,j)%W
					endif
					if (abs(tcoeff(i,j)%E) > 1.0e-7_wp) then
						ndrho = ndrho + T(i+1,j) * tcoeff(i,j)%E
					endif
					T(i,j) = (b(i,j)-ndrho)/tcoeff(i,j)%P
				end do
			end do
		end do

		residual = resid(tcoeff,T,b)
		if(residual < 1e-6_wp) then
			exit
		end if
	end do
	call cpu_time(timestop)
	write(*,*) "Hybrid 1 Finished in ",k-1," iterations(",timestop-timestart,"s) with residual: ", residual

!===========================================
!=--------------Hybrid-2-------------------=
!===========================================


	T = 0.0_wp
	call cpu_time(timestart)
	do k = 1, 100
		do j = 1, psize
			do i = 1, psize
				A(i,4) = b(i,j)
				if (abs(tcoeff(i,j)%N) > 1.0e-7_wp) then
					A(i,4) = A(i,4) - T(i,j+1) * tcoeff(i,j)%N
				endif
				if (abs(tcoeff(i,j)%S) > 1.0e-7_wp) then
					A(i,4) = A(i,4) - T(i,j-1) * tcoeff(i,j)%S
				endif
				A(i,1) = tcoeff(i,j)%W
				A(i,3) = tcoeff(i,j)%E
				A(i,2) = tcoeff(i,j)%P
			end do
			call TDMA(A(:,1:3),T(:,j),A(:,4))
		end do

		do i = 1, psize
			do j = 1, psize
				A(j,4) = b(i,j)
				if (abs(tcoeff(i,j)%E) > 1.0e-7_wp) then
					A(j,4) = A(j,4) - T(i+1,j) * tcoeff(i,j)%E
				endif
				if (abs(tcoeff(i,j)%W) > 1.0e-7_wp) then
					A(j,4) = A(j,4) - T(i-1,j) * tcoeff(i,j)%W
				endif
				A(j,1) = tcoeff(i,j)%S
				A(j,3) = tcoeff(i,j)%N
				A(j,2) = tcoeff(i,j)%P
			end do
			call TDMA(A(:,1:3),T(i,:),A(:,4))
		end do
		residual = resid(tcoeff,T,b)
		if(residual < 1e-2_wp) then
			exit
		end if
	end do

	do k = 1, 100000
		do m = 1, 16
			do i = 1, psize
				do j = 1, psize
					ndrho = 0.0_wp
					if (abs(tcoeff(i,j)%N) > 1.0e-7_wp) then
						ndrho = ndrho + T(i,j+1) * tcoeff(i,j)%N
					endif
					if (abs(tcoeff(i,j)%S) > 1.0e-7_wp) then
						ndrho = ndrho + T(i,j-1) * tcoeff(i,j)%S
					endif
					if (abs(tcoeff(i,j)%W) > 1.0e-7_wp) then
						ndrho = ndrho + T(i-1,j) * tcoeff(i,j)%W
					endif
					if (abs(tcoeff(i,j)%E) > 1.0e-7_wp) then
						ndrho = ndrho + T(i+1,j) * tcoeff(i,j)%E
					endif
					T(i,j) = (b(i,j)-ndrho)/tcoeff(i,j)%P
				end do
			end do
		end do

		residual = resid(tcoeff,T,b)
		if(residual < 1e-6_wp) then
			exit
		end if
	end do
	call cpu_time(timestop)
	write(*,*) "Hybrid 2 Finished in ",k-1," iterations(",timestop-timestart,"s) with residual: ", residual


	
	call setup(device='svg',fileName='fig-%n.svg',figSize=[400,300])
	call figure()
	call subplot(1,1,1)
	call xylim([0.0_wp,psize*L],[0.0_wp,psize*L])
	call contourf(linspace(0.0_wp,L*psize,psize),linspace(0.0_wp,L*psize,psize),T,10)
	call colorbar(T,10)
	call ticks()
	call labels('x','y','')	
	call show()
contains

subroutine parsebc(bc, pcoeff, b)
	type(boundc),intent(in)::bc
	type(coeff5),intent(inout)::pcoeff
	real(wp),intent(out)::b
	select case(bc%bctype)
		case(0)
!			b = -kval*bc%bcvalue
			b = bc%bcvalue
			pcoeff%P = 1.0_wp
			pcoeff%W = 0.0_wp
			pcoeff%E = 0.0_wp
			pcoeff%S = 0.0_wp
			pcoeff%N = 0.0_wp			
		case(1)
			b = bc%bcvalue
			pcoeff%P = 0.0_wp
	end select
end subroutine

function resid(tcoeff,T,b) result(v)
	type(coeff5),dimension(psize,psize),intent(in)::tcoeff
	real(wp),dimension(psize,psize),intent(in)::T
	real(wp),dimension(psize,psize),intent(in)::b
	real(wp)::rho
	real(wp)::v
	integer::i,j
	v = 0.0_wp

	do i = 1, psize
		do j = 1, psize
			rho = 0.0_wp
			if (abs(tcoeff(i,j)%N) > 1.0e-7_wp) then
				rho = rho + T(i,j+1) * tcoeff(i,j)%N
			endif
			if (abs(tcoeff(i,j)%S) > 1.0e-7_wp) then
				rho = rho + T(i,j-1) * tcoeff(i,j)%S
			endif
			if (abs(tcoeff(i,j)%W) > 1.0e-7_wp) then
				rho = rho + T(i-1,j) * tcoeff(i,j)%W
			endif
			if (abs(tcoeff(i,j)%E) > 1.0e-7_wp) then
				rho = rho + T(i+1,j) * tcoeff(i,j)%E
			endif
			v = v + (rho + T(i,j)*tcoeff(i,j)%P - b(i,j))**2
		end do
	end do
end function
end program main_prg
