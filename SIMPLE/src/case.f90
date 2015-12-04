module case_mod
	use iso_fortran_env
	use kinds_mod
	use types_mod
	use coeffutil_mod
	implicit none
	integer,parameter::psize = 160
	real(wp)::L = 1.0_wp/psize
	real(wp),public::muval=0.005_wp
!	real(wp),public::muval=0.002_wp
	real(wp),public::rhoval=1.0_wp

	type(bounds)::PBounds = bounds( &
		boundc(bctype=1,bcvalue=1.0_wp), &
		boundc(bctype=1,bcvalue=2.0_wp), &
		boundc(bctype=1,bcvalue=3.0_wp), &
		boundc(bctype=1,bcvalue=4.0_wp))

	type(bounds)::UBounds = bounds( &
		boundc(bctype=0,bcvalue=1.0_wp), &
		boundc(bctype=0,bcvalue=0.0_wp), &
		boundc(bctype=0,bcvalue=0.0_wp), &
		boundc(bctype=0,bcvalue=0.0_wp))

	type(bounds)::VBounds = bounds( &
		boundc(bctype=0,bcvalue=0.0_wp), &
		boundc(bctype=0,bcvalue=0.0_wp), &
		boundc(bctype=0,bcvalue=0.0_wp), &
		boundc(bctype=0,bcvalue=0.0_wp))
contains

function cellvol(coeff,i,j) result(v)
	type(coeff5),dimension(:,:),intent(inout)::coeff
	integer,intent(in)::i,j
	real(wp)::v

	v = L*L
	if(i == size(coeff,1) .or. i == 1) then
		v = v/2
	end if
	if(j == size(coeff,2) .or. j == 1) then
		v = v/2
	end if
end function

subroutine parsebc(bc, coeff)
	type(boundc),intent(in)::bc
	type(coeff5),intent(inout)::coeff
	select case(bc%bctype)
		case(0)
			coeff = coeff5(0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,bc%bcvalue)
		case(1)
			coeff = coeff5(0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp)
			if(bc%bcvalue < 1.5_wp) then !N
				coeff%S = 1.0_wp
			elseif(bc%bcvalue < 2.5_wp) then !E
				coeff%W = 1.0_wp
			elseif(bc%bcvalue < 3.5_wp) then !W
				coeff%E = 1.0_wp
			else !S
				coeff%N = 1.0_wp
			endif
	end select
end subroutine

subroutine pressurebc(coeff)
	type(coeff5),dimension(:,:),intent(inout)::coeff
	integer::i,j
	i = size(coeff,1)
	j = size(coeff,2)
	coeff(1,:) = coeff5(0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp)
	coeff(i,:) = coeff5(0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp)
	coeff(:,1) = coeff5(0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp)
	coeff(:,j) = coeff5(0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp)

	coeff(1,:)%N = coeff(1,:)%N + 1.0_wp
	coeff(i,:)%S = coeff(i,:)%S + 1.0_wp
	coeff(:,1)%E = coeff(:,1)%E + 1.0_wp
	coeff(:,j)%W = coeff(:,j)%W + 1.0_wp

	coeff(1,1) = coeff5(0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp) ! fix node
end subroutine

subroutine momentum(coeff, u)
	type(coeff5),dimension(:,:),intent(inout)::coeff
	real(wp),dimension(:,:,:),intent(in)::u
	type(coeff4)::uface,Pe
	logical,parameter::upwind = .true.
	integer::i,j
	do i = 1,size(coeff,1)
		do j = 1, size(coeff,2)
			uface = vface(u,i,j)
			Pe%n = uface%n*L/muval
			Pe%s = uface%s*L/muval
			Pe%e = uface%e*L/muval
			Pe%w = uface%w*L/muval
			if(upwind) then
				coeff(i,j)%N = muval + max(0.0_wp,-Pe%n)
				coeff(i,j)%S = muval + max(0.0_wp, Pe%s)
				coeff(i,j)%E = muval + max(0.0_wp,-Pe%e)
				coeff(i,j)%W = muval + max(0.0_wp, Pe%w)
			else
				coeff(i,j)%N = muval*max(0.0_wp,(1.0_wp-abs(Pe%n)/10.0_wp)**5) + max(0.0_wp,-Pe%n)
				coeff(i,j)%S = muval*max(0.0_wp,(1.0_wp-abs(Pe%s)/10.0_wp)**5) + max(0.0_wp, Pe%s)
				coeff(i,j)%E = muval*max(0.0_wp,(1.0_wp-abs(Pe%e)/10.0_wp)**5) + max(0.0_wp,-Pe%e)
				coeff(i,j)%W = muval*max(0.0_wp,(1.0_wp-abs(Pe%w)/10.0_wp)**5) + max(0.0_wp, Pe%w)
			endif
			coeff(i,j)%P = coeff(i,j)%N + coeff(i,j)%S + coeff(i,j)%E + coeff(i,j)%W + &
				L*(uface%e - uface%w + uface%n - uface%s)
		end do
	end do
end subroutine

function vface(field, i, j) result(v)
	real(wp),dimension(:,:,:),intent(in)::field
	integer,intent(in)::i,j
	type(coeff4)::v

	if(i == size(field,1)) then
		v%e = field(i,j,1)
	else
		v%e = (field(i+1,j,1)+field(i,j,1))/2.0_wp
	endif
	if(i == 1) then
		v%w = field(i,j,1)
	else
		v%w = (field(i,j,1) + field(i-1,j,1))/2.0_wp
	endif

	if(j == size(field,2)) then
		v%n = field(i,j,2)
	else
		v%n = (field(i,j+1,2) + field(i,j,2))/2.0_wp
	endif
	if(j == 1) then
		v%s = field(i,j,2)
	else
		v%s = (field(i,j,2) + field(i,j-1,2))/2.0_wp
	endif
end function

function face(field, i, j) result(v)
	real(wp),dimension(:,:),intent(in)::field
	integer,intent(in)::i,j
	type(coeff4)::v

	if(i == size(field,1)) then
		v%e = field(i,j)
	else
		v%e = (field(i+1,j)+field(i,j))/2.0_wp
	endif
	if(i == 1) then
		v%w = field(i,j)
	else
		v%w = (field(i,j) + field(i-1,j))/2.0_wp
	endif

	if(j == size(field,2)) then
		v%n = field(i,j)
	else
		v%n = (field(i,j+1) + field(i,j))/2.0_wp
	endif
	if(j == 1) then
		v%s = field(i,j)
	else
		v%s = (field(i,j) + field(i,j-1))/2.0_wp
	endif
end function

function cellgrad(field, i, j) result(v)
	real(wp),dimension(:,:),intent(in)::field
	integer,intent(in)::i,j
	type(coeff4)::v

	if(i == size(field,1)) then
		v%e = 0.0_wp
	else
		v%e = (field(i+1,j)-field(i,j))/L
	endif
	if(i == 1) then
		v%w = 0.0_wp
	else
		v%w = (field(i,j) - field(i-1,j))/L
	endif

	if(j == size(field,2)) then
		v%n = 0.0_wp
	else
		v%n = (field(i,j+1)-field(i,j))/L
	endif
	if(j == 1) then
		v%s = 0.0_wp
	else
		v%s = (field(i,j) - field(i,j-1))/L
	endif
end function

subroutine umomentum(coeff, u, p, pgrad)
	type(coeff5),dimension(:,:),intent(inout)::coeff
	real(wp),dimension(:,:,:),intent(in)::u
	real(wp),dimension(:,:),intent(in)::p
	real(wp),dimension(:,:,:),intent(out)::pgrad
	integer::i,j
	call momentum(coeff, u)
	call gradient(p,pgrad)
	do i = 1, size(coeff,1)
		do j = 1, size(coeff,2)
			coeff(i,j)%b = -cellvol(coeff,i,j)*pgrad(i,j,1)
		end do
	end do
end subroutine

subroutine vmomentum(coeff, u, p, pgrad)
	type(coeff5),dimension(:,:),intent(inout)::coeff
	real(wp),dimension(:,:,:),intent(in)::u
	real(wp),dimension(:,:),intent(in)::p
	real(wp),dimension(:,:,:),intent(out)::pgrad
	integer::i,j
	call momentum(coeff, u)
	call gradient(p,pgrad)
	do i = 2, size(coeff,1)-1
		do j = 2, size(coeff,2)-1
			coeff(i,j)%b = -cellvol(coeff,i,j)*pgrad(i,j,2)
		end do
	end do
end subroutine

subroutine pcorrection(coeff, u, p, d, pgrad)
	type(coeff5),dimension(:,:),intent(inout)::coeff
	real(wp),dimension(:,:,:),intent(in)::u
	real(wp),dimension(:,:),intent(in)::p
	type(coeff4)::F,pg,dface,rcpg,uface
	real(wp),dimension(:,:,:)::pgrad
	real(wp),dimension(:,:)::d
	integer::i,j

	call momentum(coeff, u)
	d = -L**2/coeff%P
	call gradient(p, pgrad)

	do i = 1, size(coeff,1)
		do j = 1, size(coeff,2)
			dface = face(d,i,j)
			pg = cellgrad(p,i,j)
			rcpg = vface(pgrad,i,j)
			uface = vface(u,i,j)
			F%e = L*(uface%e + dface%e * (pg%e - rcpg%e))
			F%w = L*(uface%w + dface%w * (pg%w - rcpg%w))
			F%n = L*(uface%n + dface%n * (pg%n - rcpg%n))
			F%s = L*(uface%s + dface%s * (pg%s - rcpg%s))
			coeff(i,j)%N = dface%n
			coeff(i,j)%S = dface%s
			coeff(i,j)%E = dface%e
			coeff(i,j)%W = dface%w
			coeff(i,j)%P = coeff(i,j)%N + coeff(i,j)%S + coeff(i,j)%E + coeff(i,j)%W
			coeff(i,j)%b = (F%e - F%w + F%n - F%s)
		end do
	end do
end subroutine

subroutine ucorrection(temp, u, ucorr, pcorr)
	real(wp),dimension(:,:,:),intent(inout)::u,ucorr
	real(wp),dimension(:,:),intent(in)::pcorr
	type(coeff5),dimension(:,:),intent(out)::temp
	integer::i,j

	call momentum(temp, u)

	do i = 2, size(pcorr,1)-1
		do j = 2, size(pcorr,2)-1
			ucorr(i,j,1) = - L**2/temp(i,j)%P * (pcorr(i+1,j)-pcorr(i-1,j)) / (2.0_wp*L)
			ucorr(i,j,2) = - L**2/temp(i,j)%P * (pcorr(i,j+1)-pcorr(i,j-1)) / (2.0_wp*L)
		end do
	end do
end subroutine

subroutine gradient(field, grad)
	real(wp),dimension(:,:),intent(in)::field
	real(wp),dimension(:,:,:),intent(out)::grad
	integer::i,j

	do i = 2, size(field,1)-1
		do j = 2, size(field,2)-1
			grad(i,j,1) = (value(field,i+1,j) - value(field,i-1,j))/(L*2.0_wp)
			grad(i,j,2) = (value(field,i,j+1) - value(field,i,j-1))/(L*2.0_wp)
		end do
	end do

	grad(1,:,1) = (field(2,:)-field(1,:))/L
	grad(size(field,1),:,1) = (field(size(field,1),:)-field(size(field,1)-1,:))/L
	grad(:,1,2) = (field(:,2)-field(:,1))/L
	grad(:,size(field,2),2) = (field(:,size(field,1))-field(:,size(field,1)-1))/L
end subroutine

function continuity(field) result(v)
	real(wp),dimension(:,:,:),intent(in)::field
	real(wp)::v
	type(coeff4)::u
	integer::i,j

	v = 0.0_wp
	do i = 1, size(field,1)
		do j =1, size(field,2)
			u = vface(field,i,j)
			v = v + (u%n - u%s + u%e - u%w)
		end do
	end do
end function

subroutine boundary(coeff,bound)
	type(coeff5),dimension(:,:),intent(inout)::coeff
	type(bounds)::bound
	integer::i,j

	do j = 1, size(coeff,1)
		do i = 1, size(coeff,2)
			if(i==size(coeff,1)) then
				call parsebc(bound%E, coeff(i,j))
			else if(i==1) then
				call parsebc(bound%W, coeff(i,j))
			end if

			if(j==size(coeff,2)) then
				call parsebc(bound%N, coeff(i,j))
			else if(j==1) then
				call parsebc(bound%S, coeff(i,j))
			end if
		end do
	end do
end subroutine
end module
