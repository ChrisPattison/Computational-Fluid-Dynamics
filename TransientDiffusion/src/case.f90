module case_mod
	use kinds_mod
	implicit none
	real(wp),public::kval=1.0_wp
	real(wp),public::qval=0.0_wp
	real(wp),public::rhoCp = 1.0_wp
	integer,public::lbcType = 1
	integer,public::rbcType = 1
	real(wp),public::lbc = 0.0_wp
	real(wp),public::rbc = 0.0_wp
contains
subroutine CoeffSetup(A,grid,dt)
	real(wp),dimension(:,:),intent(out)::A
	real(wp),dimension(:),intent(in)::grid
	real(wp),intent(in)::dt
	integer::i
	
	
	do i = 2, size(grid)-1
		A(i,1) = kval/(grid(i)-grid(i-1))
		A(i,2) = -((grid(i+1)-grid(i-1))/2 * rhoCp / dt + kval/(grid(i+1)-grid(i)) + kval/(grid(i)-grid(i-1)))
		A(i,3) = kval/(grid(i+1)-grid(i))
	end do

	i = 1
	A(i,:) = 0.0_wp
	select case(lbcType)
		case(0)
			A(i,2) = 1.0_wp		
		case(3)
			A(i,2) = 1.0_wp		
		case(1)
			A(i,2) = -((grid(i+1)-grid(i))/2 * rhoCp / dt + kval/(grid(i+1)-grid(i)))
			A(i,3) = kval/(grid(i+1)-grid(i))
	end select
	
	i = size(grid)
	A(i,:) = 0.0_wp
	select case(rbcType)
		case(0)
			A(i,2) = 1.0_wp
		case(3)
			A(i,2) = 1.0_wp
		case(1)
			A(i,1) = kval/(grid(i)-grid(i-1))
			A(i,2) = -((grid(i)-grid(i-1))/2 * rhoCp / dt + kval/(grid(i)-grid(i-1)))
	end select
end subroutine

subroutine ExplCoeffSetup(A,grid,dt)
	real(wp),dimension(:,:),intent(out)::A
	real(wp),dimension(:),intent(in)::grid
	real(wp),intent(in)::dt
	integer::i
	
	
	do i = 2, size(grid)-1
		A(i,1) = kval/(grid(i)-grid(i-1))
		A(i,2) = -(-(grid(i+1)-grid(i-1))/2 * rhoCp / dt + kval/(grid(i+1)-grid(i)) + kval/(grid(i)-grid(i-1)))
		A(i,3) = kval/(grid(i+1)-grid(i))
	end do

	i = 1
	A(i,:) = 0.0_wp
	select case(lbcType)
		case(0)
			A(i,2) = 1.0_wp		
		case(3)
			A(i,2) = 1.0_wp		
		case(1)
			A(i,2) = -(-(grid(i+1)-grid(i))/2 * rhoCp / dt + kval/(grid(i+1)-grid(i)))
			A(i,3) = kval/(grid(i+1)-grid(i))
	end select
	
	i = size(grid)
	A(i,:) = 0.0_wp
	select case(rbcType)
		case(0)
			A(i,2) = 1.0_wp
		case(3)
			A(i,2) = 1.0_wp
		case(1)
			A(i,1) = kval/(grid(i)-grid(i-1))
			A(i,2) = -(-(grid(i)-grid(i-1))/2 * rhoCp / dt + kval/(grid(i)-grid(i-1)))
	end select
end subroutine

subroutine Step(b,t,grid,dt)
	real(wp),dimension(:),intent(out)::b
	real(wp),dimension(:),intent(in)::t,grid
	real(wp),intent(in)::dt
	integer::i
	
	do i = 2, size(grid)-1
		b(i) = -(qval + t(i)*rhoCp/dt)*(grid(i+1)-grid(i-1))/2
	end do
	
	i = 1
	select case(lbcType)
		case(0)
			b(i) = lbc
		case(1)
			b(i) = -(qval + t(i)*rhoCp/dt)*(grid(i+1)-grid(i))/2 - lbc
	end select
	
	i=size(grid)
	select case(rbcType)
		case(0)
			b(i) = rbc
		case(1)
			b(i) = -(qval + t(i)*rhoCp/dt)*(grid(i)-grid(i-1))/2 - rbc
	end select
end subroutine 

subroutine ManStep(b,t,grid,time,dt)
	real(wp),dimension(:),intent(out)::b
	real(wp),dimension(:),intent(in)::t,grid
	real(wp),intent(in)::time,dt
	integer::i
	
	do i = 2, size(grid)-1
		b(i) = -(((source(grid(i-1),time+dt) + source(grid(i),time+dt)*2 + source(grid(i+1),time+dt))/4 &
			+ qval) * dt +t(i)*rhoCp/dt)*(grid(i+1)-grid(i-1))/2
	end do
	
	i = 1
	select case(lbcType)
		case(0)
			b(i) = lbc
		case(1)
			b(i) = -(((source(grid(i),time+dt)*3 + source(grid(i+1),time+dt))/4 &
				+ qval) * dt + t(i)*rhoCp/dt)*(grid(i+1)-grid(i))/2 - lbc
		case(3)
			b(i) = solution(grid(i),time+dt)
	end select
	
	i=size(grid)
	select case(rbcType)
		case(0)
			b(i) = rbc
		case(1)
			b(i) = -(((source(grid(i),time+dt)*3 + source(grid(i-1),time+dt))/4 &
				+ qval) * dt + t(i)*rhoCp/dt)*(grid(i)-grid(i-1))/2 - rbc
		case(3)
			b(i) = solution(grid(i),time+dt)
	end select
end subroutine

subroutine ManExplStep(b,grid,time,dt)
	real(wp),dimension(:),intent(out)::b,grid
	real(wp),intent(in)::time,dt
	integer::i
	
	do i = 2, size(grid)-1
		b(i) = -(((source(grid(i-1),time+dt) + source(grid(i),time+dt)*2 + source(grid(i+1),time+dt))/4 &
			+ qval) * dt)*(grid(i+1)-grid(i-1))/2
	end do
	
	b(1) = 0
	b(size(b)) = 0	
!	i = 1
!	select case(lbcType)
!		case(0)
!			b(i) = lbc
!		case(1)
!			b(i) = -(((source(grid(i),time+dt)*3 + source(grid(i+1),time+dt))/4 &
!				+ qval) * dt + t(i)*rhoCp/dt)*(grid(i+1)-grid(i))/2 - lbc
!	end select
!	
!	i=size(grid)
!	select case(rbcType)
!		case(0)
!			b(i) = rbc
!		case(1)
!			b(i) = -(((source(grid(i),time+dt)*3 + source(grid(i-1),time+dt))/4 &
!				+ qval) * dt + t(i)*rhoCp/dt)*(grid(i)-grid(i-1))/2 - rbc
!	end select
end subroutine

elemental function source(x,t) result(v)
	real(wp),intent(in)::x,t
	real(wp)::v
!	v = (kval*cos(t)-sin(t))*sin(x)-qval
	v = -(x-x**2)*sin(t)-2*kval*cos(t)-qval
!	v = cos(t) - qval
!	v = 0
!	v = exp(t)
end function

elemental function solution(x,t) result(v)
	real(wp),intent(in)::x,t
	real(wp)::v
!	v = sin(x)*cos(t)
	v = (x-x**2)*cos(t)
!	v = sin(t)
!	v = 0
!	v = exp(t)
end function

function cellVol(grid,i) result(v)
	real(wp),dimension(:),intent(in)::grid
	integer,intent(in)::i
	real(wp)::v
	
	if(i.eq.1) then
		v = (grid(i+1)-grid(i))/2
	else if (i.eq.size(grid)) then
		v = (grid(i)-grid(i-1))/2
	else
		v = (grid(i+1)-grid(i-1))/2
	end if
end function
end module
