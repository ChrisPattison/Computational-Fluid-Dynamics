module case_mod
	use kinds_mod
	use util_mod
	implicit none

	type::boundc
		integer::bctype
		real(wp)::bcvalue
	end type

    type(boundc),parameter,public::leftbc = boundc(bctype=0,bcvalue=0.0_wp)
	type(boundc),parameter,public::rightbc = boundc(bctype=0,bcvalue=1.0_wp)
	real(wp),parameter,public::qval = 0.0_wp
	real(wp),parameter,public::kval = 1.0_wp
	real(wp),parameter,public::rhoCpUval =  100.0_wp
!	real(wp),parameter,public::rhoCpUval =  0.000001_wp

contains

	subroutine ConvectionDiffusion(grid,A,b,inttype)
		real(wp),dimension(:),intent(in)::grid
		real(wp),dimension(:,:),intent(out)::A
		real(wp),dimension(:),intent(out)::b
		integer,intent(in)::inttype
		integer::i

		do i=1,size(b)
			A(i,1) = kval/nodeDist(grid,i,i-1)
			A(i,2) = -(merge(0.0_wp, A(i,1), isnan(A(i,1))) + merge(0.0_wp, A(i,3), isnan(A(i,3))))
			A(i,3) = kval/nodeDist(grid,i+1,i)
			call Convection(A(i,:),inttype)
			b(i) = -qval*cellVol(grid,i)
		end do

		A(1,1) = 0.0_wp
		A(size(b),3) = 0.0_wp

		i = 1
		select case(leftbc%bctype)
			case(0)
				A(i,:) = 0.0_wp
				A(i,2) = 1.0_wp
				b(i) = leftbc%bcvalue
			case(1)
				b(i) = b(i) + leftbc%bcvalue
		end select

		i = size(b)
		select case(rightbc%bctype)
			case(0)
				A(i,:) = 0.0_wp
				A(i,2) = 1.0_wp
				b(i) = rightbc%bcvalue
			case(1)
				b(i) = b(i) + rightbc%bcvalue
		end select

	end subroutine

	subroutine Convection(A, inttype)
		real(wp),dimension(:),intent(inout)::A
		integer,intent(in)::inttype
		select case(inttype)
		case(0) !No Convection
		case(1) !Upwinding
			A(1) = A(1) + max(rhoCpUval, 0.0_wp)
			A(2) = A(2) - abs(rhoCpUval)
			A(3) = A(3) + max(-rhoCpUval, 0.0_wp)
		case(2) !Centeral Difference
			A(1) = A(1) + rhoCpUval/2
			A(2) = A(2)
			A(3) = A(3) - rhoCpUval/2
		end select
	end subroutine

	function cellVol(grid,i) result(v)
		real(wp),dimension(:),intent(in)::grid
		integer,intent(in)::i
		real(wp)::v
		
		if(i==1) then
			v = (grid(i+1)-grid(i))/2
		else if (i.eq.size(grid)) then
			v = (grid(i)-grid(i-1))/2
		else
			v = (grid(i+1)-grid(i-1))/2
		end if
	end function

	function nodeDist(grid,i,j) result(v)
		real(wp),dimension(:),intent(in)::grid
		integer,intent(in)::i,j
		real(wp)::v
		
		if(i<1 .or. j<1 .or. i>size(grid) .or. j>size(grid)) then
			v = 0.0_wp
		else
			v = abs(grid(j)-grid(i))
		end if
	end function
end module
