program main_prg
	use kinds_mod
	use plplotlib_mod
	use case_mod
	implicit none
	
	real(wp),dimension(:),allocatable::x,t,b
	real(wp),dimension(:,:),allocatable::A
	integer,dimension(:),allocatable::probsize
	real(wp),dimension(:),allocatable::gl2error,error,h
	integer,parameter::bench = 2000
	integer::k
	integer::dist = 3
	
	probsize = floor(linspace(50.0_wp**(1.0_wp/dist),1.0e6_wp**(1.0_wp/dist),bench)**dist)

	allocate(gl2error(bench))
	call setup(device='svg', fileName='plot-%n.svg',figSize=[400,300])
!	call figure()
!	call subplot(1,1,1)

	do k = 1, size(probsize)
		x = linspace(0.0_wp,2.0_wp,probsize(k))
		allocate(A(probsize(k),3))
		allocate(t(probsize(k)))
		allocate(b(probsize(k)))
		allocate(error(probsize(k)))
		call HeatConduction(A,b,x)
		call TDMA(A,t,b)
		gl2error(k) = l2error(x,t)
		error = abs(test(x)-t)/test(x)
!		if(k.eq.size(probsize)) then
!			call printvector(x)
!			write(*,*) ""
!			call printvector(t)
!		endif
		deallocate(A,x,t,b,error)
	end do
	
	allocate(h(size(gl2error)))
	h = 2.0_wp/probsize

	call printvector(h)
	write(*,*) "" 
	call printvector(gl2error)
	call show()
contains
	subroutine HeatConduction(A,b,grid)
		real(wp),dimension(:,:),intent(out)::A
		real(wp),dimension(:),intent(out)::b
		real(wp),dimension(:),intent(in)::grid
		integer::i,psize
		psize = size(grid)
		
		A = 0.0_wp
		do i=2,psize-1
			A(i,2) = -(kmean(grid(i-1),grid(i))/(grid(i)-grid(i-1)) &
					+ kmean(grid(i),grid(i+1))/(grid(i+1)-grid(i)))
			A(i,1) = kmean(grid(i-1),grid(i))/(grid(i)-grid(i-1))
			A(i,3) = kmean(grid(i),grid(i+1))/(grid(i+1)-grid(i))
			b(i) = (source(grid(i))-qval)*((grid(i+1)-grid(i-1)))/2
		end do

		
		select case(tleftbctype)
			case(0)
				A(1,2) = 1
				b(1) = tleftbc
			case(1)
				A(1,2) = kmean(grid(1),grid(2))/(grid(2)-grid(1))
				A(1,3) = kmean(grid(1),grid(2))/(grid(2)-grid(1))
				b(1) = (source(grid(psize)) - qval) * (grid(2)-grid(1))/2 + source(grid(1)) + tleftbc
		end select

		select case(trightbctype)
			case(0)
				A(psize,2) = 1
				b(psize) = trightbc
			case(1)
				A(psize,2) = -kmean(grid(psize-1),grid(psize))/(grid(psize)-grid(psize-1))
				A(psize,1) = kmean(grid(psize-1),grid(psize))/(grid(psize)-grid(psize-1))
				b(psize) = (source(grid(psize)) - qval)*(grid(psize)-grid(psize-1))/2 + source(grid(psize)) - trightbc
		end select
		
	end subroutine

    function kmean(xl, xr) result(v)
        real(wp),intent(in)::xl, xr
        real(wp)::v
        v = 2.0_wp*kval(xl)*kval(xr)/(kval(xl)+kval(xr))
    end function kmean

	subroutine TDMA(A,x,b)
		real(wp),dimension(:,:),intent(inout)::A
		real(wp),dimension(:),intent(in)::b
		real(wp),dimension(:),intent(out)::x
		integer::i,j
		x = b
		do i=1,size(x)-1
			x(i+1) = x(i+1) - x(i)  * A(i+1,1)/A(i,2)
			A(i+1,1:2) = A(i+1,1:2) - A(i,2:3) * A(i+1,1)/A(i,2)
		end do
		do i=size(x),2,-1
			x(i-1) = x(i-1) - x(i) * A(i-1,3)/A(i,2)
			A(i-1,2:3) = A(i-1,2:3) - A(i,1:2) * A(i-1,3)/A(i,2)
		end do
		x(:) = x(:)/A(:,2)
		do i=1,size(x)
			A(i,:) = A(i,:)/A(i,2)
		end do
	end subroutine

	function l2error(x,t) result(ressum)
		real(wp),dimension(:)::x,t
        real(wp)::ressum
        integer::i
		ressum = sqrt(sum((t-test(x))**2)/size(t))
	end function

    function l2resid(A,x,b) result(ressum)
        real(wp),dimension(:,:),intent(in)::A
        real(wp),dimension(:),intent(in)::b,x
        real(wp)::ressum
        integer::i
        ressum = 0
        do i=2,size(X)-1
            ressum = ressum + (A(i,1)*x(i-1)+A(i,2)*x(i)+A(i,3)*x(i+1)-b(i))**2
        end do
        ressum = ressum + (A(1,2)*x(1)+A(1,3)*x(2)-b(1))**2
        i = size(x)
        ressum = ressum + (A(i,1)*x(i-1)+A(i,2)*x(i)-b(i))**2
    end function l2resid

	subroutine printvector(x)
		real(wp),dimension(:),intent(in)::x
		integer::i
		do i = 1, size(x)
			write(*,'(*(E12.4))') x(i)
		end do
	end subroutine printvector

	subroutine printmatrix(A)
		real(wp),dimension(:,:),intent(in)::A
		integer::i
		do i = 1, size(A,1)
			write(*,'(*(E12.4))') A(i,:)
		end do
	end subroutine

	subroutine printexpmatrix(A,b)
		real(wp),dimension(:,:),intent(in)::A
		real(wp),dimension(:),intent(in)::b
		integer::i
		do i = 1, size(A,1)
			write(*,'(*(E12.4))') A(i,:),b(i)
		end do
	end subroutine	
end program main_prg
