module util_mod
use kinds_mod
implicit none
contains
subroutine GaussJordan(A,x,b)
	real(wp),dimension(:,:),intent(inout)::A
	real(wp),dimension(:),intent(in)::b
	real(wp),dimension(:),intent(out)::x
	integer::i,j
	
	x = b
	if (size(A,1) /= size(A,2) .or. size(A,1) /= size(b) .or. size(b) /= size(x)) then
		return
	end if
	
	do i=1,size(x)
		do j=i+1,size(x)
			if(A(j,i) /= 0) then
				x(j) = x(j) - x(i)*A(j,i)/A(i,i)
				A(j,:) = A(j,:) - A(i,:)*A(j,i)/A(i,i)
			endif
		end do
	end do
	do i=size(x),1,-1
		x(i) = x(i)/A(i,i)
		A(i,:) = A(i,:)/A(i,i)
		do j=i-1,1,-1
			if(A(j,i) /= 0) then
				x(j) = x(j) - x(i)*A(j,i)
				A(j,:) = A(j,:) - A(i,:)*A(j,i)
			endif
		end do
	end do
end subroutine

subroutine STDMA(A,x,b)
	real(wp),dimension(:,:),intent(in)::A
	real(wp),dimension(:),intent(in)::b
	real(wp),dimension(:),intent(out)::x
	real(wp),allocatable,dimension(:,:)::coeff
	
	allocate(coeff(size(A,1),size(A,2)))
	coeff = A
	call TDMA(coeff,x,b)
	deallocate(coeff)
end subroutine

subroutine TDMA(A,x,b)
	real(wp),dimension(:,:),intent(inout)::A !Nx3
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
end subroutine

subroutine TDMVM(A,x,b)
	real(wp),dimension(:,:),intent(in)::A
	real(wp),dimension(:),intent(in)::x
	real(wp),dimension(:),intent(out)::b
	integer::i
	
	i = 1
	b(i) = sum(A(i,2:3)*x(i:i+1))
	i = size(x)
	b(i) = sum(A(i,1:2)*x(i-1:i))
	
	do i = 2, size(x)-1
		b(i) = sum(A(i,:)*x(i-1:i+1))
	end do
end subroutine

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
end module
