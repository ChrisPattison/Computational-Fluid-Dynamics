module coeffutil_mod
use kinds_mod
use types_mod
implicit none
contains

subroutine gaussseidel(coeff,R,a)
	type(coeff5),dimension(:,:),intent(in)::coeff
	real(wp),dimension(:,:),intent(inout)::R
	real(wp),intent(in)::a
	integer::i,j
	real(wp)::ndrho
	
	do i = 1, size(coeff,1)
		do j = 1, size(coeff,2)
			ndrho = value(R,i,j+1)*coeff(i,j)%N + value(R,i,j-1)*coeff(i,j)%S + &
				value(R,i-1,j)*coeff(i,j)%W + value(R,i+1,j) *coeff(i,j)%E
			R(i,j) = (1.0_wp-a)*R(i,j) + a*(coeff(i,j)%b+ndrho)/coeff(i,j)%P
		end do
	end do
end subroutine

subroutine gsiter(coeff, R, a, iter)
	type(coeff5),dimension(:,:),intent(in)::coeff
	real(wp),dimension(:,:),intent(inout)::R
	real(wp),intent(in)::a
	integer,intent(in)::iter
	integer::k

	do k = 1, iter
		call gaussseidel(coeff,R, a)
	end do
end subroutine

function resid(coeff,R) result(v)
	type(coeff5),dimension(:,:),intent(in)::coeff
	real(wp),dimension(:,:),intent(in)::R
	real(wp)::rho
	real(wp)::v
	integer::i,j
	v = 0.0_wp

	do i = 1, size(coeff,1)
		do j = 1, size(coeff,2)
			rho = value(R,i,j+1)*coeff(i,j)%N + value(R,i,j-1)*coeff(i,j)%S + &
				value(R,i-1,j)*coeff(i,j)%W + value(R,i+1,j)*coeff(i,j)%E
			v = v + (rho - R(i,j)*coeff(i,j)%P - coeff(i,j)%b)**2
		end do
	end do
end function

subroutine residvec(coeff,R,res)
	type(coeff5),dimension(:,:),intent(in)::coeff
	real(wp),dimension(:,:),intent(in)::R
	real(wp)::rho
	real(wp),dimension(:,:),intent(out)::res
	integer::i,j

	do i = 1, size(R,1)
		do j = 1, size(R,2)
			rho = value(R,i,j+1)*coeff(i,j)%N + value(R,i,j-1)*coeff(i,j)%S + &
				value(R,i-1,j)*coeff(i,j)%W + value(R,i+1,j)*coeff(i,j)%E
			res(i,j) = rho - R(i,j)*coeff(i,j)%P -coeff(i,j)%b
		end do
	end do
end subroutine

subroutine solve(coeff, R, a, kmax, minres)
	type(coeff5),dimension(:,:),intent(inout)::coeff
	real(wp),dimension(:,:),intent(inout)::R
	real(wp),intent(in)::a
    integer,intent(in)::kmax
	real(wp),intent(in)::minres
	integer::k
	
	do k = 1, kmax
		call gaussseidel(coeff,R,a)
		
		if(modulo(k,5) == 1) then
			if(resid(coeff,R) < minres) then
				exit
			end if
		end if
	end do
end subroutine

function value(R,i,j) result(v)
	real(wp),dimension(:,:),intent(in)::R
	integer,intent(in)::i,j
	real(wp)::v

	if(i <= ubound(R,1).and.i >= lbound(R,1) .and. j <= ubound(R,2) .and. j >= lbound(R,2)) then
		v = R(i,j)
	else
		v = 0.0_wp
	end if
end function

subroutine zero(coeff)
	type(coeff5),dimension(:,:),intent(inout)::coeff
	integer::i,j
	
	do i = 1, size(coeff,1)
		do j = 1, size(coeff,2)
			coeff(i,j)%N = 0.0_wp
			coeff(i,j)%S = 0.0_wp
			coeff(i,j)%E = 0.0_wp
			coeff(i,j)%W = 0.0_wp
			coeff(i,j)%P = 0.0_wp
			coeff(i,j)%b = 0.0_wp
		end do
	end do          
end subroutine
end module coeffutil_mod