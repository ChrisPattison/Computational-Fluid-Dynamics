module case_mod
	use kinds_mod
	implicit none
	real(wp),dimension(3),public::x0 = [-10,-10,30]
	real(wp),public::sigma = 10
	real(wp),public::rho = 28
	real(wp),public::beta = 8.0_wp/3.0_wp
	real(wp),parameter::h=1.0e-3_wp
contains
function dxdt(x,t) result(dx)
	real(wp),dimension(3),intent(in)::x
	real(wp),intent(in)::t
	real(wp),dimension(3)::dx
	
	dx(1) = sigma * (x(2)-x(1))
	dx(2) = x(1) * (rho-x(3))-x(2)
	dx(3) = x(1) * x(2) - beta * x(3)
end function

function dr(x) result(J)
	real(wp),dimension(3),intent(in)::x
	real(wp),dimension(3,3)::J
	
	J(1,1) = 1 -sigma
	J(1,2) = sigma
	J(1,3) = 0
	
	J(2,1) = rho-x(3)
	J(2,2) = 1 - 1
	J(2,3) = -x(1)
	
	J(3,1) = x(2)
	J(3,2) = x(1)
	J(3,3) = 1 -beta
	
	J = J 
end function
end module
