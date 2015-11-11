module types_mod
	use kinds_mod
	implicit none
	type::coeff3
		real(wp)::L,R
		real(wp)::P
	end type

	type::coeff5
		real(wp)::E,W
		real(wp)::N,S
		real(wp)::P
	end type

	type::boundc
		integer::bctype
		real(wp)::bcvalue
	end type
contains


end module types_mod