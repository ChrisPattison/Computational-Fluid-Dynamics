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
		real(wp)::b
	end type

	type::coeff4
		real(wp)::e,w
		real(wp)::n,s
	end type

	type,extends(coeff5)::coeff9
		real(wp)::NN,EE,SS,WW
	end type

	type::boundc
		integer::bctype
		real(wp)::bcvalue
	end type

	type::bounds
		type(boundc)::N
		type(boundc)::E
		type(boundc)::W
		type(boundc)::S
	end type
contains


end module types_mod
