module case_mod
	use kinds_mod
	implicit none
    real(wp),parameter,public::tleftbc = 1.0_wp
	real(wp),parameter,public::trightbc = exp(2.0_wp)! 0.0_wp
	real(wp),parameter,public::qval = 0.0_wp
	integer,parameter,public::tleftbctype = 0
	integer,parameter,public::trightbctype = 0
contains
	pure elemental function kval(x) result(v)
		real(wp)::v
		real(wp),intent(in)::x
		v = 1.0_wp
	end function

	pure elemental function test(x) result(v)
		real(wp)::v
		real(wp),intent(in)::x
		v = exp(x) 
	end function

	pure elemental function source(x) result(v)
		real(wp)::v
		real(wp),intent(in)::x
		v = kval(x)*exp(x)-qval
	end function
end module
