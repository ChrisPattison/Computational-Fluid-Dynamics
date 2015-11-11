module case_mod
	use kinds_mod
	implicit none
	real(wp),parameter,public::kcoef0 = 1.0_wp
	real(wp),parameter,public::kcoef1 = 10.0_wp
    real(wp),parameter,public::tleftbc = 0.0_wp
	real(wp),parameter,public::trightbc = 0.0_wp
	real(wp),parameter,public::qval = 1.0e3_wp
	real(wp),parameter,public::tcoef0 = 70.0_wp
contains
	elemental function kval(T,pos) result(v)
	real(wp)::v
	real(wp),intent(in)::T,pos
	v = kcoef0 + kcoef1 * exp(-T/Tcoef0)
!	v = merge(kcoef0 + kcoef1 * exp(-T/Tcoef0), 10.0_wp, pos.lt.1.0_wp)
	end function
end module
