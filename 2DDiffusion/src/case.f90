module case_mod
	use kinds_mod
	use types_mod
	implicit none
	real(wp)::L = 1.0_wp
	real(wp),public::kval=1.0_wp
	real(wp),public::qval=0.0_wp
	type(boundc)::NBound = boundc(bctype=0,bcvalue=1.0_wp)
	type(boundc)::EBound = boundc(bctype=0,bcvalue=0.0_wp)
	type(boundc)::WBound = boundc(bctype=0,bcvalue=0.0_wp)
	type(boundc)::SBound = boundc(bctype=0,bcvalue=0.0_wp)
contains
end module
