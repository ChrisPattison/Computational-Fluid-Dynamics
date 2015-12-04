program main_prg
	use kinds_mod
	use case_mod
	use util_mod
	use types_mod
	implicit none

	real(wp),allocatable,dimension(:,:,:)::U,UCorr,Uold
	real(wp),allocatable,dimension(:,:)::P,PCorr
	real(wp),allocatable,dimension(:,:)::temps
	real(wp),allocatable,dimension(:,:,:)::tempv
	type(coeff5),allocatable,dimension(:,:)::coeff
	logical,parameter::printsol = .true.
	integer::k
	integer,parameter::sol_unit=output_unit
	real(wp)::pres,ures,cont,starttime,stoptime
!	real(wp),parameter::amoment=0.7_wp, apcorr=0.3_wp
!	integer,parameter::momentumiter=2,pcorriter=5
	real(wp),parameter::amoment=0.7_wp, apcorr=0.1_wp
	integer,parameter::momentumiter=1,pcorriter=5

	call cpu_time(starttime)

	allocate(U(psize,psize,2))
	allocate(Uold(psize,psize,2))
	allocate(UCorr(psize,psize,2))
	allocate(P(psize,psize))
	allocate(PCorr(psize,psize))
	allocate(coeff(psize,psize))
	allocate(temps(psize,psize))
	allocate(tempv(psize,psize,2))

	P = 0.0_wp
	U = 0.0_wp
	PCorr = 0.0_wp
	call zero(coeff)

	write(error_unit,*) "Iter   |  P  Resid  |  U  Resid  |  Continuity"

	do k = 1,50000
		Uold = U
		call umomentum(coeff,U,P,tempv)
		call boundary(coeff, UBounds)
		call gsiter(coeff, U(:,:,1), amoment, momentumiter)

		call vmomentum(coeff,U,P,tempv)
		call boundary(coeff, VBounds)
		call gsiter(coeff, U(:,:,2), amoment, momentumiter)

		call pcorrection(coeff,U,p,temps,tempv)
		call pressurebc(coeff)
		PCorr = 0.0_wp
		call gsiter(coeff, PCorr, 1.0_wp, pcorriter)
		PCorr = apcorr*PCorr
		
		call ucorrection(coeff, U, UCorr, PCorr)
		P = P + PCorr
		U = U + UCorr
		pres = norm2(PCorr)
		ures = norm2(Uold-U)
		cont = continuity(u)
		if(pres>1.0e98_wp.or.ures>1.0e98_wp.or.cont>1.0e98_wp.or.pres/=pres.or.ures/=ures.or.cont/=cont) then
			exit
		endif
		if((pres<1e-6_wp.and.ures<1e-6_wp).and.k>1) then
			exit
		endif
		write(error_unit,'(I7,3(E13.4))') k, pres, ures, merge(cont,0.0_wp,abs(cont)>1.0e-40_wp)
	end do
	call cpu_time(stoptime)
	write(error_unit,'(A,I5,A,F6.1,A)') "Done in ",k," iterations and ", stoptime-starttime,"s"
	if(printsol) then
		call writesolution(sol_unit,U,P)
	endif
	!close(sol_unit)
contains
subroutine writesolution(unit,U,P)
	integer,intent(in)::unit
	real(wp),intent(in),dimension(:,:,:)::U
	real(wp),intent(in),dimension(:,:)::P
	integer::i

	write(unit,'(*(E12.4))') [(real(i-1,wp)/real(size(U,1)-1,wp),i=1,size(U,1))]
	write(unit,'(*(E12.4))') [(real(i-1,wp)/real(size(U,2)-1,wp),i=1,size(U,2))]

	write(unit,*) size(U,2)
	do i = 1, size(U,2)
		write(unit,'(*(E12.4))') min(max(merge(U(:,i,1), 0.0_wp, abs(U(:,i,1))>1e-20_wp),-1.0e98_wp),1.0e98_wp)
	end do

	write(unit,*) size(U,2)
	do i = 1, size(U,2)
		write(unit,'(*(E12.4))') min(max(merge(U(:,i,2), 0.0_wp, abs(U(:,i,2))>1e-20_wp),-1.0e98_wp),1.0e98_wp)
	end do

	write(unit,*) size(P,2)
	do i = 1, size(P,2)
		write(unit,'(*(E12.4))') min(max(merge(P(:,i), 0.0_wp, abs(P(:,i))>1e-20_wp),-1.0e98_wp),1.0e98_wp)
	end do
end subroutine
end program main_prg
