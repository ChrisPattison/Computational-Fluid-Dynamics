program main_prg
	use kinds_mod
	use plplotlib_mod
	use case_mod
	use util_mod
	use list_mod
	implicit none
	
	real(wp),dimension(:),allocatable::tupw,tlin,b,fgrid
	type(list)::grid
	real(wp),dimension(:,:),allocatable::A,At
	integer,parameter::psize = 100
	integer::i,k = 0
	logical::allPos

	fgrid = linspace(0.0_wp, 1.0_wp, psize)
	call grid%fill(fgrid)
	deallocate(fgrid)

	call setup(device='svg',fileName='fig-%n.svg',figSize=[1200,900])
	call figure()
	call subplot(1,1,1)
	call xylim([0.0_wp,1.0_wp],[-0.5_wp,1.0_wp])
	
	allPos = .false.
	do while(.not. allPos)
	
		if(allocated(fgrid)) then
			deallocate(A,tupw,tlin,b,fgrid)
		end if		
		k = k + 1

		call grid%freeze(fgrid)
		allocate(tlin(size(fgrid)))
		allocate(tupw(size(fgrid)))
		allocate(b(size(fgrid)))
		allocate(A(size(fgrid),3))
		tlin = 0.0_wp
		tupw = 0.0_wp
		
		call ConvectionDiffusion(fgrid,A,b,1)
		call STDMA(A,tupw,b)
		if(modulo(k,2)==1) then
			call plot(fgrid,tupw, lineColor='b', lineWidth=2.0_wp, markStyle=' ', markSize=0.25_wp)
		endif

		call ConvectionDiffusion(fgrid,A,b,1)
		call STDMA(A,tlin,b)
		if(modulo(k,2)==1) then
			call plot(fgrid,tlin, lineColor='r', lineWidth=1.5_wp, markStyle=' ', markSize=0.25_wp)
		endif
				
		allPos = .true.
		do i = size(fgrid)-1, 1, -1
			if(A(i,3) < 0) then
				allPos = .false.
				call grid%insert(i,(fgrid(i)+fgrid(i+1))/2.0)
			end if
		end do

		write(*,*) k, size(fgrid), grid%length
		if(k==10) then
			exit
		end if
		exit
	end do

	call printvector(abs(tlin-tupw))

	call plot(fgrid,tupw, lineColor='b', lineWidth=2.0_wp, markStyle=' ', markSize=0.25_wp)
	call plot(fgrid,tlin, lineColor='r', lineWidth=1.5_wp, markStyle=' ', markSize=0.25_wp)
	call ticks()
	call labels('x [m]','T [K]','')

	call show()
end program main_prg
