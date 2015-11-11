program main_prg
	use kinds_mod
	use plplotlib_mod
	use case_mod
	use util_mod
	implicit none

	real(wp),dimension(:),allocatable::dt
	integer::steps
!	integer,parameter::steps=int(200/dt)
	real(wp),dimension(:),allocatable::grid,b,time,ImplT
	real(wp),dimension(:,:),allocatable::A,ExplA,Ac,T
!	integer,parameter::points = int(5e2)
	integer,parameter::points = 50
	integer,parameter::bench = 20
	real(wp),dimension(bench,2)::error
	real(wp)::cnw = 1.5_wp
	real(wp),dimension(:),allocatable::terror
	integer::i,j,k
	
	dt = 10**linspace(log10(1e-6_wp),log10(1e-4_wp),bench)
!	dt = 10**linspace(log10(5e-4_wp),log10(5e-5_wp),bench)
	allocate(A(points,3))
	allocate(ExplA(points,3))
	allocate(Ac(points,3))
	allocate(b(points))
	allocate(ImplT(points))
	grid = linspace(0.0_wp,2.0_wp*PI,points)
!	grid = linspace(0.0_wp,1.0_wp,points)
	
	write(*,*) ' ', ' dt ', ' Impl Euler Error ', ' Crank-Nicolson Error '
	
	do j = bench, 1, -1
!		steps=int(dt(bench)/dt(j))
		steps=int(1.0_wp/dt(j))
!		steps=int(2.0_wp*PI/dt(j))
		allocate(T(points,steps))
		allocate(time(steps))
	
		T(:,1) = solution(grid,0.0_wp)
		time(1) = 0.0_wp
		call CoeffSetup(A,grid,dt(j))
		call ExplCoeffSetup(ExplA,grid,dt(j))
		
!===================================================		
!---------------Explicit Euler----------------------
!===================================================		
!		do k = 2, steps
!			time(k) = time(k-1) + dt(j)
!			call ManExplStep(b,grid,time(k-1),dt(j))
!			call TDMVM(A,T(:,k-1),T(:,k))
!			T(:,k) = (/((T(i,k) + b(i))*dt(j)/rhoCp/cellVol(grid,i),i=2,size(b)-1)/)
!		end do
!		error(j,2) = sqrt(sum((solution(grid,maxval(time)) - T(:,size(T,2)))**2)/size(T,1))

!===================================================		
!---------------Crank-Nicolson----------------------
!===================================================
		do k = 2, steps
			time(k) = time(k-1) + dt(j)

			call ManStep(b,T(:,k-1),grid,time(k),dt(j))
			Ac = A
			call TDMA(Ac,ImplT,b)
			
			call ManExplStep(b,grid,time(k-1),dt(j))
			call TDMVM(ExplA,T(:,k-1),T(:,k))
			
			T(:,k) = (/((T(i,k) + b(i))*dt(j)/rhoCp/cellVol(grid,i),i=1,size(b))/)
			T(:,k) = ((2-cnw)*T(:,k) + cnw*ImplT)/2
		end do
		error(j,2) = sqrt(sum((solution(grid,maxval(time)) - T(:,size(T,2)))**2)/size(T))
!		error(j,2) = sqrt(sum((/(solution(grid,time(i)) - T(:,i), i=1,size(time))/)**2)/size(T,2))

!===================================================		
!---------------Implicit Euler----------------------
!===================================================		
		
		do k = 2, steps
			time(k) = time(k-1) + dt(j)
			call ManStep(b,T(:,k-1),grid,time(k),dt(j))
			Ac = A			
			call TDMA(Ac,T(:,k),b)
		end do
		error(j,1) = sqrt(sum((solution(grid,maxval(time)) - T(:,size(T,2)))**2)/size(T))
!		error(j,1) = sqrt(sum((/(solution(grid,time(i)) - T(:,i), i=1,size(time))/)**2)/size(T,2))



!===================================================		
	
	
		if(j.ne.1) then
			deallocate(T,time)
		end if
		
		write(*,*) j, dt(j), "||", error(j,1), "|", error(j,2)!, "|", error(j,3)
	end do

!	call printvector(T(:,size(T,2)))

	call setup(device='pngqt', fileName='plot-%n.png',figSize=[800,600])
!	call setup(device='svg', fileName='plot-%n.svg',figSize=[800,600])

	call figure()
	call subplot(1,1,1)



!================
!-----ErrorVt---
!================

	allocate(terror(size(T,1)))
	terror = (/(sqrt(sum((solution(grid,time(j)) - T(:,j))**2)/size(T(:,1:j))), j=1,size(T,2))/)
	call xylim(mixval(time),log10(mixval(terror)))
	call plot(time,log10(terror),lineColor='r',lineWidth=3.0_wp)
	call ticks(logy=.true.)
	call labels('t [s]', 'Error [K]', '')

!================
!-----ErrorVdt---
!================

!	dt = log10(dt)
!	error = log10(error)
!	call xylim(mixval(dt),mixval(error))
!	call plot(dt,error(:,1),lineColor='r',lineWidth=3.0_wp)
!	call plot(dt,error(:,2),lineColor='b',lineWidth=2.5_wp)
!	call ticks(logx=.true.,logy=.true.)
!	call labels('dt [s]', 'Error [K]', '')

!================
!-----Solution---
!================

!	call xylim(mixval(grid),[0.0_wp,1.0_wp])
!	call plot(grid,T(:,size(T,2)),lineColor='r',lineWidth=3.0_wp)
!	call plot(grid,T(:,size(T,2)/2),lineColor='o',lineWidth=3.0_wp)
!	call plot(grid,T(:,1),lineColor='g',lineWidth=3.0_wp)
!	call plot(grid,solution(grid,time(size(time))),lineColor='b',lineWidth=2.0_wp)
!	call ticks()
!	call labels('x [m]', 'T [K]', '')

!================
!-----ErrorVx----
!================

!	call xylim(mixval(grid),mixval(abs(solution(grid,maxval(time)) - T(:,size(T,2)))))
!	call plot(grid,abs(solution(grid,maxval(time)) - T(:,size(T,2))),lineColor='r',lineWidth=3.0_wp)
!	call ticks()
!	call labels('x [m]', '[K]', '')

	call show()
contains
end program main_prg
