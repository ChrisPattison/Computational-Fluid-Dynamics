program main_prg
	use kinds_mod
	use plplotlib_mod
	use case_mod
	use util_mod
	implicit none
	
	integer,parameter::steps=10.0_wp/h
!	integer,parameter::steps=10.0_wp/h
	real(wp),dimension(3,steps)::x
	real(wp),dimension(3,steps)::xck
	real(wp),dimension(steps)::t
	integer,dimension(2,steps)::ii
	real(wp),parameter::implrelax = 0.3_wp
	real(wp),dimension(3)::xresid
	real(wp),dimension(3,3)::J
	integer,parameter::iimax = 100
	real(wp),parameter::iiresid = 1e-6_wp
	integer::k,i

	write(*,*) 'Time steps:', steps
	x(:,1) = x0
	xck(:,1) = x0
	t(1) = 0
	ii(:,1) = 0
	do k=2,steps

		t(k) = t(k-1) + h 
		x(:,k) = x(:,k-1) + implrelax * dxdt(x(:,k-1), t(k)) * h
		do i=1,iimax
			xresid = x(:,k)
			
			J = dr(x(:,k))
			call GaussJordan(J,x(:,k),-(x(:,k) - x(:,k-1) + dxdt(x(:,k),t(k))))
			x(:,k) = xresid + x(:,k)
			if(sqrt(sum(((xresid-x(:,k))/xresid)**2)) .le. iiresid) then
				exit
			end if
		end do
		if(i.ge.iimax) then
			write(*,*) "Implicit inner iteration (Newton) limit reached"
		end if
		ii(2,k) =  i

!==================================================================


		t(k) = t(k-1) + h 
		x(:,k) = x(:,k-1) + implrelax * dxdt(x(:,k-1), t(k)) * h
		do i=1,iimax
			xresid = x(:,k)
			x(:,k) = x(:,k-1) + implrelax * dxdt(x(:,k),t(k)) * h
			if(sqrt(sum(((xresid-x(:,k))/xresid)**2)) .le. iiresid) then
				exit
			end if
		end do
		if(i.ge.iimax) then
			write(*,*) "Implicit inner iteration (Point) limit reached"
		end if
		ii(1,k) =  i



		x(:,k) = x(:,k-1) + (dxdt(x(:,k),t(k)) + dxdt(x(:,k-1),t(k))) * h / 2.0_wp

!		do i=1,iimax
!			xresid = xck(:,k)
!			xck(:,k) = xck(:,k-1) + implrelax * dxdt(x(:,k),t(k)) * h
!			if(sqrt(sum(((xresid-x(:,k))/xresid)**2)) .le. 1e5) then
!				exit
!			end if
!		end do
!		xck(:,k) = xck(:,k-1) + (dxdt(x(:,k),t(k)) + dxdt(x(:,k-1),t(k))) * h / 2.0_wp
!		if(i.ge.iimax) then
!			write(*,*) "Crank-Nicolson inner iteration limit reached. resid:",sqrt(sum(((xresid-x(:,k))/xresid)**2))
!		end if
	end do

	call setup(device='svg', fileName='plot-%n.svg',figSize=[800,600])

	call figure()
	call subplot(1,1,1)
	call xylim(mixval(t),mixval(x))
	call plot(t,x(1,:),lineColor='r',lineWidth=3.0_wp)
	call plot(t,x(2,:),lineColor='g',lineWidth=3.0_wp)
	call plot(t,x(3,:),lineColor='b',lineWidth=3.0_wp)
	call ticks()
	call labels('t', '', '')

	call figure()
	call subplot(1,1,1)
	call xylim(mixval(t),mixval(real(ii,wp)))
	call plot(t,real(ii(1,:),wp),lineColor='r',lineWidth=3.0_wp)
	call plot(t,real(ii(2,:),wp),lineColor='b',lineWidth=2.0_wp)
	call ticks()
	call labels('t [s]', 'iterations', '')

!	call figure()
!	call subplot(1,1,1)
!	call xyzlim(mixval([x(1,:),xck(1,:)]),mixval([x(2,:),xck(2,:)]),mixval([x(3,:),xck(3,:)]),45.0_wp,45.0_wp)
!	call plot3(x(1,:),x(2,:),x(3,:),lineColor='r')
!	call plot3(xck(1,:),xck(2,:),xck(3,:),lineColor='b')
!	call box('x','y','z')
	
	call show()
contains
end program main_prg
