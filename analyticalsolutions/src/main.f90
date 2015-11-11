program main_prg
	use kinds_mod
	use plplotlib_mod
	implicit none

	real(wp),dimension(:),allocatable::x,y,yfit,ycubicfit
	real(wp),dimension(:),allocatable::xu,yu
	real(wp),dimension(:,:),allocatable::t
	real(wp),dimension(:,:,:),allocatable::u

	call setup(device='svg',fileName='plot-%n.svg',figSize=[400,300])

	x = linspace(0.0_wp,2.0_wp,101)
	allocate(y(size(x)))
	call heatconduction1d(x,y)

	call figure()
	call subplot(1,1,1)
	call xylim(mixval(x),mixval(y))
	call plot(x,y,lineColor='r',lineWidth=2.0_wp)
	call ticks()
	call labels('x','y','')

	allocate(yfit(size(x)))
	allocate(ycubicfit(size(x)))
	call vortexvel(x,y)
	call vortexfitvel(x,ycubicfit,5,[0.0_wp,0.35_wp,0.8_wp,1.2_wp,1.7_wp,2.0_wp])
	call vortexfitvel(x,yfit,7,[0.0_wp])
	
	call figure()
	call subplot(1,1,1)
	call xylim(mixval(x),mixval(y))
	call plot(x,y,lineColor='r',lineWidth=2.0_wp)
	call plot(x,yfit,lineColor='b',lineWidth=1.25_wp)
	call plot(x,ycubicfit,lineColor='g',lineWidth=1.25_wp)
	call ticks()
	call labels('x','y','')

	x = linspace(-1.0_wp,1.0_wp,101)
	y = linspace(-1.0_wp,1.0_wp,101)
	allocate(t(size(x),size(y)))
	call heatconduction2d(x,y,t)

	call figure()
	call subplot(1,1,1)
	call xylim(mixval(x),mixval(y))
	call contourf(x,y,t,10)
	call colorbar(t,10)
	call ticks()
	call labels('x','y','')


	xu = linspace(-3.0_wp,3.0_wp,41)
	yu = linspace(-3.0_wp,3.0_wp,41)
	allocate(u(3,size(xu),size(yu)))
	call vortex2d(xu,yu,u)
	u(3,:,:) = sqrt(u(1,:,:)**2+u(2,:,:)**2)

	call figure()
	call subplot(1,1,1)
	call xylim(mixval(xu),mixval(yu))
	call quiver(xu,yu,u(1,:,:),u(2,:,:),c=u(3,:,:),scaling=1.5_wp)
	call colorbar(u(3,:,:),10)
	call ticks()
	call labels('x','y','')

	call show()

contains

	subroutine heatconduction1d(x,t)
		real(wp),dimension(:),intent(in)::x
		real(wp),dimension(:),intent(out)::t
		real(wp)::tl=1,tr=3,q=5,k=1,L
		integer::i
		L=maxval(x)-minval(x)
	
		do i=1,size(x)
			t(i) = tl + (L**2/2.0 * q/k + (tr-tl)) * t(i)/L - L**2/2.0 * q/k * (t(i)/L)**2
		end do
	end subroutine heatconduction1d

	function hc2d_l(n,L) result(v)
		real(wp),intent(in)::L
		integer,intent(in)::n
		real(wp)::v
	
		v = (2.0*n + 1)/(2*L)*PI	
	end function hc2d_l

	subroutine heatconduction2d(x,y,t)	
		real(wp),dimension(:),intent(in)::x,y
		real(wp),dimension(:,:),intent(out)::t
		real(wp)::q=5,L
		integer::i,xi,yi
		L=maxval(x)
	
		do xi=1,size(x)
			do yi=1,size(y)
				t(xi,yi) = q*L**2/2*(1-x(xi)**2/L**2) - 2*q*L**2*&
					&sum((/( ((-1.0)**i*cosh(hc2d_l(i,L)*y(yi))*cos(hc2d_l(i,L)*x(xi))&
					&/((hc2d_l(i,L)*L)**3*cosh(hc2d_l(i,L)))),i=0,40)/))
			end do
		end do
	end subroutine heatconduction2d

	subroutine vortex2d(x,y,u)
		real(wp),dimension(:),intent(in)::x,y
		real(wp),dimension(:,:,:),intent(out)::u
		real(wp)::vel=2,trans=1.0,r
	
		integer::xi,yi
		
		do xi=1,size(x)
			do yi=1,size(y)
				r = sqrt(y(yi)**2+x(xi)**2)
				u(1,xi,yi) = cos(atan2(y(yi),x(xi)) - PI/2) * merge(r*vel/trans,trans*vel/r,r.le.trans)
				u(2,xi,yi) = sin(atan2(y(yi),x(xi)) - PI/2) * merge(r*vel/trans,trans*vel/r,r.le.trans)
			end do
		end do
	end subroutine vortex2d

	subroutine vortexfitvel(x,y,order,samplem)
		real(wp),dimension(:),intent(in)::x,samplem
		real(wp),dimension(:),intent(out)::y
		real(wp)::vel=2.0,trans=1.0
		integer,intent(in)::order
		real(wp),allocatable,dimension(:)::coeff,sample
		real(wp),allocatable,dimension(:,:)::A
		integer::i,j
		
		allocate(coeff(order+1))
		allocate(A(order+1,order+1))
		allocate(sample(order+1))

		if (size(samplem).ne.order+1) then
			sample = linspace(minval(x)+minval(x)*0.05_wp,maxval(x),order+1)			
		else
			sample = samplem
		end if
		
		coeff = merge(sample*vel/trans,trans*vel/sample,sample.le.trans)
		do i=1, size(sample)
			A(:,i) = sample**(i-1)
		end do
		
		do i=1,size(sample)
			do j=i+1,size(sample)			
				coeff(j) = coeff(j) - coeff(i)*A(j,i)/A(i,i)
				A(j,:) = A(j,:) - A(i,:)*A(j,i)/A(i,i)
			end do
		end do

		do i=size(sample),1,-1
			coeff(i) = coeff(i)/A(i,i)
			A(i,:) = A(i,:)/A(i,i)
			do j=i-1,1,-1
				coeff(j) = coeff(j) - coeff(i)*A(j,i)
				A(j,:) = A(j,:) - A(i,:)*A(j,i)
			end do
		end do		
		
		y = sum(reshape( (/(x**(i-1)*coeff(i),i=1,size(coeff))/), [size(y),size(coeff)]),2)
	end subroutine vortexfitvel

	subroutine vortexvel(x,y)
		real(wp)::vel=2.0,trans=1.0
		real(wp),dimension(:),intent(in)::x
		real(wp),dimension(:),intent(out)::y
		y = merge(x*vel/trans,trans*vel/x,x.le.trans)
	end subroutine
	
	subroutine printmatrix(A)
		real(wp),dimension(:,:),intent(in)::A
		integer:: i
		do i = 1, size(A,1)
			write(*,'(*(F7.3))') A(i,:)
		end do
	end subroutine printmatrix
end program main_prg
