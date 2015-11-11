program main_prg
	use kinds_mod
	use plplotlib_mod
	use case_mod
	implicit none
	
	real(wp),dimension(:),allocatable::x,t,b,conductivity
	real(wp),dimension(:,:),allocatable::A
	integer,parameter::probsize = 101
	integer::k
	
	x = linspace(0.0_wp,2.0_wp,probsize)
	allocate(A(probsize,3))
	allocate(t(probsize))
	allocate(b(probsize))
	allocate(conductivity(probsize))
	
	t=(tleftbc+trightbc)/2
	k=0
	call HeatConduction(A,t,b,x)
    do while(l2resid(A,t,b).gt.1e-6)
        call TDMA(A,t,b)
        call HeatConduction(A,t,b,x)
        k = k + 1
        write(*,*) k,l2resid(A,t,b)
        if(k.ge.5e2) then
        	exit
        end if
    end do
	write(*,*) 'Picard finished with iterations:',k
		
	call setup(device='svg', fileName='plot-%n.svg',figSize=[800,600])
	
  	call figure()
	call subplot(1,1,1)
	call xylim(mixval(x),mixval(t))
	call plot(x,t,lineColor='r',lineWidth=3.0_wp)
	call ticks()
	call labels('x [m]','Temperature T [k]','')
	
	call figure()
	call subplot(2,1,1)
	conductivity = kval(T,x)
	call xylim(mixval(x),mixval(conductivity))
	call plot(x,conductivity,lineColor='r',lineWidth=3.0_wp)
	call ticks()
	call labels('x [m]', 'Conductivity k [W/m-K]', '')

	T = linspace(0.0_wp, 50.0_wp,probsize)
 	x = 0
	call subplot(2,1,2)
	conductivity = kval(t,x)
	call xylim(mixval(t),mixval(conductivity))
	call plot(t,conductivity,lineColor='r',lineWidth=3.0_wp)
	call ticks()
	call labels('T [K]', 'Conductivity k [W/m-K]', '')
	call show()
contains
	subroutine HeatConduction(A,t,b,grid)
		real(wp),dimension(:,:),intent(out)::A
		real(wp),dimension(:),intent(out)::b
		real(wp),dimension(:),intent(in)::grid
        real(wp),dimension(:),intent(in)::t
		integer::i,psize
		psize = size(grid)
		
		A = 0.0_wp
		do i=2,psize-1
			A(i,2) = -(kmean(t(i-1),t(i),grid(i-1),grid(i))/(grid(i)-grid(i-1)) &
					+ kmean(t(i),t(i+1),grid(i),grid(i+1))/(grid(i+1)-grid(i)))
			A(i,1) = kmean(t(i-1),t(i),grid(i-1),grid(i))/(grid(i)-grid(i-1))
			A(i,3) = kmean(t(i),t(i+1),grid(i),grid(i+1))/(grid(i+1)-grid(i))
			b(i) = -qval*(grid(i+1)-grid(i-1))
		end do
		A(1,2) = 1
		b(1) = tleftbc
		A(psize,2) = 1
		b(psize) = trightbc
	end subroutine

!    subroutine HeatJacobi(J,T,grid)
!        real(wp),dimension(:,:),intent(out)::J
!        real(wp),dimension(:),intent(in)::T,grid
!        real(wp)::dxe,dxw,dkedt,dkwdt
!        integer::i
!
!        do i = 2, size(T)-1
!            dxw = grid(i)-grid(i-1)
!            dxe = grid(i+1)-grid(i)
!
!            dkwdt = dk(t(i-1),t(i))
!            J(i,1) = dkwdt*t(i-1)/dxw &
!               + kmean(t(i-1),t(i)) &
!               - dkwdt * T(i) / dxw
!
!            dkwdt = dk(t(i),t(i-1))
!            dkedt = dk(t(i),t(i+1))
!            J(i,2) = dkwdt * t(i-1)/dxw &
!               - (dkwdt/dxw+dkedt/dxe)*T(i) &
!               - (kmean(t(i),t(i-1))/dxw + kmean(t(i+1),t(i))/dxe) &
!               + dkedt * t(i+1)/dxe
!
!            dkedt = dk(t(i+1),t(i))
!            J(i,3) = dkedt*t(i+1)/dxe &
!               + kmean(t(i+1),t(i)) &
!               - dkedt * T(i) / dxe
!        end do
!
!        i=1
!        J(i,1) = 0
!        dkedt = dk(t(i),t(i+1))
!        dxe = grid(i+1)-grid(i)
!        J(i,2) = dkedt * t(i-1)/dxe &
!           - dkedt/dxe * T(i) &
!           - kmean(t(i),t(i-1))/dxe
!        dkedt = dk(t(i+1),t(i))
!        J(i,3) = dkedt*t(i+1)/dxe &
!           + kmean(t(i+1),t(i)) &
!           - dkedt * T(i) / dxe
!
!        i=size(t)
!        dxw = grid(i)-grid(i-1)
!        dkwdt = dk(t(i-1),t(i))
!        J(i,1) = dkwdt*t(i-1)/dxw &
!           + kmean(t(i-1),t(i)) &
!           - dkwdt * T(i) / dxw
!
!        dkwdt = dk(t(i),t(i-1))
!        J(i,2) = dkwdt * t(i-1)/dxw &
!           - dkwdt/dxw * T(i) &
!           - kmean(t(i),t(i-1))/dxw
!        J(i,3) = 0
!    end subroutine
!
!    function dk(t1,t2) result(v)!evaluates dk/dT1 at T1
!        real(wp),intent(in)::t1,t2
!        real(wp)::v
!        real(wp)::k1,k2
!        k1 = kval(t1)
!        k2 = kval(t2)
!        v = 2.0_wp * (1.0_wp/k1 + 1.0_wp/k2)**(-2) * k1**(-2) * kcoef1/tcoef0*exp(-t1/tcoef0)
!    end function dk

!    function kval(T) result(v)
!        real(wp)::v
!        real(wp),intent(in)::T
!        v = kcoef0 + kcoef1 * exp(-T/Tcoef0)
!    end function kval

    function kmean(Tl,Tr, xl, xr) result(v)
        real(wp),intent(in)::Tl,Tr, xl, xr
        real(wp)::v
        v = 2.0_wp/(1.0_wp/kval(Tl,xl)+1.0_wp/kval(Tr,xr))
!		v = (kval(Tl)  + kval(Tr))/2.0_wp
    end function kmean

	subroutine TDMA(A,x,b)
		real(wp),dimension(:,:),intent(inout)::A
		real(wp),dimension(:),intent(in)::b
		real(wp),dimension(:),intent(out)::x
		integer::i,j
		x = b
		do i=1,size(x)-1
			x(i+1) = x(i+1) - x(i)  * A(i+1,1)/A(i,2)
			A(i+1,1:2) = A(i+1,1:2) - A(i,2:3) * A(i+1,1)/A(i,2)
		end do
		do i=size(x),2,-1
			x(i-1) = x(i-1) - x(i) * A(i-1,3)/A(i,2)
			A(i-1,2:3) = A(i-1,2:3) - A(i,1:2) * A(i-1,3)/A(i,2)
		end do
		x(:) = x(:)/A(:,2)
		do i=1,size(x)
			A(i,:) = A(i,:)/A(i,2)
		end do
	end subroutine

    function l2resid(A,x,b) result(ressum)
        real(wp),dimension(:,:),intent(in)::A
        real(wp),dimension(:),intent(in)::b,x
        real(wp)::ressum
        integer::i
        ressum = 0
        do i=2,size(X)-1
            ressum = ressum + (A(i,1)*x(i-1)+A(i,2)*x(i)+A(i,3)*x(i+1)-b(i))**2
        end do
        ressum = ressum + (A(1,2)*x(1)+A(1,3)*x(2)-b(1))**2
        i = size(x)
        ressum = ressum + (A(i,1)*x(i-1)+A(i,2)*x(i)-b(i))**2
    end function l2resid

!    subroutine resid(A,x,b,R)
!        real(wp),dimension(:,:),intent(in)::A
!        real(wp),dimension(:),intent(in)::b,x
!        real(wp),dimension(:),intent(out)::R
!        integer::i
!        do i=2,size(X)-1
!            R(i) =  (A(i,1)*x(i-1)+A(i,2)*x(i)+A(i,3)*x(i+1)-b(i))
!        end do
!        R(1) (A(1,2)*x(1)+A(1,3)*x(2)-b(1))**2
!        i = size(x)
!        R(i) = (A(i,1)*x(i-1)+A(i,2)*x(i)-b(i))**2
!    end subroutine resid

	subroutine printvector(x)
		real(wp),dimension(:),intent(in)::x
		integer::i
		do i = 1, size(x)
			write(*,'(*(E12.4))') x(i)
		end do
	end subroutine printvector

	subroutine printmatrix(A)
		real(wp),dimension(:,:),intent(in)::A
		integer::i
		do i = 1, size(A,1)
			write(*,'(*(E12.4))') A(i,:)
		end do
	end subroutine

	subroutine printexpmatrix(A,b)
		real(wp),dimension(:,:),intent(in)::A
		real(wp),dimension(:),intent(in)::b
		integer::i
		do i = 1, size(A,1)
			write(*,'(*(E12.4))') A(i,:),b(i)
		end do
	end subroutine	
end program main_prg
