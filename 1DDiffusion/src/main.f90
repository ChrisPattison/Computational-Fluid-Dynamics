program main_prg
	use kinds_mod
	use plplotlib_mod
	implicit none
	
	real(wp),dimension(:),allocatable::x,t,b
	real(wp),dimension(:,:),allocatable::A
	real(wp)::starttime,endtime,GJRT,MGRT,GSRT,OGSRT
	integer,parameter::probsize = 1001
	integer::k,i
	real(wp),dimension(30)::sorval
	integer,dimension(30)::soriter
	
	
	x = linspace(0.0_wp,2.0_wp,probsize)
	allocate(A(probsize,probsize))
	allocate(t(probsize))
	allocate(b(probsize))

	sorval = linspace(1.0_wp,2.0_wp, size(sorval))
	do i=1,size(soriter)
		t=0
		call heatConduction(A,b,x)
		call OptimizedGaussSeidel(a,t,b,sorval(i),1e-7_wp,k)
		soriter(i) = k
		write(*,*) sorval(i),k
	end do
	write(*,*) "============================================="
		
	t=0
	call HeatConduction(A,b,x)
	call cpu_time(starttime)
	call GaussElimination(A,t,b)
	call cpu_time(endtime)
	write(*,*) "Gauss-Jordan done."
	GJRT = endtime-starttime
	write(*,*) "Gauss-Jordan took (seconds):", GJRT
	write(*,*) "============================================="

	t=0
	call HeatConduction(A,b,x)
	call cpu_time(starttime)
	call MultiGrid(x,t,4,k)
	call cpu_time(endtime)
	write(*,*) "Multigrid done with iterations:",k
	MGRT = endtime-starttime
	write(*,*) "MultiGrid took (seconds):", MGRT
	write(*,*) "============================================="
		
	t=0
	call HeatConduction(A,b,x)
	call cpu_time(starttime)
	call GaussSeidel(A,t,b,1.9_wp,1e-7_wp,k)
	call cpu_time(endtime)
	GSRT = endtime-starttime
	write(*,*) "Gauss-Siedel done with iterations:",k
	write(*,*) "Gauss-Siedel took (seconds):", GSRT
	write(*,*) "============================================="
	
	t=0
	call HeatConduction(A,b,x)
	call cpu_time(starttime)
	call OptimizedGaussSeidel(A,t,b,1.9_wp,1e-7_wp,k)
	call cpu_time(endtime)
	OGSRT = endtime-starttime
	write(*,*) "Optimized Gauss-Siedel done with iterations:",k
	write(*,*) "Optimized Gauss-Siedel took (seconds):", OGSRT

	call setup(device='svg', fileName='plot-%n.svg',figSize=[800,600])
	
	call figure()
	call subplot(1,1,1)
	call xylim(mixval(x),mixval(t))
	call plot(x,t,lineColor='r',lineWidth=3.0_wp)
	call ticks()
	call labels('x','T','')
	
	call figure()
	call subplot(1,1,1)
	call xylim(mixval(sorval),mixval(real(soriter,wp)))
	call plot(sorval,real(soriter,wp),lineColor='r',lineWidth=3.0_wp)
	call ticks()
	call labels('Relaxation Factor', 'Iterations','')
	call show()
contains
	subroutine HeatConduction(A,b,grid)
		real(wp),dimension(:,:),intent(out)::A
		real(wp),dimension(:),intent(out)::b
		real(wp),dimension(:),intent(in)::grid
		real(wp),parameter::q=1.0_wp,tl=0.0_wp,tr=1.0_wp,k=1.0_wp
		integer::i,psize
		psize = size(grid)
		
		A = 0
		do i=2,psize-1
			A(i,i) = -(grid(i+1)-grid(i-1))
			A(i,i-1) = grid(i)-grid(i-1)
			A(i,i+1) = grid(i+1)-grid(i)
			A(i,:) = A(i,:)/((grid(i)-grid(i-1))*(grid(i+1)-grid(i)))
			b(i) = -q*(grid(i+1)-grid(i-1))/k
		end do
		A(1,1) = 1
		b(1) = tl
		A(psize,psize) = 1
		b(psize) = tr
	end subroutine
	
	subroutine MultiGrid(grid,T,levels,totaliter)
		integer,intent(in)::levels
		real(wp),dimension(:),intent(in)::grid
		real(wp),dimension(:),intent(out)::T
		integer,intent(out)::totaliter
		real(wp),allocatable,dimension(:,:)::A
		real(wp),allocatable,dimension(:)::b
		real(wp),allocatable,dimension(:)::subgrid,subT
		real(wp),allocatable,dimension(:)::supergrid,superT
		integer,allocatable,dimension(:)::leveldof
		integer::i,k
		totaliter = 0
		
		allocate(leveldof(levels))
		leveldof = floor(linspace(sqrt(real(size(grid),wp)),real(size(grid),wp),levels))
		!write(*,*) "Multigrid levels:"
		!write(*,*) leveldof
		
		allocate(superT(leveldof(1)))
		allocate(supergrid(leveldof(1)))
		supergrid = linspace(minval(grid),maxval(grid),leveldof(1))
		superT = 0
		
		do i=1,levels
			!write(*,*) "---------------------------------------------"
			!write(*,*) "Grid Level:",i
			allocate(A(leveldof(i),leveldof(i)))
			allocate(b(leveldof(i)))
			allocate(subT(leveldof(i)))
			allocate(subgrid(leveldof(i)))
			if(i.eq.levels) then
				subgrid = grid
			else
				subgrid = linspace(minval(grid),maxval(grid),leveldof(i))
			endif
			call map(supergrid,subgrid,superT,subT)
			deallocate(supergrid, superT)

			call HeatConduction(A,b,subgrid)
			call OptimizedGaussSeidel(A,subT,b,1.7_wp,1e-7_wp,k)
			!write(*,*) "Done with iterations:",k
			totaliter = totaliter + k
			deallocate(A,b)
			if(i.eq.levels) then
				exit
			endif
			
			allocate(supergrid(leveldof(i)))
			allocate(superT(leveldof(i)))
			supergrid = subgrid !TODO: Alias instead of copy
			superT = subT
			deallocate(subgrid, subT)
		end do
		T = subT
	end subroutine
	
	subroutine GaussElimination(A,x,b)
		real(wp),dimension(:,:),intent(inout)::A
		real(wp),dimension(:),intent(in)::b
		real(wp),dimension(:),intent(out)::x
		integer::i,j
		
		x = b
		if (size(A,1).ne.size(A,2).or.size(A,1).ne.size(b).or.size(b).ne.size(x)) then
			return
		end if
		
		do i=1,size(x)
			do j=i+1,size(x)
				if(A(j,i).ne.0) then
					x(j) = x(j) - x(i)*A(j,i)/A(i,i)
					A(j,:) = A(j,:) - A(i,:)*A(j,i)/A(i,i)
				endif
			end do
		end do
		do i=size(x),1,-1
			x(i) = x(i)/A(i,i)
			A(i,:) = A(i,:)/A(i,i)
			do j=i-1,1,-1
				if(A(j,i).ne.0) then
					x(j) = x(j) - x(i)*A(j,i)
					A(j,:) = A(j,:) - A(i,:)*A(j,i)
				endif
			end do
		end do
	end subroutine

	subroutine GaussSeidel(A,x,b,w,r,totaliter) !A must be tridiagonal
		real(wp),dimension(:,:),intent(inout)::A
		real(wp),dimension(:),intent(in)::b
		real(wp),intent(in)::r
		real(wp),dimension(:),intent(inout)::x
		real(wp),intent(in)::w
		integer,intent(out)::totaliter
		integer::i,iter
		integer,parameter::maxiter=2e6
		real(wp)::resid
		resid = 1e20
		totaliter = 0
		if (size(A,1).ne.size(A,2).or.size(A,1).ne.size(b).or.size(b).ne.size(x)) then
			write(*,*) "GS Matrix Size Mismatch"
			totaliter = -1
			return
		endif		
		!x = (x(1)+x(size(x)))/2
		
		do while(resid.ge.r.and.totaliter.lt.maxiter)
			do iter=1,250
				x(1)=(b(1)-A(1,2)*x(2))/A(1,1)
				do i=2,size(x)-1
					x(i)=(1-w)*x(i) + w*(b(i)-(A(i,i-1)*x(i-1)+A(i,i+1)*x(i+1)))/A(i,i)
				end do
				x(size(x))=(b(size(x))-A(size(x),size(x)-1)*x(size(x)-1))/A(1,1)
				totaliter = totaliter + 1
			end do
			resid = sum((matmul(A,x)-b)**2)
		end do
	end subroutine

	subroutine OptimizedGaussSeidel(A,x,b,w,r,totaliter) !A must be tridiagonal
		real(wp),dimension(:,:),intent(inout)::A
		real(wp),dimension(:),intent(in)::b
		real(wp),intent(in)::r
		real(wp),dimension(:),intent(inout)::x
		real(wp),intent(in)::w
		integer,parameter::maxiter=2e6
		integer,intent(out)::totaliter
		integer::i,iter
		real(wp)::resid
		resid = 1e20
		totaliter = 0
		if (size(A,1).ne.size(A,2).or.size(A,1).ne.size(b).or.size(b).ne.size(x)) then
			write(*,*) "GS Matrix Size Mismatch"
			totaliter = -1
			return
		endif		
		!x = (x(1)+x(size(x)))/2
		
		do while(resid.ge.r.and.totaliter.lt.maxiter)
			call invertdiagonal(A)
			do iter=1,250
				x(1)=(b(1)-A(1,2)*x(2))*A(1,1)
				do i=2,size(x)-1
					x(i)=(1-w)*x(i) + w*(b(i)-(A(i,i-1)*x(i-1)+A(i,i+1)*x(i+1)))*A(i,i)
				end do
				x(size(x))=(b(size(x))-A(size(x),size(x)-1)*x(size(x)-1))*A(1,1)
				totaliter = totaliter + 1
			end do
			call invertdiagonal(A)			
			resid = sum((matmul(A,x)-b)**2)
		end do
	end subroutine

	subroutine invertdiagonal(A)
		real(wp),dimension(:,:),intent(inout)::A
		integer::i
		do i=1,size(A,1)
			A(i,i) = 1/A(i,i)
		end do
	end subroutine

	subroutine printvector(x)
		real(wp),dimension(:),intent(in)::x
		integer::i
		do i = 1, size(x)
			write(*,'(*(F7.3))') x(i)
		end do
	end subroutine printvector

	subroutine printmatrix(A)
		real(wp),dimension(:,:),intent(in)::A
		integer::i
		do i = 1, size(A,1)
			write(*,'(*(F7.3))') A(i,:)
		end do
	end subroutine

	subroutine map(grid,subgrid,values,subgridvalues) !TODO: Fix O(N) version
		real(wp),dimension(:),intent(in)::grid,values,subgrid
		real(wp),dimension(:),intent(out)::subgridvalues
		integer::i,j=1
		
		do i=1,size(subgrid)
			do j=2,size(grid)
				if(grid(j).ge.subgrid(i).and.grid(j-1).le.subgrid(i)) then
					subgridvalues(i) = values(j-1)+(values(j)-values(j-1))*(subgrid(i)-grid(j-1))/(grid(j)-grid(j-1))
					exit
				endif
			end do
		end do
	end subroutine

!	subroutine map(grid,subgrid,values,subgridvalues)
!		real(wp),dimension(:),intent(in)::grid,values,subgrid
!		real(wp),dimension(:),intent(out)::subgridvalues
!		integer::i,j=1
!		
!		do i=2,size(grid)
!			do while(subgrid(j)<grid(i))
!				subgridvalues(j) = values(i-1)+(values(i)-values(i-1))*(subgrid(j)-grid(i-1))/(grid(i)-grid(i-1))
!				j = j+1
!			end do
!		end do
!		subgridvalues(j:) = values(size(values))
!		
!		call printmatrix(reshape([grid,values],[size(grid),2]))
!		write(*,*) "Mapped to"
!		call printmatrix(reshape([subgrid,subgridvalues],[size(subgrid),2]))
!	end subroutine
	
end program main_prg
