! Program for 2D Poisson equation with  Jacobi, Gauss Sidel and SOR

program poisson

implicit none



integer nx,ny,i,j, iterations,it, choose

double precision dx,dy,xmax,xmin,ymax,ymin,l2_target,L

double precision pi,b(111,111),p(111,111),pn(111,111),l2norm,iter_diff,denominator



double precision x(111), y(111),p_an(111,111)

double precision error(10000000),omega




pi=-1.0d0
pi=acos(pi)

nx=41
ny=41


omega = 2.0d0/(1 + (pi/nx))

write(*,*) omega
! omega=1.5d0

xmin=0.0d0
xmax=1.0d0

ymin=-0.5d0
ymax=0.5d0

l2_target = 2.0d-7

dx=(xmax-xmin)/(nx-1)
dy=(ymax-ymin)/(ny-1)

! Finally a solution for this shitty 0 problem
x(1)=0.0d0
do i=1,nx-1
x(i+1)=x(i)+dx
enddo

y(1)=-0.5d0
do i=1,ny-1
y(i+1)=y(i)+dy
enddo


L = xmax-xmin

do i=1,nx
do j=1,ny
 b(i,j) = -2*(pi/L)**2*sin(pi*x(i)/L)*cos(pi*y(j)/L)
 enddo
 enddo

do i=1,nx
do j=1,ny
		P(i,j)=0.0d0
		pn(i,j)=0.0d0
enddo
enddo



do i=1,nx
do j=1,ny
 p_an(i,j) = sin(x(i)*pi/L)*cos(y(j)*pi/L)
 enddo
 enddo


it=1
l2norm=1.0d0

do while(l2norm>l2_target)
! do it=1,402
pn=p


	do i=2,nx-1
		do j=2,ny-1
! 
		
		!Jacobi Method
		 p(i,j)= (( dy**2 * ( pn(i+1,j)+pn(i-1,j) ) + dx**2 * ( pn(i,j+1)+pn(i,j-1) ) - (b(i,j)*(dx**2)*(dy**2))) / (2*(dx**2 + dy**2)))



		! !Gauss Seidel Method
		! ! All Tn have become T in the RHS for G-S method and also in the boundary conditions
		 ! p(i,j)= (( dy**2 * ( p(i+1,j)+p(i-1,j) ) + dx**2 * ( p(i,j+1)+p(i,j-1) ) - (b(i,j)*(dx**2)*(dy**2))) / (2*(dx**2 + dy**2)))

		 ! p(i,j)= .25 * (p(i,j-1) + p(i,j+1) + p(i-1,j) + p(i+1,j) - (b(i,j)*(dx**2)))


		 !Successive Over relaxation method

		p(i,j)= (1-omega)*p(i,j) + omega*.25 * (p(i,j-1) + p(i,j+1) + p(i-1,j) + p(i+1,j) - (b(i,j)*(dx**2)))



		enddo
	enddo




! Compute l2norm for each iteration
	l2norm=0.0d0 
	iter_diff=0.0d0
	denominator=0.0d0
	do i=1,nx
	do j=1,ny
	denominator=denominator + (p(i,j)*p(i,j))
	iter_diff= iter_diff + (p(i,j)-pn(i,j))**2
	l2norm=(iter_diff/denominator)**0.5
	enddo
	enddo
	 
	error(it)=l2norm

	iterations= iterations+1
	it=it+1

	write(*,*)  iterations
enddo



 open(unit=21, file="Poisson.txt",action="write",status="replace")

do i=1,nx
    write(21,'(1600F14.7)')(p(i,j),j=1,ny)
enddo


	l2norm=0.0d0 
	iter_diff=0.0d0
	denominator=0.0d0
	do i=1,nx
	do j=1,ny
	denominator=denominator + (p(i,j)**2)
	iter_diff= iter_diff + (p(i,j)-p_an(i,j))**2
	l2norm=(iter_diff/denominator)**0.5
	enddo
	enddo
 

! Save the file for error 
 open(unit=22, file="error_gauss-seidel_41.txt",action="write",status="replace")
do i=1,iterations
write(22,*) error(i),i
enddo
write(*,*) iterations,l2norm,omega


end program poisson





        ! p(:,nx)=0  !p = 0 at y = 2 Correct, top 
        ! p(1,:)=p(2,:) !dp/dx = 0 at x = 0,left
        ! p(:,1)=p(:,2) !dp/dy = 0 at y = 0, bottom
        ! p(nx,:)=p(nx-1,:) !dp/dy = 0 at x = 2, right 