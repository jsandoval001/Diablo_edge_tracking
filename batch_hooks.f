C******************************************************************************|
C these are the hooks required by the 'basic' version of diablo, even
C when not running the 'ensemble' flavor    
C
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine init_march
	include 'header_batch'
	real version, current_version
	real*8 zbqlu01, dummy
	integer i,j,k,v
	
C	Initialization for Random Number Generator
	call zbqlini(0)
	dummy = 0.2
	
C	Initialization for Measurements
C	ADD THESE TO INPUT.DAT !!!
	meas_type_adjoint(1) = 'P'
	meas_var(1) = 0.001
	meas_type_adjoint(2) = 'SH1'
	meas_var(2) = 0.001
	meas_type_adjoint(3) = 'SH3'
	meas_var(3) = 0.001

C	Intialize other variables needed for measurements
	ave_var = 0.0 ! Gaussian Averaging Over Measurements - ADD TO INPUT.DAT
	smpl_freq = 100 ![Hz]						 ADD TO INPUT.DAT
	smpl_freq = 1/(delta_t*smpl_freq)
	if (smpl_freq.eq.0) then
	 smpl_freq = 1
	end if
	write(6,*) 'Measurement Freq [Hz]: ',1./(smpl_freq*delta_t)
      write(6,*) 'Measurement Freq [Index]: ', smpl_freq	
				
	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine init_adjoint
	include 'header_batch'
	real version, current_version
	integer i,j,k	
	
c	do j=0,nxm
c	  do k=0,tnkz
c	    do i=0,nkx
c	      cu1(i,k,j) = 0.0;
c	      cu1(i,k,j) = 0.0;
c		cu1(i,k,j) = 0.0;
c	    end do
c	  end do
c	end do
		
c	open(13,file='input_adjoint.dat',form='formatted',status='old')  
	
c	current_version = 1.0
c	read(13,*)

	
c	close(13)
	
C	write(6,*) 'Real-Time Scaling: ',scaling
	
	return
	end 

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine quasi_rev_per(TEMP5,I,J,K)
c If direction = -1 then add to RHS of CN coefficients.
c If direction = +1 then add to implicit solve of 
	include 'header_batch'
	integer J,K,I
	real*8 TEMP5

	if (verbosity.gt.2) then
	  write(*,*) 'This is Quasi Reversible Periodic !'
	end if
	
	TEMP5 = TEMP5 + reg_mixed*(KX2(I)+KY2(J)+KZ2(K))
	
	return
	end



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine save_stats_adjoint(final)
c	include 'header_batch'
	logical final
	integer i,j,k
	
c	if (num_per_dir.eq.3) then
c	  write(6,*) 'Saving stats for ensemble...'
	  
c	  open(unit=23,file="field.dat",status="unknown")
c	  write(23,*)
C	  close(23)
	  
	  
c	  call system('mv field.dat ../post_process/matlab/field.dat',i)
	  
c	else
c	  write(6,*) 'Non-periodic cases not supported'
c	end if
	
	
	return
	end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine quasi_reversible_xz
	include 'header_batch'
	integer J,K,I
	real*8 small 
	real*8 TEMP1, KZ4(0:TNKZ), KX4(0:NKX)
	parameter (small = 1.0e-20)

	if (verbosity.gt.2) then
	  write(*,*) 'This is Quasi Reversible XZ-Plane !'
	end if
	
c	RK Constant from Main Program
      TEMP1 = -NU * H_BAR(RK_STEP) / 2.0
		
	DO I=0,NKX
	  KX4(I) = KX2(I)*KX2(I)
	END DO
	DO K=0,TNKZ
	  KZ4(K) = KZ2(K)*KZ2(K)
	END DO
	
C Add the explicit part of the hyperviscosity terms
C Compute the horizontal viscous terms in fourier space
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CR1(I,K,J) 
     &          + small*TEMP1*(KX4(I)+KZ4(K))*CU1(I,K,J)
            CR3(I,K,J)=CR3(I,K,J)
     &          + small*TEMP1*(KX4(I)+KZ4(K))*CU3(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CR2(I,K,J)=CR2(I,K,J)
     &          + small*TEMP1*(KX4(I)+KZ4(K))*CU2(I,K,J)
          END DO
        END DO
      END DO
	
	return
	end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine quasi_reversible_y
	include 'header_batch'
	integer J,K,I
	real*8 TEMP1,   alpha_m3, alpha_m2, alpha_m1
	real*8 alpha_0, alpha_1,  alpha_2,  alpha_3
	real*8 small 
	parameter (small = 1.0e-20)
		
	if (verbosity.gt.2) then
	  write(*,*) 'This is Quasi Reversible Y-Direction !'
	end if
	
c	RK Constant from Main Program
      TEMP1 = - NU * H_BAR(RK_STEP) / 2.0

C Add the explicit part of the hyperviscosity terms
C Compute the horizontal viscous terms in fourier space
	if (NY.gt.5) then

C Lower Boudary Finite Difference Method: 4th-order Forward Finite
c Difference Approximation
	J=2
	  alpha_3 =  24./(dy(J+3)+dy(J+2)+dy(J+1)+dy(J))/
     &               (dy(J+2)+dy(J+1)+dy(J))/(dy(J+1)+dy(J))/dy(J)
	  alpha_2 = -24./(dy(J+3)+dy(J+2)+dy(J+1))/
     &                               (dy(J+2)+dy(J+1))/dy(J+1)/dy(J)
     	  alpha_1 =  24./(dy(J+1)+dy(J))/dy(J+1)/
     &                                      dy(J+2)/(dy(J+3)+dy(J+2))
	  alpha_0 = -24./dy(J+2)/(2.*dy(j+2)*dy(j+1)+(dy(j+1)**2.)
     &          +(dy(j+2)**2)+dy(j)*dy(j+1)+dy(j)*dy(j+2))/dy(j+3)
	  alpha_m1 = 24./dy(j+3)/(3.*dy(j+3)*(dy(j+2)**2.)
     &              +dy(j+3)*(dy(j+1)**2.)+dy(j)*(dy(j+3)**2.)
     &           +4.*dy(j+3)*dy(j+2)*dy(j+1)+2.*(dy(j+3)**2.)
     &              *dy(j+1)+(dy(j+2)**2.)*dy(j)+2.*dy(j+3)*dy(j+2)
     &              *dy(j)+dy(j+2)*(dy(j+1)**2.)+dy(j+3)*dy(j+1)
     &              *dy(j)+(dy(j+3)**3.)+dy(j+2)*dy(j+1)*dy(j)
     &          +3.*(dy(j+3)**2.)*dy(j+2)+2.*(dy(j+2)**2.)*dy(j+1)
     &             +(dy(j+2)**3.))   



	  	  
        DO K=0,TNKZ
          DO I=0,NKX
		CR3(I,K,j)=CR3(I,K,j) + small*TEMP1*
     &                ( alpha_m1*CU3(I,K,J-1) + 
     &                   alpha_0*CU3(I,K,J) + 
     &                   alpha_1*CU3(I,K,J+1) + 
     &                   alpha_2*CU3(I,K,J+2) + 
     &                   alpha_3*CU3(I,K,J+3)   )

	    END DO
        END DO

c	  alpha_3 = 0.0
c	  alpha_2 = 0.0
c        alpha_1 = 0.0
c        alpha_0 = 0.0
c        alpha_m1 = 0.0
		  	  
        DO K=0,TNKZ
          DO I=0,NKX
		CR1(I,K,j)=CR1(I,K,j) + small*TEMP1*
     &                ( alpha_m1*CU1(I,K,J-1) + 
     &                   alpha_0*CU1(I,K,J) + 
     &                   alpha_1*CU1(I,K,J+1) + 
     &                   alpha_2*CU1(I,K,J+2) + 
     &                   alpha_3*CU1(I,K,J+3)   )
	    end do
	  end do

	
      DO J=3,NYM-1
	
C 4th Order Central Finite Difference Approximation fo the 4th Derivative
C	  alpha_2  ----> f(x(j+2))
C	  alpha_1  ----> f(x(j+1))
C	  alpha_0  ----> f(x(j  ))
C	  alpha_m1 ----> f(x(j-1))
C	  alpha_m2 ----> f(x(j-2))

	  alpha_2 =  24./(dy(j-1)+dy(j+1)+dy(j+2)+dy(j))/
     &           (dy(j+2)+dy(j+1)+dy(j))/(dy(j+2)+dy(j+1))/dy(j+2)
	  alpha_1 = -24./(dy(j-1)+dy(j)+dy(j+1))/
     &                               (dy(j)+dy(j+1))/dy(j+1)/dy(j+2)
	  alpha_0 =  24./(dy(j-1)+dy(j))/dy(j)/
     &                                      dy(j+1)/(dy(j+2)+dy(j+1))
	  alpha_m1 = -24./dy(j-1)/dy(j)/((dy(j)**2.)+(dy(j+1)**2.)
     &            +dy(j+2)*dy(j)+2.*dy(j+1)*dy(j)+dy(j+2)*dy(j+1))
	  alpha_m2 =  24./dy(j-1)/(3.*(dy(j-1)**2.)*dy(j)+(dy(j-1)**3.)
     &          +(dy(j)**3.)+dy(j+2)*dy(j+1)*dy(j-1)+(dy(j+1)**2.)
     &          *dy(j)+dy(j+2)*dy(j+1)*dy(j)+(dy(j+1)**2.)*dy(j-1)
     &           +dy(j+2)*(dy(j)**2.)+2.*dy(j+1)*(dy(j-1)**2.)
     &        +4.*dy(j+1)*dy(j-1)*dy(j)+2.*dy(j+2)*dy(j-1)*dy(j)
     &        +2.*dy(j+1)*(dy(j)**2.)+dy(j+2)*(dy(j-1)**2.)
     &        +3.*dy(j-1)*(dy(j)**2.))



        DO K=0,TNKZ
          DO I=0,NKX     
		CR3(I,K,J)=CR3(I,K,J) + small*TEMP1*
     &                ( alpha_m2*CU3(I,K,J-2) + 
     &                  alpha_m1*CU3(I,K,J-1) + 
     &                  alpha_0*CU3(I,K,J)   + 
     &                  alpha_1*CU3(I,K,J+1) + 
     &                  alpha_2*CU3(I,K,J+2)   )
          end do
	  end do

c	  alpha_2 = 0.0
c        alpha_1 = 0.0
c        alpha_0 = 0.0
c        alpha_m1 = 0.0
c        alpha_m2 = 0.0

        DO K=0,TNKZ
          DO I=0,NKX 
		CR1(I,K,J)=CR1(I,K,J) + small*TEMP1*
     &                ( alpha_m2*CU1(I,K,J-2) + 
     &                  alpha_m1*CU1(I,K,J-1) + 
     &                  alpha_0*CU1(I,K,J)   + 
     &                  alpha_1*CU1(I,K,J+1) + 
     &                  alpha_2*CU1(I,K,J+2)   )
     
	    end do
	  end do
	  
      END DO

C Upper Boudary Finite Difference Method: 4th-order Backwards Finite
c Difference Approximation
	J = NYM
	  alpha_1  =  24./(dy(j)+dy(j-2)+dy(j+1)+dy(j-1))/
     &             (dy(j)+dy(j+1)+dy(j-1))/(dy(j+1)+dy(j))/dy(j+1)
	  alpha_0  = -24./(dy(j-1)+dy(j-2)+dy(j))/(dy(j-1)+dy(j))/
     &                                                   dy(j+1)/dy(j)
	  alpha_m1 =  24./(dy(j-1)+dy(j-2))/dy(j-1)/dy(j)/(dy(j+1)+dy(j))
	  alpha_m2 = -24./dy(j-2)/dy(j-1)/((dy(j)**2)+2.*dy(j-1)*dy(j)
     &                +dy(j+1)*dy(j-1)+(dy(j-1)**2.)+dy(j)*dy(j+1))
	  alpha_m3 =  24./dy(j-2)/(2.*dy(j)*(dy(j-2)**2.)
     &        +4.*dy(j)*dy(j-1)*dy(j-2)+2.*dy(j+1)*dy(j-1)*dy(j-2)
     &        +3.*dy(j-1)**2.*dy(j-2)+(dy(j)**2.)*dy(j-1)
     &          +(dy(j)**2.)*dy(j-2)+2.*dy(j)*(dy(j-1)**2.)
     &          +(dy(j-2)**3.)+3.*dy(j-1)*(dy(j-2)**2.)
     &          +dy(j+1)*dy(j)*dy(j-1)+(dy(j-1)**3.)
     &          +dy(j+1)*(dy(j-1)**2.)+dy(j+1)*dy(j)*dy(j-2)
     &          +dy(j+1)*(dy(j-2)**2.))



	  DO K=0,TNKZ
          DO I=0,NKX    
		CR3(I,K,j)=CR3(I,K,j) + small*TEMP1*
     &                ( alpha_1*CU3(I,K,j+1) + 
     &                  alpha_0*CU3(I,K,j)   + 
     &                  alpha_m1*CU3(I,K,j-1) + 
     &                  alpha_m2*CU3(I,K,j-2) + 
     &                  alpha_m3*CU3(I,K,j-3)   )
	    END DO
        END DO

c	  alpha_1 = 0.0
c        alpha_0 = 0.0
c        alpha_m1 = 0.0
c        alpha_m2 = 0.0
c        alpha_m3 = 0.0

        DO K=0,TNKZ
          DO I=0,NKX    
		CR1(I,K,j)=CR1(I,K,j) + small*TEMP1*
     &                ( alpha_1*CU1(I,K,j+1) + 
     &                  alpha_0*CU1(I,K,j)   + 
     &                  alpha_m1*CU1(I,K,j-1) + 
     &                  alpha_m2*CU1(I,K,j-2) + 
     &                  alpha_m3*CU1(I,K,j-3)   )
	    END DO
        END DO

      
c For the 4th derivative of the U2 velocity, much care needs to be taken
c when performing the finite differnce derivatives.  Consider two things:
c     1.) Does the velocity exist ?
c     2.) Does dy(j) exist ?
c Below are custom finite difference methods for solving for the 4th
c derivative on a base grid.

C Lower Boudary Finite Difference Method: 4th-order Forward Finite
c Difference Approximation
	J=2
	  alpha_3 =  24./(dyf(J+3)+dyf(J+2)+dyf(J+1)+dyf(J))/
     &               (dyf(J+2)+dyf(J+1)+dyf(J))/(dyf(J+1)+dyf(J))/dyf(J)
	  alpha_2 = -24./(dyf(J+3)+dyf(J+2)+dyf(J+1))/
     &                               (dyf(J+2)+dyf(J+1))/dyf(J+1)/dyf(J)
     	  alpha_1 =  24./(dyf(J+1)+dyf(J))/dyf(J+1)/
     &                                      dyf(J+2)/(dyf(J+3)+dyf(J+2))
	  alpha_0 = -24./dyf(J+2)/(2.*dyf(j+2)*dyf(j+1)+(dyf(j+1)**2.)
     &          +(dyf(j+2)**2)+dyf(j)*dyf(j+1)+dyf(j)*dyf(j+2))/dyf(j+3)
	  alpha_m1 = 24./dyf(j+3)/(3.*dyf(j+3)*(dyf(j+2)**2.)
     &              +dyf(j+3)*(dyf(j+1)**2.)+dyf(j)*(dyf(j+3)**2.)
     &           +4.*dyf(j+3)*dyf(j+2)*dyf(j+1)+2.*(dyf(j+3)**2.)
     &              *dyf(j+1)+(dyf(j+2)**2.)*dyf(j)+2.*dyf(j+3)*dyf(j+2)
     &              *dyf(j)+dyf(j+2)*(dyf(j+1)**2.)+dyf(j+3)*dyf(j+1)
     &              *dyf(j)+(dyf(j+3)**3.)+dyf(j+2)*dyf(j+1)*dyf(j)
     &          +3.*(dyf(j+3)**2.)*dyf(j+2)+2.*(dyf(j+2)**2.)*dyf(j+1)
     &             +(dyf(j+2)**3.))  


        DO K=0,TNKZ
          DO I=0,NKX
		CR2(I,K,j)=CR2(I,K,j) + small*TEMP1*
     &                ( alpha_m1*CU2(I,K,J-1) + 
     &                  alpha_0*CU2(I,K,J) + 
     &                  alpha_1*CU2(I,K,J+1) + 
     &                  alpha_2*CU2(I,K,J+2) + 
     &                  alpha_3*CU2(I,K,J+3)   )

	    END DO
        END DO

	DO J=3,NYM-1
C 4th Order Central Finite Difference Approximation fo the 4th Derivative
C	  alpha0 ----> f(x(j+2))
C	  alpha1 ----> f(x(j+1))
C	  alpha2 ----> f(x(j  ))
C	  alpha3 ----> f(x(j-1))
C	  alpha4 ----> f(x(j-2))

	  alpha_2 =  24./(dyf(j-1)+dyf(j+1)+dyf(j+2)+dyf(j))/
     &           (dyf(j+2)+dyf(j+1)+dyf(j))/(dyf(j+2)+dyf(j+1))/dyf(j+2)
	  alpha_1 = -24./(dyf(j-1)+dyf(j)+dyf(j+1))/
     &                               (dyf(j)+dyf(j+1))/dyf(j+1)/dyf(j+2)
	  alpha_0 =  24./(dyf(j-1)+dyf(j))/dyf(j)/
     &                                      dyf(j+1)/(dyf(j+2)+dyf(j+1))
	  alpha_m1 = -24./dyf(j-1)/dyf(j)/((dyf(j)**2.)+(dyf(j+1)**2.)
     &            +dyf(j+2)*dyf(j)+2.*dyf(j+1)*dyf(j)+dyf(j+2)*dyf(j+1))
	  alpha_m2 = 24./dyf(j-1)/(3.*(dyf(j-1)**2.)*dyf(j)+(dyf(j-1)**3.)
     &          +(dyf(j)**3.)+dyf(j+2)*dyf(j+1)*dyf(j-1)+(dyf(j+1)**2.)
     &          *dyf(j)+dyf(j+2)*dyf(j+1)*dyf(j)+(dyf(j+1)**2.)*dyf(j-1)
     &           +dyf(j+2)*(dyf(j)**2.)+2.*dyf(j+1)*(dyf(j-1)**2.)
     &        +4.*dyf(j+1)*dyf(j-1)*dyf(j)+2.*dyf(j+2)*dyf(j-1)*dyf(j)
     &        +2.*dyf(j+1)*(dyf(j)**2.)+dyf(j+2)*(dyf(j-1)**2.)
     &        +3.*dyf(j-1)*(dyf(j)**2.))


        DO K=0,TNKZ
          DO I=0,NKX
		CR2(I,K,J)=CR2(I,K,J) + small*TEMP1*
     &                ( alpha_m2*CU2(I,K,J-2) + 
     &                  alpha_m1*CU2(I,K,J-1) + 
     &                  alpha_0*CU2(I,K,J)   + 
     &                  alpha_1*CU2(I,K,J+1) + 
     &                  alpha_2*CU2(I,K,J+2)   )
	    END DO
        END DO
      END DO	

C Different Difference Methods are necessary for the boundary, because
c there are not enough ghost cells to perform the 2nd order central
c finite difference approximation of the 4th dirivative.  Thus, 1st
c order forward/backward finite difference approximations are used.

C Upper Boudary Finite Difference Method: 4th-order Backwards Finite
c Difference Approximation
	J = NYM
	  alpha_1  =  24./(dyf(j)+dyf(j-2)+dyf(j+1)+dyf(j-1))/
     &             (dyf(j)+dyf(j+1)+dyf(j-1))/(dyf(j+1)+dyf(j))/dyf(j+1)
	  alpha_0  = -24./(dyf(j-1)+dyf(j-2)+dyf(j))/(dyf(j-1)+dyf(j))/
     &                                                   dyf(j+1)/dyf(j)
	  alpha_m1 =  24./(dyf(j-1)+dyf(j-2))/dyf(j-1)/
     &                                          dyf(j)/(dyf(j+1)+dyf(j))
	  alpha_m2 = -24./dyf(j-2)/dyf(j-1)/((dyf(j)**2)+2.*dyf(j-1)
     &         *dyf(j)+dyf(j+1)*dyf(j-1)+(dyf(j-1)**2.)+dyf(j)*dyf(j+1))
	  alpha_m3 =  24./dyf(j-2)/(2.*dyf(j)*(dyf(j-2)**2.)
     &        +4.*dyf(j)*dyf(j-1)*dyf(j-2)+2.*dyf(j+1)*dyf(j-1)*dyf(j-2)
     &        +3.*dyf(j-1)**2.*dyf(j-2)+(dyf(j)**2.)*dyf(j-1)
     &          +(dyf(j)**2.)*dyf(j-2)+2.*dyf(j)*(dyf(j-1)**2.)
     &          +(dyf(j-2)**3.)+3.*dyf(j-1)*(dyf(j-2)**2.)
     &          +dyf(j+1)*dyf(j)*dyf(j-1)+(dyf(j-1)**3.)
     &          +dyf(j+1)*(dyf(j-1)**2.)+dyf(j+1)*dyf(j)*dyf(j-2)
     &          +dyf(j+1)*(dyf(j-2)**2.))

        DO K=0,TNKZ
          DO I=0,NKX    
		CR2(I,K,j)=CR2(I,K,j) + small*TEMP1*
     &                ( alpha_1*CU2(I,K,j+1) + 
     &                  alpha_0*CU2(I,K,j)   + 
     &                  alpha_m1*CU2(I,K,j-1) + 
     &                  alpha_m2*CU2(I,K,j-2) + 
     &                  alpha_m3*CU2(I,K,j-3)   )
	    END DO
        END DO
	  	
	end if
	
	return
	end


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine septa_implicit_solve_u
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	include 'header_batch'
	real*8 TEMP1, alpha_m2, alpha_m1, alpha_0, alpha_1, alpha_2 
	real*8 KZ4(0:TNKZ), KX4(0:NKX)
	real*8 small 
	parameter (small = 1.0e-20)
		
	if (verbosity.gt.2) then
	  write(*,*) 'Septa-Implicit Solve for Velocity !!!'
	end if
	
c	RK Constant from Main Program
      TEMP1 = - NU * H_BAR(RK_STEP) / 2.0
	
	DO I=0,NKX
	  KX4(I) = KX2(I)*KX2(I)
	END DO
	DO K=0,TNKZ
	  KZ4(K) = KZ2(K)*KZ2(K)
	END DO
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! U2 IMPLICIT SOLVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Initialize the matrix to zeros to be used for implicit solves
C Note that the system size is NY+1, but only 1..NY are used
      DO J=0,NY+1
        DO I=0,NKX
          MATL3_C(I,J)=0.
          MATL2_C(I,J)=0.
          MATL_C(I,J) =0.
          MATD_C(I,J) =0.
          MATU_C(I,J) =0.
          MATU2_C(I,J)=0.
          MATU3_C(I,J)=0.
        END DO
      END DO

C Build implicit matrix for U2
      DO K=0,TNKZ
        DO J=2,NY
	  
C 4th Order Central Finite Difference Approximation fo the 4th Derivative
C	  alpha_2 ----> f(x(j+2))
C	  alpha_1 ----> f(x(j+1))
C	  alpha_0 ----> f(x(j  ))
C	  alpha_m1 ---> f(x(j-1))
C	  alpha_m2 ---> f(x(j-2))

	  if (J.eq.2) then
	  
	  alpha_3 =  24./(dyf(J+3)+dyf(J+2)+dyf(J+1)+dyf(J))/
     &               (dyf(J+2)+dyf(J+1)+dyf(J))/(dyf(J+1)+dyf(J))/dyf(J)
	  alpha_2 = -24./(dyf(J+3)+dyf(J+2)+dyf(J+1))/
     &                               (dyf(J+2)+dyf(J+1))/dyf(J+1)/dyf(J)
     	  alpha_1 =  24./(dyf(J+1)+dyf(J))/dyf(J+1)/
     &                                      dyf(J+2)/(dyf(J+3)+dyf(J+2))
	  alpha_0 = -24./dyf(J+2)/(2.*dyf(j+2)*dyf(j+1)+(dyf(j+1)**2.)
     &          +(dyf(j+2)**2)+dyf(j)*dyf(j+1)+dyf(j)*dyf(j+2))/dyf(j+3)
	  alpha_m1 = 24./dyf(j+3)/(3.*dyf(j+3)*(dyf(j+2)**2.)
     &              +dyf(j+3)*(dyf(j+1)**2.)+dyf(j)*(dyf(j+3)**2.)
     &           +4.*dyf(j+3)*dyf(j+2)*dyf(j+1)+2.*(dyf(j+3)**2.)
     &              *dyf(j+1)+(dyf(j+2)**2.)*dyf(j)+2.*dyf(j+3)*dyf(j+2)
     &              *dyf(j)+dyf(j+2)*(dyf(j+1)**2.)+dyf(j+3)*dyf(j+1)
     &              *dyf(j)+(dyf(j+3)**3.)+dyf(j+2)*dyf(j+1)*dyf(j)
     &          +3.*(dyf(j+3)**2.)*dyf(j+2)+2.*(dyf(j+2)**2.)*dyf(j+1)
     &             +(dyf(j+2)**3.)) 
	  alpha_m2 = 0.0
	  alpha_m3 = 0.0

	  else if(j.eq.nym) then
	  
	  alpha_3 = 0.0
	  alpha_2 = 0.0
	  alpha_1  =  24./(dyf(j)+dyf(j-2)+dyf(j+1)+dyf(j-1))/
     &             (dyf(j)+dyf(j+1)+dyf(j-1))/(dyf(j+1)+dyf(j))/dyf(j+1)
	  alpha_0  = -24./(dyf(j-1)+dyf(j-2)+dyf(j))/(dyf(j-1)+dyf(j))/
     &                                                   dyf(j+1)/dyf(j)
	  alpha_m1 =  24./(dyf(j-1)+dyf(j-2))/dyf(j-1)/
     &                                          dyf(j)/(dyf(j+1)+dyf(j))
	  alpha_m2 = -24./dyf(j-2)/dyf(j-1)/((dyf(j)**2)+2.*dyf(j-1)
     &         *dyf(j)+dyf(j+1)*dyf(j-1)+(dyf(j-1)**2.)+dyf(j)*dyf(j+1))
	  alpha_m3 =  24./dyf(j-2)/(2.*dyf(j)*(dyf(j-2)**2.)
     &        +4.*dyf(j)*dyf(j-1)*dyf(j-2)+2.*dyf(j+1)*dyf(j-1)*dyf(j-2)
     &        +3.*dyf(j-1)**2.*dyf(j-2)+(dyf(j)**2.)*dyf(j-1)
     &          +(dyf(j)**2.)*dyf(j-2)+2.*dyf(j)*(dyf(j-1)**2.)
     &          +(dyf(j-2)**3.)+3.*dyf(j-1)*(dyf(j-2)**2.)
     &          +dyf(j+1)*dyf(j)*dyf(j-1)+(dyf(j-1)**3.)
     &          +dyf(j+1)*(dyf(j-1)**2.)+dyf(j+1)*dyf(j)*dyf(j-2)
     &          +dyf(j+1)*(dyf(j-2)**2.))
     	  
	  else
	  
	  alpha_3 = 0.0
	  alpha_2 =  24./(dyf(j-1)+dyf(j+1)+dyf(j+2)+dyf(j))/
     &           (dyf(j+2)+dyf(j+1)+dyf(j))/(dyf(j+2)+dyf(j+1))/dyf(j+2)
	  alpha_1 = -24./(dyf(j-1)+dyf(j)+dyf(j+1))/
     &                               (dyf(j)+dyf(j+1))/dyf(j+1)/dyf(j+2)
	  alpha_0 =  24./(dyf(j-1)+dyf(j))/dyf(j)/
     &                                      dyf(j+1)/(dyf(j+2)+dyf(j+1))
	  alpha_m1 = -24./dyf(j-1)/dyf(j)/((dyf(j)**2.)+(dyf(j+1)**2.)
     &            +dyf(j+2)*dyf(j)+2.*dyf(j+1)*dyf(j)+dyf(j+2)*dyf(j+1))
	  alpha_m2 = 24./dyf(j-1)/(3.*(dyf(j-1)**2.)*dyf(j)+(dyf(j-1)**3.)
     &          +(dyf(j)**3.)+dyf(j+2)*dyf(j+1)*dyf(j-1)+(dyf(j+1)**2.)
     &          *dyf(j)+dyf(j+2)*dyf(j+1)*dyf(j)+(dyf(j+1)**2.)*dyf(j-1)
     &           +dyf(j+2)*(dyf(j)**2.)+2.*dyf(j+1)*(dyf(j-1)**2.)
     &        +4.*dyf(j+1)*dyf(j-1)*dyf(j)+2.*dyf(j+2)*dyf(j-1)*dyf(j)
     &        +2.*dyf(j+1)*(dyf(j)**2.)+dyf(j+2)*(dyf(j-1)**2.)
     &        +3.*dyf(j-1)*(dyf(j)**2.))
        alpha_m3 = 0.0
	  
	  end if

	    DO I=0,NKX
	      MATL3_C(I,J) = -temp1*small*alpha_m3
	      MATL2_C(I,J) = -temp1*small*alpha_m2
            MATL_C(I,J)= -TEMP1/(DYF(J-1)*DY(J)) -temp1*small*alpha_m1
            MATD_C(I,J)=1.
! Vertical Viscous term
     &         +TEMP1/(DYF(J)*DY(J)) + TEMP1/(DYF(J-1)*DY(J)) 
     &         -temp1*small*alpha_0
! Horizontal Viscous terms:
     &         +TEMP1 * (KX2(I)+KZ2(K)) - temp1*small*(KX4(I) + KZ4(K))
            MATU_C(I,J)= -TEMP1/(DYF(J)*DY(J)) -temp1*small*alpha_1
	      MATU2_C(I,J) = -temp1*small*alpha_2
	      MATU3_C(I,J) = -temp1*small*alpha_3
            VEC_C(I,J)=CR2(I,K,J)
          END DO
        END DO 

C Else, we are running in serial mode
C Set the boundary conditions for U2
        CALL APPLY_SEPTA_BC_2_LOWER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
        CALL APPLY_SEPTA_BC_2_UPPER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
C Now, solve the penta-diagonal system for U2(i,:,k)
	  CALL septa_complex(MATL3_C,MATL2_C,MATL_C,MATD_C,MATU_C,
     &                     MATU2_C,MATU3_C,VEC_C,NY+1,NKX)
	  
        DO J=1,NY+1
          DO I=0,NKX
            CU2(I,K,J)=VEC_C(I,J)
          END DO
        END DO    
      END DO 
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! U1 IMPLICIT SOLVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Solve for U1 and U3
C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)
      DO J=0,NY+1
        DO I=0,NKX
          MATL3_C(I,J)=0.
          MATL2_C(I,J)=0.
          MATL_C(I,J) =0.
          MATD_C(I,J) =0.
          MATU_C(I,J) =0.
          MATU2_C(I,J)=0.
          MATU3_C(I,J)=0.
        END DO
      END DO

      DO K=0,TNKZ
        DO J=JSTART,JEND
C 4th Order Central Finite Difference Approximation fo the 4th Derivative
C	  alpha_2  ----> f(x(j+2))
C	  alpha_1  ----> f(x(j+1))
C	  alpha_0  ----> f(x(j  ))
C	  alpha_m1 ----> f(x(j-1))
C	  alpha_m2 ----> f(x(j-2))


	  if (J.eq.2) then
	  
	  alpha_3 =  24./(dy(J+3)+dy(J+2)+dy(J+1)+dy(J))/
     &               (dy(J+2)+dy(J+1)+dy(J))/(dy(J+1)+dy(J))/dy(J)
	  alpha_2 = -24./(dy(J+3)+dy(J+2)+dy(J+1))/
     &                               (dy(J+2)+dy(J+1))/dy(J+1)/dy(J)
     	  alpha_1 =  24./(dy(J+1)+dy(J))/dy(J+1)/
     &                                      dy(J+2)/(dy(J+3)+dy(J+2))
	  alpha_0 = -24./dy(J+2)/(2.*dy(j+2)*dy(j+1)+(dy(j+1)**2.)
     &          +(dy(j+2)**2)+dy(j)*dy(j+1)+dy(j)*dy(j+2))/dy(j+3)
	  alpha_m1 = 24./dy(j+3)/(3.*dy(j+3)*(dy(j+2)**2.)
     &              +dy(j+3)*(dy(j+1)**2.)+dy(j)*(dy(j+3)**2.)
     &           +4.*dy(j+3)*dy(j+2)*dy(j+1)+2.*(dy(j+3)**2.)
     &              *dy(j+1)+(dy(j+2)**2.)*dy(j)+2.*dy(j+3)*dy(j+2)
     &              *dy(j)+dy(j+2)*(dy(j+1)**2.)+dy(j+3)*dy(j+1)
     &              *dy(j)+(dy(j+3)**3.)+dy(j+2)*dy(j+1)*dy(j)
     &          +3.*(dy(j+3)**2.)*dy(j+2)+2.*(dy(j+2)**2.)*dy(j+1)
     &             +(dy(j+2)**3.)) 
	  alpha_m2 = 0.0
	  alpha_m3 = 0.0


	  else if(j.eq.nym) then
	  
	  alpha_3 = 0.0
	  alpha_2 = 0.0
	  alpha_1  =  24./(dy(j)+dy(j-2)+dy(j+1)+dy(j-1))/
     &             (dy(j)+dy(j+1)+dy(j-1))/(dy(j+1)+dy(j))/dy(j+1)
	  alpha_0  = -24./(dy(j-1)+dy(j-2)+dy(j))/(dy(j-1)+dy(j))/
     &                                                   dy(j+1)/dy(j)
	  alpha_m1 =  24./(dy(j-1)+dy(j-2))/dy(j-1)/dy(j)/(dy(j+1)+dy(j))
	  alpha_m2 = -24./dy(j-2)/dy(j-1)/((dy(j)**2)+2.*dy(j-1)*dy(j)
     &                +dy(j+1)*dy(j-1)+(dy(j-1)**2.)+dy(j)*dy(j+1))
	  alpha_m3 =  24./dy(j-2)/(2.*dy(j)*(dy(j-2)**2.)
     &        +4.*dy(j)*dy(j-1)*dy(j-2)+2.*dy(j+1)*dy(j-1)*dy(j-2)
     &        +3.*dy(j-1)**2.*dy(j-2)+(dy(j)**2.)*dy(j-1)
     &          +(dy(j)**2.)*dy(j-2)+2.*dy(j)*(dy(j-1)**2.)
     &          +(dy(j-2)**3.)+3.*dy(j-1)*(dy(j-2)**2.)
     &          +dy(j+1)*dy(j)*dy(j-1)+(dy(j-1)**3.)
     &          +dy(j+1)*(dy(j-1)**2.)+dy(j+1)*dy(j)*dy(j-2)
     &          +dy(j+1)*(dy(j-2)**2.))

	  else
	  
	  alpha_3 = 0.0
	  alpha_2 =  24./(dy(j-1)+dy(j+1)+dy(j+2)+dy(j))/
     &           (dy(j+2)+dy(j+1)+dy(j))/(dy(j+2)+dy(j+1))/dy(j+2)
	  alpha_1 = -24./(dy(j-1)+dy(j)+dy(j+1))/
     &                               (dy(j)+dy(j+1))/dy(j+1)/dy(j+2)
	  alpha_0 =  24./(dy(j-1)+dy(j))/dy(j)/
     &                                      dy(j+1)/(dy(j+2)+dy(j+1))
	  alpha_m1 = -24./dy(j-1)/dy(j)/((dy(j)**2.)+(dy(j+1)**2.)
     &            +dy(j+2)*dy(j)+2.*dy(j+1)*dy(j)+dy(j+2)*dy(j+1))
	  alpha_m2 =  24./dy(j-1)/(3.*(dy(j-1)**2.)*dy(j)+(dy(j-1)**3.)
     &          +(dy(j)**3.)+dy(j+2)*dy(j+1)*dy(j-1)+(dy(j+1)**2.)
     &          *dy(j)+dy(j+2)*dy(j+1)*dy(j)+(dy(j+1)**2.)*dy(j-1)
     &           +dy(j+2)*(dy(j)**2.)+2.*dy(j+1)*(dy(j-1)**2.)
     &        +4.*dy(j+1)*dy(j-1)*dy(j)+2.*dy(j+2)*dy(j-1)*dy(j)
     &        +2.*dy(j+1)*(dy(j)**2.)+dy(j+2)*(dy(j-1)**2.)
     &        +3.*dy(j-1)*(dy(j)**2.))
        alpha_m3 = 0.0

	  end if

c	  alpha_3 = 0.0
c        alpha_2 = 0.0
c        alpha_1 = 0.0
c        alpha_0 = 0.0
c        alpha_m1 = 0.0
c        alpha_m2 = 0.0
c	  alpha_m3 = 0.0
	    
          DO I=0,NKX
	      MATL3_C(I,J) = - temp1*small*alpha_m3
	      MATL2_C(I,J) = - temp1*small*alpha_m2
            MATL_C(I,J) = -TEMP1/(DY(J)*DYF(J)) - temp1*small*alpha_m1
            MATD_C(I,J) = 1.
! Vertical Viscous terms:
     &         -TEMP1*(-1./(DY(J+1)*DYF(J))-1./(DY(J)*DYF(J)))
     &         - temp1*small*alpha0
! Horizontal Viscous term:
     &         +TEMP1*(KX2(I)+KZ2(K)) - temp1*small*(KX4(I) + KZ4(K))
            MATU_C(I,J) = -TEMP1/(DY(J+1)*DYF(J)) - temp1*small*alpha_1
	      MATU2_C(I,J) = - temp1*small*alpha_2
	      MATU3_C(I,J) = - temp1*small*alpha_3
            VEC_C(I,J)=CR1(I,K,J)
          END DO
        END DO


C Else, we are running in serial mode
C Set the boundary conditions for U1
        CALL APPLY_SEPTA_BC_1_LOWER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
        CALL APPLY_SEPTA_BC_1_UPPER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
C Now, solve the penta-diagonal system for U1(:,k,:)
	  CALL septa_complex(MATL3_C,MATL2_C,MATL_C,MATD_C,MATU_C,
     &                     MATU2_C,MATU3_C,VEC_C,NY+1,NKX)
	 
C Now, solve the penta-diagonal system for U1(i,k,:)
        DO J=0,NY+1
          DO I=0,NKX
            CU1(I,K,J)=VEC_C(I,J)
          END DO
        END DO
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! U3 IMPLICIT SOLVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO K=0,TNKZ
        DO J=JSTART,JEND
C 4th Order Central Finite Difference Approximation fo the 4th Derivative
C	  alpha_2  ----> f(x(j+2))
C	  alpha_1  ----> f(x(j+1))
C	  alpha_0  ----> f(x(j  ))
C	  alpha_m1 ----> f(x(j-1))
C	  alpha_m2 ----> f(x(j-2))

	  if (J.eq.2) then
	  
	  alpha_3 =  24./(dy(J+3)+dy(J+2)+dy(J+1)+dy(J))/
     &               (dy(J+2)+dy(J+1)+dy(J))/(dy(J+1)+dy(J))/dy(J)
	  alpha_2 = -24./(dy(J+3)+dy(J+2)+dy(J+1))/
     &                               (dy(J+2)+dy(J+1))/dy(J+1)/dy(J)
     	  alpha_1 =  24./(dy(J+1)+dy(J))/dy(J+1)/
     &                                      dy(J+2)/(dy(J+3)+dy(J+2))
	  alpha_0 = -24./dy(J+2)/(2.*dy(j+2)*dy(j+1)+(dy(j+1)**2.)
     &          +(dy(j+2)**2)+dy(j)*dy(j+1)+dy(j)*dy(j+2))/dy(j+3)
	  alpha_m1 = 24./dy(j+3)/(3.*dy(j+3)*(dy(j+2)**2.)
     &              +dy(j+3)*(dy(j+1)**2.)+dy(j)*(dy(j+3)**2.)
     &           +4.*dy(j+3)*dy(j+2)*dy(j+1)+2.*(dy(j+3)**2.)
     &              *dy(j+1)+(dy(j+2)**2.)*dy(j)+2.*dy(j+3)*dy(j+2)
     &              *dy(j)+dy(j+2)*(dy(j+1)**2.)+dy(j+3)*dy(j+1)
     &              *dy(j)+(dy(j+3)**3.)+dy(j+2)*dy(j+1)*dy(j)
     &          +3.*(dy(j+3)**2.)*dy(j+2)+2.*(dy(j+2)**2.)*dy(j+1)
     &             +(dy(j+2)**3.)) 
	  alpha_m2 = 0.0
	  alpha_m3 = 0.0

	  else if(j.eq.nym) then
	  
	  alpha_3 = 0.0
	  alpha_2 = 0.0
	  alpha_1  =  24./(dy(j)+dy(j-2)+dy(j+1)+dy(j-1))/
     &             (dy(j)+dy(j+1)+dy(j-1))/(dy(j+1)+dy(j))/dy(j+1)
	  alpha_0  = -24./(dy(j-1)+dy(j-2)+dy(j))/(dy(j-1)+dy(j))/
     &                                                   dy(j+1)/dy(j)
	  alpha_m1 =  24./(dy(j-1)+dy(j-2))/dy(j-1)/dy(j)/(dy(j+1)+dy(j))
	  alpha_m2 = -24./dy(j-2)/dy(j-1)/((dy(j)**2)+2.*dy(j-1)*dy(j)
     &                +dy(j+1)*dy(j-1)+(dy(j-1)**2.)+dy(j)*dy(j+1))
	  alpha_m3 =  24./dy(j-2)/(2.*dy(j)*(dy(j-2)**2.)
     &        +4.*dy(j)*dy(j-1)*dy(j-2)+2.*dy(j+1)*dy(j-1)*dy(j-2)
     &        +3.*dy(j-1)**2.*dy(j-2)+(dy(j)**2.)*dy(j-1)
     &          +(dy(j)**2.)*dy(j-2)+2.*dy(j)*(dy(j-1)**2.)
     &          +(dy(j-2)**3.)+3.*dy(j-1)*(dy(j-2)**2.)
     &          +dy(j+1)*dy(j)*dy(j-1)+(dy(j-1)**3.)
     &          +dy(j+1)*(dy(j-1)**2.)+dy(j+1)*dy(j)*dy(j-2)
     &          +dy(j+1)*(dy(j-2)**2.))
     	  
	  else
	  
	  alpha_3 = 0.0
	  alpha_2 =  24./(dy(j-1)+dy(j+1)+dy(j+2)+dy(j))/
     &           (dy(j+2)+dy(j+1)+dy(j))/(dy(j+2)+dy(j+1))/dy(j+2)
	  alpha_1 = -24./(dy(j-1)+dy(j)+dy(j+1))/
     &                               (dy(j)+dy(j+1))/dy(j+1)/dy(j+2)
	  alpha_0 =  24./(dy(j-1)+dy(j))/dy(j)/
     &                                      dy(j+1)/(dy(j+2)+dy(j+1))
	  alpha_m1 = -24./dy(j-1)/dy(j)/((dy(j)**2.)+(dy(j+1)**2.)
     &            +dy(j+2)*dy(j)+2.*dy(j+1)*dy(j)+dy(j+2)*dy(j+1))
	  alpha_m2 =  24./dy(j-1)/(3.*(dy(j-1)**2.)*dy(j)+(dy(j-1)**3.)
     &          +(dy(j)**3.)+dy(j+2)*dy(j+1)*dy(j-1)+(dy(j+1)**2.)
     &          *dy(j)+dy(j+2)*dy(j+1)*dy(j)+(dy(j+1)**2.)*dy(j-1)
     &           +dy(j+2)*(dy(j)**2.)+2.*dy(j+1)*(dy(j-1)**2.)
     &        +4.*dy(j+1)*dy(j-1)*dy(j)+2.*dy(j+2)*dy(j-1)*dy(j)
     &        +2.*dy(j+1)*(dy(j)**2.)+dy(j+2)*(dy(j-1)**2.)
     &        +3.*dy(j-1)*(dy(j)**2.))
        alpha_m3 = 0.0
	  end if

	  alpha_3 = 0.0
        alpha_2 = 0.0
        alpha_1 = 0.0
        alpha_0 = 0.0
        alpha_m1 = 0.0
        alpha_m2 = 0.0
	  alpha_m3 = 0.0
    
          DO I=0,NKX
	      MATL3_C(I,J) = - temp1*small*alpha_m3
	      MATL2_C(I,J) = - temp1*small*alpha_m2
		MATL_C(I,J)=-TEMP1/(DY(J)*DYF(J)) - temp1*small*alpha_m1
            MATD_C(I,J)=1.
! Vertical Viscous terms:
     &         -TEMP1*(-1./(DY(J+1)*DYF(J))-1./(DY(J)*DYF(J)))
     &         - temp1*small*alpha_0
! Horizontal Viscous term:
     &         +TEMP1*(KX2(I)+KZ2(K)) - temp1*small*(KX4(I) + KZ4(K))
            MATU_C(I,J)= -TEMP1/(DY(J+1)*DYF(J)) - temp1*small*alpha_1
	      MATU2_C(I,J) = -temp1*small*alpha_2
	      MATU3_C(I,J) = -temp1*small*alpha_3
            VEC_C(I,J)=CR3(I,K,J)
          END DO
        END DO


C Else, we are running in serial mode
C Set the boundary conditions for U3
	  CALL APPLY_SEPTA_BC_3_LOWER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
	  CALL APPLY_SEPTA_BC_3_UPPER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
C Now, solve the penta-diagonal system for U3(i,:,k)
	  CALL septa_complex(MATL3_C,MATL2_C,MATL_C,MATD_C,MATU_C,
     &                     MATU2_C,MATU3_C,VEC_C,NY+1,NKX)

C Now, solve the penta-diagonal system for U3(i,k,:)
        DO J=0,NY+1
          DO I=0,NKX
            CU3(I,K,J)=VEC_C(I,J)
          END DO
        END DO
      END DO
	
	return
	end


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE penta_real(A,B,C,D,E,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Uses the Thomas algorithm to solve Ax=b for pentadiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x (in storage G).
C The indexing should be done by ROW, ie.
C [ c1  d1  e1   0   0   0 ...
C [ b2  c2  d2  e2   0   0 ...
C [ a3  b3  c3  d3  e3   0 ...
C [  0  a4  b4  c4  d4  e4 ...
C [  0   0  a5  b5  c5  d5 ...
C [  0   0   0  a6  b6  c6 ...

      INTEGER I, J, NX, NY
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY)
      REAL*8 D(0:NX,0:NY), E(0:NX,0:NY), G(0:NX,0:NY)

C	FORWARD SUBSTITUTION
      DO J=0,NY-2
        DO I=0,NX
C         FIRST LINE
          B(I,J+1) = -B(I,J+1)/C(I,J)
          C(I,J+1) =  C(I,J+1) + B(I,J+1)*D(I,J)
	    D(I,J+1) =  D(I,J+1) + B(I,J+1)*E(I,J)
          G(I,J+1) =  G(I,J+1) + B(I,J+1)*G(I,J)
	    
C	    SECOND LINE
	    A(I,J+2) = -A(I,J+2)/C(I,J)
	    B(I,J+2) =  B(I,J+2) + A(I,J+2)*D(I,J)
	    C(I,J+2) =  C(I,J+2) + A(I,J+2)*E(I,J)
	    G(I,J+2) =  G(I,J+2) + A(I,J+2)*G(I,J)
        END DO
      END DO
	DO I=0,NX
	  B(I,NY) = -B(I,NY)/C(I,NY-1)
	  C(I,NY) =  C(I,NY) + B(I,NY)*D(I,NY-1)
	  G(I,NY) =  G(I,NY) + B(I,NY)*G(I,NY-1)
	END DO
C	BACKWARDS SUBSTITUTION
      DO I=0,NX
        G(I,NY)   = G(I,NY)/C(I,NY)
	  G(I,NY-1) = (G(I,NY-1)-D(I,NY-1)*G(I,NY))/C(I,NY-1)
      END DO
      DO J=NY-2,0,-1
        DO I=0,NX
          G(I,J)=(G(I,J) - D(I,J)*G(I,J+1) - E(I,J)*G(I,J+2))/C(I,J)
        END DO
      END DO
    
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE penta_complex(A,B,C,D,E,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Uses the Thomas algorithm to solve Ax=b for pentadiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x (in storage G).
C The indexing should be done by ROW, ie.
C [ c1  d1  e1   0   0   0 ...
C [ b2  c2  d2  e2   0   0 ...
C [ a3  b3  c3  d3  e3   0 ...
C [  0  a4  b4  c4  d4  e4 ...
C [  0   0  a5  b5  c5  d5 ...
C [  0   0   0  a6  b6  c6 ...

      INTEGER I, J, NX, NY
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY)
      REAL*8 D(0:NX,0:NY), E(0:NX,0:NY)
	COMPLEX*16 G(0:NX,0:NY)

C	FORWARD SUBSTITUTION
      DO J=0,NY-2
        DO I=0,NX
C         FIRST LINE
          B(I,J+1) = -B(I,J+1)/C(I,J)
          C(I,J+1) =  C(I,J+1) + B(I,J+1)*D(I,J)
	    D(I,J+1) =  D(I,J+1) + B(I,J+1)*E(I,J)
          G(I,J+1) =  G(I,J+1) + B(I,J+1)*G(I,J)
	    
C	    SECOND LINE
	    A(I,J+2) = -A(I,J+2)/C(I,J)
	    B(I,J+2) =  B(I,J+2) + A(I,J+2)*D(I,J)
	    C(I,J+2) =  C(I,J+2) + A(I,J+2)*E(I,J)
	    G(I,J+2) =  G(I,J+2) + A(I,J+2)*G(I,J)
        END DO
      END DO
	DO I=0,NX
	  B(I,NY) = -B(I,NY)/C(I,NY-1)
	  C(I,NY) =  C(I,NY) + B(I,NY)*D(I,NY-1)
	  G(I,NY) =  G(I,NY) + B(I,NY)*G(I,NY-1)
	END DO
C	BACKWARDS SUBSTITUTION
      DO I=0,NX
        G(I,NY)   = G(I,NY)/C(I,NY)
	  G(I,NY-1) = (G(I,NY-1)-D(I,NY-1)*G(I,NY))/C(I,NY-1)
      END DO
      DO J=NY-2,0,-1
        DO I=0,NX
          G(I,J)=(G(I,J) - D(I,J)*G(I,J+1) - E(I,J)*G(I,J+2))/C(I,J)
        END DO
      END DO
    
      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE septa_complex(A,B,C,D,E,F,G,H,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Uses the Thomas algorithm to solve Ax=b for pentadiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x (in storage H).
C The indexing should be done by ROW, ie.
C [ d1  e1  f1  g1   0   0   0   0   0   0 ...
C [ c2  d2  e2  f2  g2   0   0   0   0   0 ...
C [ b3  c3  d3  e3  f3  g3   0   0   0   0 ...
C [ a4  b4  c4  d4  e4  f4  g4   0   0   0 ...
C [  0  a5  b5  c5  d5  e5  f5  g5   0   0 ...
C [  0   0  a6  b6  c6  d6  e6  f6  g6   0 ...
C [  0   0   0  a7  b7  c7  d7  e7  f7  g7 ...

      INTEGER I, J, NX, NY
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY)
      REAL*8 D(0:NX,0:NY), E(0:NX,0:NY), F(0:NX,0:NY)
	REAL*8 G(0:NX,0:NY)
	complex*16 H(0:NX,0:NY)
	
C	FORWARD SUBSTITUTION
      DO J=0,NY-3
        DO I=0,NX
C         FIRST LINE
	    c(I,J+1) = -c(I,J+1)/d(I,J)
	    d(I,J+1) =  d(I,J+1) + c(I,J+1)*e(I,J)
	    e(I,J+1) =  e(I,J+1) + c(I,J+1)*f(I,J)
	    f(I,J+1) =  f(I,J+1) + c(I,J+1)*g(I,J)
	    h(I,J+1) =  h(I,J+1) + c(I,J+1)*h(I,J)

C         Second Line
	    b(I,J+2) = -b(I,J+2)/d(I,J)
	    c(I,J+2) =  c(I,J+2) + b(I,J+2)*e(I,J)
	    d(I,J+2) =  d(I,J+2) + b(I,J+2)*f(I,J)
	    e(I,J+2) =  e(I,J+2) + b(I,J+2)*g(I,J)
	    h(I,J+2) =  h(I,J+2) + b(I,J+2)*h(I,J)
   
C	    Third Line
	    a(I,J+3) = -a(I,J+3)/d(I,J)
	    b(I,J+3) =  b(I,J+3) + a(I,J+3)*e(I,J)
	    c(I,J+3) =  c(I,J+3) + a(I,J+3)*f(I,J)
	    d(I,J+3) =  d(I,J+3) + a(I,J+3)*g(I,J)
	    h(I,J+3) =  h(I,J+3) + a(I,J+3)*h(I,J)
        END DO
      END DO
	DO I=0,NX
C	  First Line Clean Up (n-2)
	  c(i,NY-1) = -c(i,NY-1)/d(i,NY-2)
	  d(i,NY-1) =  d(i,NY-1) + c(i,NY-1)*e(i,NY-2)
	  e(i,NY-1) =  e(i,NY-1) + c(i,NY-1)*f(i,NY-2)
	  h(i,NY-1) =  h(i,NY-1) + c(i,NY-1)*h(i,NY-2)

C	  Second Line Clean Up (n-2)
	  b(i,NY) = -b(i,NY)/d(i,NY-2)
	  c(i,NY) =  c(i,NY) + b(i,NY)*e(i,NY-2)
	  d(i,NY) =  d(i,NY) + b(i,NY)*f(i,NY-2)
	  h(i,NY) =  h(i,NY) + b(i,NY)*h(i,NY-2)

C	  First Line Clean Up (n-1)
	  c(i,NY) = -c(i,NY)/d(i,NY-1)
	  d(i,NY) =  d(i,NY) + c(i,NY)*e(i,NY-1)
	  h(i,NY) =  h(i,NY) + c(i,NY)*h(i,NY-1)
	END DO
	
C	BACKWARDS SUBSTITUTION
      DO I=0,NX
	  h(i,NY)   = h(i,NY)/d(i,NY)
	  h(i,NY-1) = (h(i,NY-1)-e(i,NY-1)*h(i,NY)             )/d(i,NY-1)
	  h(i,NY-2) = (h(i,NY-2)-e(i,NY-2)*h(i,NY-1)
     &                                     -f(i,NY-2)*h(i,NY))/d(i,NY-2)
      END DO
      DO J=NY-3,0,-1
        DO I=0,NX
	    h(i,j)=(h(i,j)-e(i,j)*h(i,j+1)
     &                          -f(i,j)*h(i,j+2)-g(i,j)*h(i,j+3))/d(i,j)
	  END DO
      END DO
    
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE septa_real(A,B,C,D,E,F,G,H,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Uses the Thomas algorithm to solve Ax=b for pentadiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x (in storage H).
C The indexing should be done by ROW, ie.
C [ d1  e1  f1  g1   0   0   0   0   0   0 ...
C [ c2  d2  e2  f2  g2   0   0   0   0   0 ...
C [ b3  c3  d3  e3  f3  g3   0   0   0   0 ...
C [ a4  b4  c4  d4  e4  f4  g4   0   0   0 ...
C [  0  a5  b5  c5  d5  e5  f5  g5   0   0 ...
C [  0   0  a6  b6  c6  d6  e6  f6  g6   0 ...
C [  0   0   0  a7  b7  c7  d7  e7  f7  g7 ...

      INTEGER I, J, NX, NY
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY)
      REAL*8 D(0:NX,0:NY), E(0:NX,0:NY), F(0:NX,0:NY)
	REAL*8 G(0:NX,0:NY)
	real*8 H(0:NX,0:NY)
	
C	FORWARD SUBSTITUTION
      DO J=0,NY-3
        DO I=0,NX
C         FIRST LINE
	    c(I,J+1) = -c(I,J+1)/d(I,J)
	    d(I,J+1) =  d(I,J+1) + c(I,J+1)*e(I,J)
	    e(I,J+1) =  e(I,J+1) + c(I,J+1)*f(I,J)
	    f(I,J+1) =  f(I,J+1) + c(I,J+1)*g(I,J)
	    h(I,J+1) =  h(I,J+1) + c(I,J+1)*h(I,J)

C         Second Line
	    b(I,J+2) = -b(I,J+2)/d(I,J)
	    c(I,J+2) =  c(I,J+2) + b(I,J+2)*e(I,J)
	    d(I,J+2) =  d(I,J+2) + b(I,J+2)*f(I,J)
	    e(I,J+2) =  e(I,J+2) + b(I,J+2)*g(I,J)
	    h(I,J+2) =  h(I,J+2) + b(I,J+2)*h(I,J)
   
C	    Third Line
	    a(I,J+3) = -a(I,J+3)/d(I,J)
	    b(I,J+3) =  b(I,J+3) + a(I,J+3)*e(I,J)
	    c(I,J+3) =  c(I,J+3) + a(I,J+3)*f(I,J)
	    d(I,J+3) =  d(I,J+3) + a(I,J+3)*g(I,J)
	    h(I,J+3) =  h(I,J+3) + a(I,J+3)*h(I,J)
        END DO
      END DO
	DO I=0,NX
C	  First Line Clean Up (n-2)
	  c(i,NY-1) = -c(i,NY-1)/d(i,NY-2)
	  d(i,NY-1) =  d(i,NY-1) + c(i,NY-1)*e(i,NY-2)
	  e(i,NY-1) =  e(i,NY-1) + c(i,NY-1)*f(i,NY-2)
	  h(i,NY-1) =  h(i,NY-1) + c(i,NY-1)*h(i,NY-2)

C	  Second Line Clean Up (n-2)
	  b(i,NY) = -b(i,NY)/d(i,NY-2)
	  c(i,NY) =  c(i,NY) + b(i,NY)*e(i,NY-2)
	  d(i,NY) =  d(i,NY) + b(i,NY)*f(i,NY-2)
	  h(i,NY) =  h(i,NY) + b(i,NY)*h(i,NY-2)

C	  First Line Clean Up (n-1)
	  c(i,NY) = -c(i,NY)/d(i,NY-1)
	  d(i,NY) =  d(i,NY) + c(i,NY)*e(i,NY-1)
	  h(i,NY) =  h(i,NY) + c(i,NY)*h(i,NY-1)
	END DO
	
C	BACKWARDS SUBSTITUTION
      DO I=0,NX
	  h(i,NY)   = h(i,NY)/d(i,NY)
	  h(i,NY-1) = (h(i,NY-1)-e(i,NY-1)*h(i,NY)             )/d(i,NY-1)
	  h(i,NY-2) = (h(i,NY-2)-e(i,NY-2)*h(i,NY-1)
     &                                     -f(i,NY-2)*h(i,NY))/d(i,NY-2)
      END DO
      DO J=NY-3,0,-1
        DO I=0,NX
	    h(i,j)=(h(i,j)-e(i,j)*h(i,j+1)
     &                          -f(i,j)*h(i,j+2)-g(i,j)*h(i,j+3))/d(i,j)
	  END DO
      END DO
    
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine backwards_constants(TEMP1,TEMP2,TEMP3,TEMP4,TEMP5)
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5
	INCLUDE 'header_batch'
c	write(*,*) '!!! Backwards Constants !!!'

      TEMP1 = - TEMP1
      TEMP2 = - TEMP2
      TEMP3 = - TEMP3
      TEMP4 = - TEMP4
      TEMP5 = - TEMP5

	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_SEPTA_BC_1_LOWER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I

C Bottom Wall:
      IF (U_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL2_C(I,0)=0. 
          MATL_C(I,0)=0. 
          MATD_C(I,0)=1.
          MATU_C(I,0)=0.                   
          MATU2_C(I,0)=0.                   
          VEC_C(I,0)=0.

          MATL2_C(I,1)=0. 
          MATL_C(I,1)=0. 
          MATD_C(I,1)=1.
          MATU_C(I,1)=0.                   
          MATU2_C(I,1)=0.                   
          VEC_C(I,1)=U_BC_YMIN_C1 
        END DO
      ELSE
C Neumann
	WRITE(*,*) 'YOU IDIOT: NEUMANN BCs NOT SUPPORTED !!!'
C        DO I=0,NKX
C          MATL_C(I,0)=0.
C          MATD_C(I,0)=-1.
C          MATU_C(I,0)=1.
C          VEC_C(I,0)=DY(1)*U_BC_YMIN_C1
C        END DO
      END IF

      RETURN 
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_SEPTA_BC_1_UPPER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I

C Top wall
      IF (U_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL3_C(I,NY+1)=0.
          MATL2_C(I,NY+1)=0.
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          MATU2_C(I,NY+1)=0.
          MATU3_C(I,NY+1)=0.
          VEC_C(I,NY+1)=0.

          MATL3_C(I,NY)=0.
          MATL2_C(I,NY)=0.
          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
          MATU2_C(I,NY)=0.
          MATU3_C(I,NY)=0.
          VEC_C(I,NY)=U_BC_YMAX_C1
        END DO
      ELSE
C Neumann
	WRITE(*,*) 'YOU IDIOT: NEUMANN BCs NOT SUPPORTED !!!'
C        DO I=0,NKX
C          MATL_C(I,NY+1)=-1.
C          MATD_C(I,NY+1)=1.
C          MATU_C(I,NY+1)=0.
C          VEC_C(I,NY+1)=DY(NY+1)*U_BC_YMAX_C1
C        END DO      
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_SEPTA_BC_2_LOWER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I

C Bottom Wall:
      IF (V_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL3_C(I,1)=0.d0 
          MATL2_C(I,1)=0.d0 
          MATL_C(I,1)=0.d0 
          MATD_C(I,1)=1.d0
          MATU_C(I,1)=0.d0                   
          MATU2_C(I,1)=0.d0                   
          MATU3_C(I,1)=0.d0                   
          VEC_C(I,1)=V_BC_YMIN_C1 

          MATL3_C(I,2)=0.d0 
          MATL2_C(I,2)=0.d0 
          MATL_C(I,2)=0.d0 
          MATD_C(I,2)=1.d0
          MATU_C(I,2)=0.d0                   
          MATU2_C(I,2)=0.d0                   
          MATU3_C(I,2)=0.d0                   
          VEC_C(I,2)=V_BC_YMIN_C1 
        END DO
      ELSE IF (V_BC_YMIN.EQ.1) THEN
C Neumann
	WRITE(*,*) 'YOU IDIOT: NEUMANN BCs NOT SUPPORTED !!!'
C        DO I=0,NKX
C          MATD_C(I,1)=-1.d0
C          MATU_C(I,1)=1.d0
C          MATL_C(I,1)=0.d0
C          VEC_C(I,1)=DYF(1)*V_BC_YMIN_C1
C        END DO
      END IF

C The following is only a placeholder, this row is used for U1 and U3
      DO I=0,NKX
        MATL2_C(I,0) = 0.
        MATL_C(I,0) = 0.
        MATD_C(I,0) = 1.
        MATU_C(I,0) = 0.
        MATU2_C(I,0) = 0.
        VEC_C(I,0) = 0.
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_SEPTA_BC_2_UPPER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I
C Top wall
      IF (V_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL3_C(I,NY+1)=0.
          MATL2_C(I,NY+1)=0.
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          MATU2_C(I,NY+1)=0.
          MATU3_C(I,NY+1)=0.
          VEC_C(I,NY+1)=V_BC_YMAX_C1
          
          MATL3_C(I,NY)=0.
          MATL2_C(I,NY)=0.
          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
          MATU2_C(I,NY)=0.
          MATU3_C(I,NY)=0.
          VEC_C(I,NY)=V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.EQ.1) THEN
C Neumann
	WRITE(*,*) 'YOU IDIOT: NEUMANN BCs NOT SUPPORTED !!!'
C        DO I=0,NKX
C          MATL_C(I,NY+1)=-1.
C          MATD_C(I,NY+1)=1.
C          VEC_C(I,NY+1)=DYF(NY)*V_BC_YMAX_C1
C        END DO      
      END IF
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_SEPTA_BC_3_LOWER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I

C Bottom Wall:
      IF (W_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL3_C(I,0)=0. 
          MATL2_C(I,0)=0. 
          MATL_C(I,0)=0. 
          MATD_C(I,0)=1.
          MATU_C(I,0)=0.                   
          MATU2_C(I,0)=0.                   
          MATU3_C(I,0)=0.                   
          VEC_C(I,0)=0.

          MATL3_C(I,1)=0. 
          MATL2_C(I,1)=0. 
          MATL_C(I,1)=0. 
          MATD_C(I,1)=1.
          MATU_C(I,1)=0.                   
          MATU2_C(I,1)=0.                   
          MATU3_C(I,1)=0.                   
          VEC_C(I,1)=W_BC_YMIN_C1
        END DO
      ELSE
C Neumann
	WRITE(*,*) 'YOU IDIOT: NEUMANN BCs NOT SUPPORTED !!!'
C        DO I=0,NKX
C          MATL_C(I,0)=0.
C          MATD_C(I,0)=-1.
C          MATU_C(I,0)=1.
C          VEC_C(I,0)=DY(1)*W_BC_YMIN_C1
C        END DO
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_SEPTA_BC_3_UPPER_C(MATL2_C,MATL_C,MATD_C,MATU_C,
     &                                                    MATU2_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I

C Top wall
      IF (W_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL3_C(I,NY+1)=0.
          MATL2_C(I,NY+1)=0.
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
	    MATU2_C(I,NY+1)=0.
	    MATU3_C(I,NY+1)=0.
          VEC_C(I,NY+1)=0.

          MATL3_C(I,NY)=0.
          MATL2_C(I,NY)=0.
          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
	    MATU2_C(I,NY)=0.
	    MATU3_C(I,NY)=0.
          VEC_C(I,NY)=W_BC_YMAX_C1
        END DO
      ELSE
C Neumann
	WRITE(*,*) 'YOU IDIOT: NEUMANN BCs NOT SUPPORTED !!!'
C        DO I=0,NKX
C          MATL_C(I,NY+1)=-1.
C          MATD_C(I,NY+1)=1.
C          MATU_C(I,NY+1)=0.
C          VEC_C(I,NY+1)=DY(NY+1)*W_BC_YMAX_C1
C        END DO      
      END IF

      RETURN
      END

	
