C******************************************************************************|
C channel.f, the channel-flow solvers for diablo.                  VERSION 0.9
C These solvers were written primarily by John R. Taylor (spring 2005).
C******************************************************************************|
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Initialize any constants here
      INCLUDE 'header'
      INTEGER J, N

! Defined starting and ending indices in the wall-bounded direction
        JSTART=2
!        JSTART_TH(1)=1     
        JEND=NY-1
!        JEND_TH(1)=NY
	
	WRITE(*,*) 'RHO_TYPE: ',RHO_TYPE

	IF (RHO_TYPE .EQ. 1) THEN
! shear layer holmboe instabilitities no normal flux
	  JSTART_TH(1) = 1
	  JEND_TH(1) = NY
	ELSE IF (RHO_TYPE .EQ. 2) THEN
! constant gradient fixed density
	  JSTART_TH(1) = 2
	  JEND_TH(1) = NY-1	
	ELSE
	  STOP 'Undefined background density profile. Use type 1 or 2.'
	END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_1
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Main time-stepping algorithm for the channel-flow case.
C This algorithm uses Crank-Nicolson for all terms involving vertical
C derivatives (viscous and nonlinear) and 3rd order Runge-Kutta for the
C rest of the terms
C INPUTS  (in Fourier space):  CUi, CP, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, CP, and (if k<3) CFi at (k)
C Each RK step, there are 14 FFT calls. 11 storage variables are used.     
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      implicit none
      INCLUDE 'header'

      INTEGER I,J,K,N      
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, UBULK

C Communicate the information between ghost cells
      IF (USE_MPI) THEN
        CALL GHOST_CHAN_MPI
      END IF

C Define the constants that are used in the time-stepping
C For reference, see Numerical Renaissance
      TEMP1=NU * H_BAR(RK_STEP) / 2.0
      TEMP2=H_BAR(RK_STEP) / 2.0
      TEMP3=ZETA_BAR(RK_STEP) * H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)
      TEMP5=BETA_BAR(RK_STEP) * H_BAR(RK_STEP)

C First, we will compute the explicit RHS terms and store in Ri
C Note, Momentum equation and hence the RHS is evaluated at the
C corresponding velocity points.

C Store the old velocity in the RHS vector
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CU1(I,K,J)
            CR3(I,K,J)=CU3(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY 
        DO K=0,TNKZ
          DO I=0,NKX
            CR2(I,K,J)=CU2(I,K,J)
          END DO
        END DO
      END DO

C Add the R-K term from the rk-1 step 
      IF (RK_STEP .GT. 1) THEN
        DO J=JSTART,JEND
          DO K=0,TNKZ
            DO I=0,NKX
              CR1(I,K,J)=CR1(I,K,J)+TEMP3*CF1(I,K,J)
              CR3(I,K,J)=CR3(I,K,J)+TEMP3*CF3(I,K,J)
            END DO
          END DO
        END DO
        DO J=2,NY
          DO K=0,TNKZ
            DO I=0,NKX
              CR2(I,K,J)=CR2(I,K,J)+TEMP3*CF2(I,K,J)
            END DO
          END DO
        END DO
      END IF
          

C Take the y-derivative of the pressure at GY points in Fourier space
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=(CP(I,K,J) - CP(I,K,J-1)) / DY(J)
          END DO
        END DO
      END DO

C Add the pressure gradient to the RHS as explicit Euler
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CR1(I,K,J)-TEMP4*(CIKX(I)*CP(I,K,J))
            CR3(I,K,J)=CR3(I,K,J)-TEMP4*(CIKZ(K)*CP(I,K,J))
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CR2(I,K,J)=CR2(I,K,J)-TEMP4*CS1(I,K,J)
          END DO
        END DO
      END DO

C Here, add the constant, forcing pressure gradient (if needed)
 
C Now compute the term R-K term Ai
C Compile terms of Ai in CFi which will be saved for next time step
C First, store the horizontal viscous terms in CFi
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=-NU * KX2(I) * CU1(I,K,J) 
     &            - NU * KZ2(K) * CU1(I,K,J)
            CF3(I,K,J)=-NU * KX2(I) * CU3(I,K,J) 
     &            - NU * KZ2(K) * CU3(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=-NU * KX2(I) * CU2(I,K,J) 
     &            - NU * KZ2(K) * CU2(I,K,J)
          END DO 
        END DO
      END DO

! Do for each scalar
      DO N=1,N_TH
! If a scalar contributes to the denisty, RI_TAU is not equal to zero and
! add the buoyancy term as explicit R-K.  Don't add the 0,0 mode which 
! corresponds to a plane average.  The plane averaged density balances
! the hydrostratic pressure component.
      DO J=2,NY
        DO K=1,TNKZ
          DO I=1,NKX
! Use second order interpolation
             CF2(I,K,J)=CF2(I,K,J)-RI_TAU(N)*
     &      (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
          END DO
        END DO
        K=0
        DO I=1,NKX 
             CF2(I,K,J)=CF2(I,K,J)-RI_TAU(N)*
     &      (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
        END DO
        I=0
        DO K=1,TNKZ
             CF2(I,K,J)=CF2(I,K,J)-RI_TAU(N)*
     &      (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
        END DO
      END DO

! Now, compute the RHS vector for the scalar equations
! Since TH is defined at horizontal velocity points, the
! scalar update equation will be very similar to the horizontal
! velocity update.

! We will store the RHS scalar terms in CRTH, RTH
! The k-1 term for the R-K stepping is saved in FTH, CFTH

! First, build the RHS vector, use CRTH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
          CRTH(I,K,J,N)=CTH(I,K,J,N)
         ENDDO
       END DO
      END DO
! Add term from k-2 step to free up CFTH variable
      IF (RK_STEP .GT. 1) THEN
        DO J=JSTART_TH(N),JEND_TH(N)
          DO K=0,TNKZ
            DO I=0,NKX
              CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP3*CFTH(I,K,J,N)
            END DO
          END DO
        END DO
       END IF

! Now compute the explicit R-K term Ai
! Compile terms of Ai in CFi which will be saved for next time step
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=-(NU/PR(N)) * KX2(I) * CTH(I,K,J,N)
     &            - (NU/PR(N)) * KZ2(K) * CTH(I,K,J,N)
          END DO
        END DO
      END DO

C End do number of passive scalars (N_TH)
      END DO

        CALL FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
        CALL FFT_XZ_TO_PHYSICAL(CU2,U2,0,NY+1)
        CALL FFT_XZ_TO_PHYSICAL(CU3,U3,0,NY+1)
! Transform THETA to physical space for computation of nonlinear terms
! Here pass the first location in memory of the array for scalar n
        DO N=1,N_TH
          CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N),0,NY+1)
        END DO
 

C Compute the nonlinear products in physical space, then transform
C back to Fourier space to compute the derivative.
C Here, we compute the horizontal derivatives of the nonlinear terms
C which will be treated with RKW3.  Nonlinear terms with vertical
C derivatives will be treated with Crank-Nicolson later
C Do terms one at a time to save on memory
C U1*U3
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=U3(I,K,J)*U1(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J) - CIKZ(K) * CS1(I,K,J) 
            CF3(I,K,J)=CF3(I,K,J) - CIKX(I) * CS1(I,K,J) 
          END DO
        END DO
      END DO

C U1*U1
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=U1(I,K,J)*U1(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J) - CIKX(I) * CS1(I,K,J) 
          END DO
        END DO
      END DO

C U3*U3
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=U3(I,K,J)*U3(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF3(I,K,J)=CF3(I,K,J) - CIKZ(K) * CS1(I,K,J) 
          END DO
        END DO
      END DO


C U1*U2
      DO J=2,NY
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=((DYF(J)*U1(I,K,J)
     &                +DYF(J-1)*U1(I,K,J-1))/(2.*DY(J))) 
     &                *U2(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=CF2(I,K,J) - CIKX(I) * CS1(I,K,J) 
          END DO
        END DO
      END DO

C U3*U2
      DO J=2,NY
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=((DYF(J)*U3(I,K,J)
     &                +DYF(J-1)*U3(I,K,J-1))/(2.*DY(J))) 
     &                *U2(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=CF2(I,K,J) - CIKZ(K) * CS1(I,K,J)
          END DO
        END DO
      END DO

! Add the vertical derivative term explicitly
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=
     &     (U1(I,K,J+1)*U2(I,K,J+1) + U1(I,K,J)*U2(I,K,J+1)
     &     - U1(I,K,J)*U2(I,K,J) - U1(I,K,J-1)*U2(I,K,J))/(2.d0*DYF(J))
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)

      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J) - CS1(I,K,J) 
         END DO
        END DO
      END DO

! Add the vertical derivative term explicitly
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=
     &     (U3(I,K,J+1)*U2(I,K,J+1) + U3(I,K,J)*U2(I,K,J+1)
     &     - U3(I,K,J)*U2(I,K,J) - U3(I,K,J-1)*U2(I,K,J))/(2.d0*DYF(J))
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)

      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF3(I,K,J)=CF3(I,K,J) - CS1(I,K,J) 
          END DO
        END DO
      END DO

! Add the vertical derivative term explicitly
      DO J=2,NY
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=
     &    (0.25d0*(U2(I,K,J)+U2(I,K,J+1))**2.d0
     &    -0.25d0*(U2(I,K,J)+U2(I,K,J-1))**2.d0)/DY(J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)

      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=CF2(I,K,J) - CS1(I,K,J)
          END DO
        END DO
      END DO


C -- At this point, we are done computing the nonlinear terms --

C Finally, Add CFi to CRi
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CR1(I,K,J) + TEMP5 * CF1(I,K,J)
            CR3(I,K,J)=CR3(I,K,J) + TEMP5 * CF3(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CR2(I,K,J)=CR2(I,K,J) + TEMP5 * CF2(I,K,J)
          END DO
        END DO
      END DO

C Convert RHS terms to physical space
      CALL FFT_XZ_TO_PHYSICAL(CR1,R1,0,NY+1)                 
      CALL FFT_XZ_TO_PHYSICAL(CR2,R2,2,NY)                 
      CALL FFT_XZ_TO_PHYSICAL(CR3,R3,0,NY+1)                 

C Compute the vertical viscous term in physical space and add to RHS
C This is the explicit part of the Crank-Nicolson term
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            R1(I,K,J)=R1(I,K,J)+TEMP1*
     &        (  ((U1(I,K,J+1) - U1(I,K,J)) / DY(J+1)  
     &           -(U1(I,K,J)   - U1(I,K,J-1)) / DY(J)) /DYF(J)  )
            R3(I,K,J)=R3(I,K,J)+TEMP1*
     &        (  ((U3(I,K,J+1) - U3(I,K,J)) / DY(J+1) 
     &           -(U3(I,K,J)   - U3(I,K,J-1)) / DY(J)) /DYF(J)  )
          END DO
        END DO
      END DO
      DO J=2,NY 
        DO K=0,NZM
          DO I=0,NXM
            R2(I,K,J)=R2(I,K,J)+TEMP1*
     &        (  ((U2(I,K,J+1) - U2(I,K,J))  / DYF(J) 
     &           -(U2(I,K,J)   - U2(I,K,J-1))/ DYF(J-1))/DY(J)  )
          END DO
        END DO
      END DO

C -- Here, we are done with computation of Velocity RHS, explicit terms --

C Now, build the explicit RHS terms for the passive scalar(s)

      DO N=1,N_TH
! Do for each scalar:

! Compute the nonlinear terms that are present in the explicit term A
! U1*TH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=TH(I,K,J,N)*U1(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CIKX(I) * CS1(I,K,J)
          END DO
        END DO
      END DO
! U3*TH 
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=TH(I,K,J,N)*U3(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CIKZ(K) * CS1(I,K,J)
          END DO
        END DO
      END DO
! We are done with the horizontal derivatives of the nonlinear terms
! Add the vertical derivative term explicitly
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=
     &     (TH(I,K,J+1,N)*U2(I,K,J+1) + TH(I,K,J,N)*U2(I,K,J+1)
     &    -TH(I,K,J,N)*U2(I,K,J)-TH(I,K,J-1,N)*U2(I,K,J))/(2.d0*DYF(J))
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CS1(I,K,J)
          END DO
        END DO
      END DO


! Add CFTH to the RHS vector CRTH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
            CRTH(I,K,J,N)=CRTH(I,K,J,N) + TEMP5 * CFTH(I,K,J,N)
          END DO
        END DO
      END DO
! Done with computation of RHS, explicit terms for the THETA equation
! Transform back to physical space

      CALL FFT_XZ_TO_PHYSICAL(CRTH(0,0,0,N),RTH(0,0,0,N),0,NY+1)    

! Compute the Explicit part of the Crank-Nicolson terms for the TH equation
! First, the vertical derivative viscous term
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZM
          DO I=0,NXM
            RTH(I,K,J,N)=RTH(I,K,J,N)+(TEMP1/PR(N))*(
     &            ((TH(I,K,J+1,N) - TH(I,K,J,N)) / DY(J+1)
     &            -(TH(I,K,J,N) - TH(I,K,J-1,N)) / DY(J)) / DYF(J) )
          END DO
        END DO
      END DO

C -- Now, timestep the passive scalar equation --
C      We solve the the passive scalar before the velocity so that
C      it is advected with the velocity from the previous R-K step
C      which we have already made divergence free 
 
! Solve the implicit equation for THETA
! Note that the system size is NY+1, but only 1..NY are used

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 
    
! Build implicit matrix
! Use quasi-second order interpolation for TH on GY points
      DO K=0,NZM
        DO J=JSTART_TH(N),JEND_TH(N)
          DO I=0,NXM
            MATL(I,J) = -(TEMP1/PR(N)) / (DY(J)*DYF(J))
            MATD(I,J) = 1. + (TEMP1/PR(N)) / (DY(J+1)*DYF(J))
     &           +(TEMP1/PR(N)) / (DY(J)*DYF(J))
            MATU(I,J)=-(TEMP1/PR(N)) / (DY(J+1)*DYF(J))
! Define RHS vector
            VEC(I,J)=RTH(I,K,J,N)
          END DO
        END DO

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN

             CALL APPLY_BC_TH_MPI(MATL,MATD,MATU,VEC,N)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
             CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          ELSE
! Else we are running in serial mode
             CALL APPLY_BC_TH_LOWER(MATL,MATD,MATU,VEC,N)
             CALL APPLY_BC_TH_UPPER(MATL,MATD,MATU,VEC,N)
             CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
! Apply the boundary conditions to our linear system
          END IF

        DO J=JSTART_TH(N),JEND_TH(N)
          DO I=0,NXM
            TH(I,K,J,N)=VEC(I,J)
          END DO
        END DO

! END do k
      END DO 

! End do number of passive scalars
        END DO

C Initialize the matrix to zeros to be used for implicit solves
C Note that the system size is NY+1, but only 1..NY are used

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 

C Build implicit matrix for U2
      DO K=0,NZM
        DO J=2,NY
          DO I=0,NXM
            MATL(I,J)= -TEMP1/(DYF(J-1)*DY(J))
            MATD(I,J)=1.+TEMP1/(DYF(J)*DY(J)) + TEMP1/(DYF(J-1)*DY(J)) 
            MATU(I,J)= -TEMP1/(DYF(J)*DY(J))
            VEC(I,J)=R2(I,K,J)
          END DO 
        END DO

        IF (USE_MPI) THEN

! First, apply the boundary conditions
          CALL APPLY_BC_U2_MPI(MATL,MATD,MATU,VEC)
! If we are using MPI, split the implicit solve into forward and
! backward sweeps for efficiency
          CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
        ELSE
C Else, we are running in serial mode
C Set the boundary conditions for U2
          CALL APPLY_BC_2_LOWER(MATL,MATD,MATU,VEC)
          CALL APPLY_BC_2_UPPER(MATL,MATD,MATU,VEC)

C Now, solve the tridiagonal system for U2(i,:,k)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
        END IF

        DO J=1,NY+1
          DO I=0,NXM
            U2(I,K,J)=VEC(I,J)
          END DO
        END DO
! End do k
      END DO 

C Solve for U1
C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 

C Build the implicit system of equations for U1 
      DO K=0,NZM
        DO J=JSTART,JEND
          DO I=0,NXM
            MATL(I,J)=-TEMP1/(DY(J)*DYF(J))
            MATD(I,J)=1.-TEMP1*(-1./(DY(J+1)*DYF(J))
     &         -1./(DY(J)*DYF(J))) 
            MATU(I,J)=-TEMP1/(DY(J+1)*DYF(J))
            VEC(I,J)=R1(I,K,J)
          END DO
        END DO

        IF (USE_MPI) THEN
! First, apply the boundary conditions
          CALL APPLY_BC_U1_MPI(MATL,MATD,MATU,VEC)
! If we are using MPI, split the implicit solve into forward and
! backward sweeps for efficiency
          CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
        ELSE
C Else, we are running in serial mode
C Set the boundary conditions for U1
          CALL APPLY_BC_1_LOWER(MATL,MATD,MATU,VEC)
          CALL APPLY_BC_1_UPPER(MATL,MATD,MATU,VEC)
C Now, solve the tridiagonal system for U1(:,k,:)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
        END IF

        DO J=JSTART-1,JEND+1
          DO I=0,NXM
            U1(I,K,J)=VEC(I,J)
          END DO
        END DO
! End do k
      END DO


! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 

C Solve for U3
C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)
C Build the implicit system of equations for U3
      DO K=0,NZM
        DO J=JSTART,JEND
          DO I=0,NXM
            MATL(I,J)=-TEMP1/(DY(J)*DYF(J))
            MATD(I,J)=1.-TEMP1*(-1./(DY(J+1)*DYF(J))
     &         -1./(DY(J)*DYF(J)))
            MATU(I,J)=-TEMP1/(DY(J+1)*DYF(J))
            VEC(I,J)=R3(I,K,J)
          END DO
        END DO

        IF (USE_MPI) THEN
! First, apply the boundary conditions
          CALL APPLY_BC_U3_MPI(MATL,MATD,MATU,VEC)
! If we are using MPI, split the implicit solve into forward and
! backward sweeps for efficiency
          CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
        ELSE
C Else, we are running in serial mode
C Set the boundary conditions for U3
          CALL APPLY_BC_3_LOWER(MATL,MATD,MATU,VEC)
          CALL APPLY_BC_3_UPPER(MATL,MATD,MATU,VEC)
C Now, solve the tridiagonal system for U3(i,:,k)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
        END IF

        DO J=JSTART-1,JEND+1
          DO I=0,NXM
            U3(I,K,J)=VEC(I,J)
          END DO
        END DO
! End do k
      END DO

C -- Done getting U1hat, U2hat, U3hat at new RK Step --
 
! Transform TH and U to Fourier Space 
      CALL FFT_XZ_TO_FOURIER(U1,CU1,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U2,CU2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U3,CU3,0,NY+1)
      DO N=1,N_TH
        CALL FFT_XZ_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N),0,NY+1)
      END DO

C Begin second step of the Fractional Step algorithm, making u divergence free
C The following subroutine projects Uhat onto divergence free space

      CALL REM_DIV_CHAN

C Now, phi is stored in CR1, use this to update the pressure field
C Note, here we divide by H_BAR since it was absorbed into PHI in REM_DIV
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CP(I,K,J)=CP(I,K,J)+CR1(I,K,J)/TEMP4
          END DO
        END DO
      END DO


      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
C Compute varphi, store in variable CR1.
C Solves for phi in computational space
C H_BAR has been absorbed into PHI, so we are solving for H_BAR*PHI

      INCLUDE 'header'
      INTEGER I,J,K
 
C First, Initialize the matrix components
      DO J=0,NY+1
        DO I=0,NKX
          MATL_C(I,J)=0.
          MATD_C(I,J)=1.
          MATU_C(I,J)=0.
          VEC_C(I,J)=(0.,0.)
        END DO
      END DO

C The 2d FFT of Ui should have been taken and stored in CUi
C Solving for phi amounts to solving a tridiagonal system
C First, construct the system to be solved
      DO K=0,TNKZ
        DO J=1,NY
          DO I=0,NKX
            MATL_C(I,J)=1./(DY(J)*DYF(J))
            MATD_C(I,J)=-KX2(I)-KZ2(K)
     &         -1./(DY(J+1)*DYF(J))-1./(DY(J)*DYF(J))
            MATU_C(I,J)=1./(DY(J+1)*DYF(J))
          END DO
        END DO

C Now, create the RHS vector
        DO J=1,NY         
          DO I=0,NKX
            VEC_C(I,J)=(CIKX(I)*CU1(I,K,J) 
     &            + (CU2(I,K,J+1)-CU2(I,K,J))/DYF(J) 
     &            + CIKZ(K)*CU3(I,K,J))
          END DO
        END DO

        IF (USE_MPI) THEN
C If we are using the MPI package...
          CALL APPLY_BC_REM_DIV_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
C First, do all forward sweeps
          CALL THOMAS_FORWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                 ,NY,NX/3)
C Now, do the backward sweeps
          CALL THOMAS_BACKWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                 ,NY,NX/3)
        ELSE
C Else we are running in serial mode
        DO I=0,NKX
          IF ((K.EQ.0).AND.(I.EQ.0)) THEN
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
C Otherwise the matrix will be singular
            MATL_C(I,1)=0. 
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)

            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          ELSE
C Use Dirichlet boundary conditions, dp/dz=0 at walls
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)

            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          END IF
        END DO
C Now solve the tridiagonal system for phi, store in CR1
        CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY,NKX)
        END IF


        DO J=1,NY
          DO I=0,NKX
            CR1(I,K,J)=VEC_C(I,J)
          END DO
        END DO

      END DO

C Now, Solve for CUi, the divergenceless velocity field
      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,J)=CU1(I,K,J)-CIKX(I)*CR1(I,K,J)
            CU3(I,K,J)=CU3(I,K,J)-CIKZ(K)*CR1(I,K,J)           
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CU2(I,K,J)=CU2(I,K,J)-(CR1(I,K,J)
     &             -CR1(I,K,J-1))/DY(J)
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE POISSON_P_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C We have CUi, need to compute CP.  Solve tridiagonal system exactly

      INCLUDE 'header'

      INTEGER I,J,K,N
      
	if (flavor.eq.'Basic') then
      WRITE(*,*) 'COMPUTING CP FROM CUI'
	end if

C First, construct the RHS vector, (dui/dxj)(duj/dxi) 
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CIKX(I)*CU1(I,K,J)
            CF2(I,K,J)=(CU2(I,K,J+1)-CU2(I,K,J))/DYF(J)
            CF3(I,K,J)=CIKZ(K)*CU3(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF3,F3,0,NY+1)
      
      DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
            F1(I,K,J)=F1(I,K,J)**2.
            F2(I,K,J)=F2(I,K,J)**2.
            F3(I,K,J)=F3(I,K,J)**2.
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(F2,CF2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(F3,CF3,0,NY+1)
      
C Now we have the diagonal terms, add to the rhs term
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=CF1(I,K,J)+CF2(I,K,J)+CF3(I,K,J)
          END DO
        END DO
      END DO

C Now get the first of the off-diagonal terms
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=(CU1(I,K,J+1)-CU1(I,K,J-1))/(2.*DYF(J))
            CF2(I,K,J)=CIKX(I)*0.5*(CU2(I,K,J)+CU2(I,K,J+1))
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)

C Compute product
      DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
            F1(I,K,J)=2.*F1(I,K,J)*F2(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)

C Add to RHS term
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX 
            CS1(I,K,J)=CS1(I,K,J)+CF1(I,K,J)
          END DO
        END DO
      END DO

C Now get the second of the off-diagonal terms
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=(CU3(I,K,J+1)-CU3(I,K,J-1))/(2.*DYF(J))
            CF2(I,K,J)=CIKZ(K)*0.5*(CU2(I,K,J)+CU2(I,K,J+1))
          END DO
        END DO
      END DO

C Convert to Physical space
      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)

C Compute product
      DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
            F1(I,K,J)=2.*F1(I,K,J)*F2(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)

C Add to RHS term
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=CS1(I,K,J)+CF1(I,K,J)
          END DO
        END DO
      END DO

C Now get the third of the off-diagonal terms
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CIKZ(K)*CU1(I,K,J)
            CF2(I,K,J)=CIKX(I)*CU3(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)
      
C Compute product
      DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
            F1(I,K,J)=2.*F1(I,K,J)*F2(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)

C Add to RHS term
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=CS1(I,K,J)+CF1(I,K,J)
          END DO
        END DO
      END DO     
     
C Finally, if the buoyancy force is active, then we need to add
C the contribution of the density to the pressure.  Note that the
C plane averaged density and the corresponding hydrostatic part of the
C pressure have been cancelled, so skip the 0,0 mode
      DO N=1,N_TH
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX  
            IF ((I.NE.0).or.(K.NE.0)) THEN
              CS1(I,K,J)=CS1(I,K,J)-RI_TAU(N)*
     &          (CTH(I,K,J+1,N)-CTH(I,K,J-1,N))/(GYF(J+1)-GYF(J-1))
            END IF
          END DO
        END DO
      END DO
      END DO

C Now, the RHS term should be stored in CS1     

C Construct the tridiagonal system in Fourier space to solve for CP
C First, zero the vectors
      DO J=0,NY+1
        DO I=0,NKX
          MATL_C(I,J)=0.d0
          MATD_C(I,J)=1.d0
          MATU_C(I,J)=0.d0
          VEC_C(I,J)=(0.,0.)
        END DO
      END DO

      DO K=0,TNKZ
        DO J=2,NYM
          DO I=0,NKX
            MATL_C(I,J)=1./(DY(J)*DYF(J))
            MATD_C(I,J)=-KX2(I)-KZ2(K)-1./(DY(J+1)*DYF(J))
     &                    -1./(DY(J)*DYF(J))
            MATU_C(I,J)=1./(DY(J+1)*DYF(J))   
            VEC_C(I,J)=-1.*CS1(I,K,J)
          END DO
        END DO
        
        IF (USE_MPI) THEN
          CALL APPLY_BC_POISSON_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
C First, do the forward sweeps
          CALL THOMAS_FORWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                 ,NY,NX/3)
C Now, do the backwared sweeps to put the solution in VEC_C
          CALL THOMAS_BACKWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                  ,NY,NX/3)
        ELSE
C Else we are running in serial mode
C Apply BCs
        DO I=0,NKX
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
          IF ((I.EQ.0).AND.(K.EQ.0)) THEN
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)
            MATD_C(I,NY)=-1.
            MATL_C(I,NY)=1.
            VEC_C(I,NY)=(0.,0.)
          ELSE
! Here, apply Neumann boundary conditions (dp/dz=0) at the walls
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)
            MATD_C(I,NY)=-1.
            MATL_C(I,NY)=1.
            VEC_C(I,NY)=(0.,0.)
          END IF
        END DO
C Now, solve for CP
        CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY,NKX)
        END IF

        DO J=1,NY
          DO I=0,NKX
            CP(I,K,J)=VEC_C(I,J)
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Initialize the scalar fields
C In this subroutine, you should initialize each scalar field for the
C particular problem of interest

      INCLUDE 'header'
      INTEGER I,J,K,N

      DO N=1,N_TH
        IF (CREATE_NEW_TH(N)) THEN

	IF (RHO_TYPE .EQ. 1) THEN
! Shear layer
       DO J=0,NY
         DO K=0,NZM
           DO I=0,NXM
             TH(I,K,J,N)=-TANH(GYF(J)*20.d0)
            END DO
          END DO
        END DO
   
	ELSE
! Linear Profile
	DO J=0,NY
         DO K=0,NZM
           DO I=0,NXM
             TH(I,K,J,N)=-GYF(J)
            END DO
          END DO
        END DO

	END IF
      
      CALL FFT_XZ_TO_FOURIER(TH(0,0,0,n),CTH(0,0,0,n),0,NY+1)

      END IF
      END DO


      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K
      REAL*8 RNUM1,RNUM2,RNUM3
      REAL CENTRE_VEL
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

C Initialize the random number generator
      CALL RANDOM_SEED(SIZE = K)
      Allocate (seed(1:K))
      seed(1:K)=10
      CALL RANDOM_SEED(PUT = seed)

C UBULK0 and KICK should be set in input.dat

C Set the laminar velocity profile in physical space
       write(*,*) 'UBULK0: ',UBULK0

C For Couette flow:
       DO J=0,NY
         DO K=0,NZM
           DO I=0,NXM
             U1(I,K,J)=gyf(J)
             U2(I,K,J)=0.
             U3(I,K,J)=0.
           END DO
         END DO
       END DO
	
	CENTRE_VEL = U1(1,1,NY/2)
	WRITE(*,*) 'centre_vel: ',CENTRE_VEL

C Zero the ghost cells
      IF (.NOT.USE_MPI) THEN
       DO K=0,NZM
         DO I=0,NXM
           U1(I,K,0)=0.
           U2(I,K,0)=0.
           U3(I,K,0)=0.
           U1(I,K,NY+1)=0.
           U2(I,K,NY+1)=0.
           U3(I,K,NY+1)=0.
         END DO
       END DO
      END IF
      
      CALL FFT_XZ_TO_FOURIER(U1,CU1,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U2,CU2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U3,CU3,0,NY+1)

      CALL SAVE_STATS_CHAN(.FALSE.)

      WRITE(*,*) 'KICK: ',KICK
      write(*,*) 'NKX,NY,TNKZ: ',NKX,NY,TNKZ
 
      DO I=1,NKX
        DO J=1,NY
          DO K=1,TNKZ
C Now, give the velocity field a random perturbation
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)

              CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK
              CU2(I,K,J)=CU2(I,K,J)+(RNUM2-0.5)*KICK
              CU3(I,K,J)=CU3(I,K,J)+(RNUM3-0.5)*KICK
           
          END DO
          IF (TNKZ.EQ.0) THEN
! Here, In the 2d case we want to add a kick to the mean in z
            K=0         
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)

              CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK
              CU2(I,K,J)=CU2(I,K,J)+(RNUM2-0.5)*KICK
              CU3(I,K,J)=CU3(I,K,J)+(RNUM3-0.5)*KICK
         
          END IF 

        END DO
      END DO

       IF (USE_MPI) THEN
         CALL GHOST_CHAN_MPI
       END IF

C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
      END IF

C Remove the divergence of the velocity field

      CALL REM_DIV_CHAN

      IF (USE_MPI) THEN
        CALL GHOST_CHAN_MPI
      END IF

C Get the pressure from the poisson equation
      CALL POISSON_P_CHAN

      IF (USE_MPI) THEN
        CALL GHOST_CHAN_MPI
      END IF

      CALL SAVE_STATS_CHAN(.FALSE.)

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE VIS_FLOW_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Convert to physical space and output the velocity and pressure 
C to be read for visualization

      INCLUDE 'header'
      CHARACTER*35 FNAME,FNAME_TH
      INTEGER      I, J, K, n,yind

      FNAME='diablo.vis'
      WRITE(6,*) 'Writing flow to ',FNAME
      OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")

      FNAME_TH='diablo_th.vis'
      WRITE(6,*) 'Writing flow to ',FNAME_TH
      OPEN(UNIT=20,FILE=FNAME_TH,STATUS="UNKNOWN",FORM="UNFORMATTED")

C Write out grid at GYF points
      open(unit=11,file='ygrid_out.txt'
     &       ,status='unknown',form='formatted')
      write(11,111) (GYF(j),j=2,NYM)
      close(11)
111   format(64(F16.8,' '))
      open(unit=12,file='xgrid_out.txt'
     &       ,status='unknown',form='formatted')
      write(12,112) (GX(i),i=0,NXM)
      close(12)
112   format(64(F16.8,' '))
      open(unit=13,file='zgrid_out.txt'
     &      ,status='unknown',form='formatted')
      write(13,113) (GZ(k),k=0,NZM)
      close(13)
113   format(64(F16.8,' '))

      IF (NUM_PER_DIR.EQ.3) THEN
C        Enter write statment here
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
C Convert to physical space
      call fft_xz_to_physical(CU1,U1,0,NY+1)
      call fft_xz_to_physical(CU2,U2,0,NY+1)
      call fft_xz_to_physical(CU3,U3,0,NY+1)
      call fft_xz_to_physical(CP,P,0,NY+1)
      do n=1,N_TH
        call fft_xz_to_physical(CTH(0,0,0,n),TH(0,0,0,n),0,NY+1)
      end do
C Interpolate the vertical velocity to GYF gridpoints
        WRITE(10) ((( 
     &    REAL(U1(I,K,J)),REAL(0.5*(U2(I,K,J)+U2(I,K,J+1)))
     &    ,REAL(U3(I,K,J)),REAL(U1(I,K,J)),
     &    REAL(0.5*(U2(I,K,J)+U2(I,K,J+1))),REAL(U3(I,K,J)),
     &    REAL(P(I,K,J))
     &    ,K=0,NZM),J=2,NYM),I=0,NXM)

        WRITE(20) ((((
     &     REAL(TH(I,K,J,n))
     &        ,n=1,N_TH),K=0,NZM),J=2,NYM), I=0,NXM)

C Output velocity field for input to LIC (Line integral convolution)
      open(61,file='lic_x.dat',form='formatted',status='unknown')
      write(61,*) NZ,NY
      open(62,file='lic_y.dat',form='formatted',status='unknown')
      write(62,*) NZ,NY
	
	yind = ny
	if (ny.ge.30) yind=30
      do j=1,NY
        write(61,161) ((real(U3(yind,k,j))),k=0,NZM)
        write(62,161) ((real(U2(yind,k,j))),k=0,NZM)
      end do
161   format(192(F8.3))
      close(61)
      close(62)

      ELSEIF (NUM_PER_DIR.EQ.1) THEN
C        Enter write statment here
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
C        Enter write statment here

      END IF
101   format(10(F16.8,' '))
      CLOSE(10)
      CLOSE(20)


C Compute the discrimenant for visualization of vortices
C Note, the velocity field will be destroyed by this calculation,
C so it is important that this is only done at the end of a simulation

C NOTE:  THIS SECTION NEEDS TO BE CHECKED, THERE MAY BE AN ERROR
C IN THE CALCULATION OF THE DISCRIMINANT

      IF ((NUM_PER_DIR.EQ.2)) THEN
C First, convert back to Fourier space
      call fft_xz_to_fourier(U1,CU1,0,NY+1)
      call fft_xz_to_fourier(U2,CU2,0,NY+1)
      call fft_xz_to_fourier(U3,CU3,0,NY+1)
C First, calculate the velocity gradient tensor at GYF points
      do j=2,NYM
        do k=0,TNKZ
          do i=0,NKX
            CA21(i,k,j)=CIKX(i)*0.5*(CU2(i,k,j+1)+CU2(i,k,j))
            CA31(i,k,j)=CIKX(i)*CU3(i,k,j)
            CA12(i,k,j)=(0.5*(CU1(i,k,j+1)+CU1(i,k,j))
     &                 - 0.5*(CU1(i,k,j)+CU1(i,k,j-1)))/DYF(j) 
            CA32(i,k,j)=(0.5*(CU3(i,k,j+1)+CU3(i,k,j))
     &                 - 0.5*(CU3(i,k,j)+CU3(i,k,j-1)))/DYF(j) 
            CA13(i,k,j)=CIKZ(k)*CU1(i,k,j)
            CA23(i,k,j)=CIKZ(k)*0.5*(CU2(i,k,j+1)+CU2(i,k,j))
C Now, the following will overwrite CUi
            CA11(i,k,j)=CIKX(i)*CU1(i,k,j)
            CA22(i,k,j)=(CU2(i,k,j+1)-CU2(i,k,j))/DYF(j)
            CA33(i,k,j)=CIKZ(k)*CU3(i,k,j)
          end do
        end do
      end do
C Transform to physicl space
      call fft_xz_to_physical(CA11,A11,0,NY+1)
      call fft_xz_to_physical(CA21,A21,0,NY+1)
      call fft_xz_to_physical(CA31,A31,0,NY+1)
      call fft_xz_to_physical(CA12,A12,0,NY+1)
      call fft_xz_to_physical(CA22,A22,0,NY+1)
      call fft_xz_to_physical(CA32,A32,0,NY+1)
      call fft_xz_to_physical(CA13,A13,0,NY+1)
      call fft_xz_to_physical(CA23,A23,0,NY+1)
      call fft_xz_to_physical(CA33,A33,0,NY+1)

C Defining S_ij=(A_ij+A_ji)/2, compute Strain_rate = S_ij S_ji
      do j=2,NYM
        do k=0,NZM
          do i=0,NXM
            Strain_rate(i,k,j)=0.5*((A12(i,k,j)+A21(i,k,j))**2.
     &                             +(A13(i,k,j)+A31(i,k,j))**2. 
     &                             +(A23(i,k,j)+A32(i,k,j))**2.)
     &                    +A11(i,k,j)**2.+A22(i,k,j)**2.+A33(i,k,j)**2. 
          end do
C Compute third invariant = -det(A)
C Overwrites A11
          do i=0,NXM
            Third_ivar(I,K,J) =
     &  - A11(I,K,J)*(A22(I,K,J)*A33(I,K,J) - A23(I,K,J)*A32(I,K,J))
     &  + A12(I,K,J)*(A21(I,K,J)*A33(I,K,J) - A23(I,K,J)*A31(I,K,J))
     &  - A13(I,K,J)*(A21(I,K,J)*A32(I,K,J) - A22(I,K,J)*A31(I,K,J))
          end do

C Defining Omega_ij=(A_ij-A_ji)/2, compute Enstrophy=-(Omega_ij Omega_ji). 
C Note that this loop overwrites A22.
          do i=0,NXM
            Enstrophy  (I,K,J) = 0.5*((A12(I,K,J) - A21(I,K,J))**2+
     &        (A13(I,K,J) - A31(I,K,J))**2+ 
     &        (A23(I,K,J) - A32(I,K,J))**2)
          end do

C Compute Second_ivar.
C Note that this loop overwrites A33.
          do i=0,NXM
            Second_ivar(I,K,J)=0.5*(Enstrophy(I,K,J)-Strain_rate(I,K,J))
          end do

C Compute Discriminant.
C Note that this loop overwrites A12.
          do i=0,NXM
            Discriminant(I,K,J) = 6.75*(Third_ivar(I,K,J))**2
     &       + (Second_ivar(I,K,J))**3
          end do
        end do 
      end do 
           
      OPEN(UNIT=20,FILE='inv.vis',STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(20) ((( 
     &    REAL(GX(I)),REAL(GYF(J))
     &    ,REAL(GZ(K)),REAL(Discriminant(i,k,j)),
     &    REAL(Enstrophy(i,k,j)),REAL(Strain_rate(i,k,j))
     &    ,K=0,NZM),J=2,NYM),I=0,NXM)
      CLOSE(20)

      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INPUT_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K,N

! Input parameters specific for Couette flow case

	U_BC_YMIN = 0
	U_BC_YMAX = 0
	V_BC_YMIN = 0
	W_BC_YMIN = 0
	V_BC_YMAX = 0
	W_BC_YMAX = 0

	U_BC_YMIN_C1 = -1.0
	U_BC_YMAX_C1 = 1.0
	V_BC_YMIN_C1 = 0.0
	W_BC_YMIN_C1 = 0.0
	V_BC_YMAX_C1 = 0.0
	W_BC_YMAX_C1 = 0.0
	
	
	IF (RHO_TYPE .EQ. 1) THEN
! Shear layer
  	  TH_BC_YMIN(1) = 1
	  TH_BC_YMAX(1) = 1
	  TH_BC_YMIN_C1(1) = 0.0
	  TH_BC_YMAX_C1(1) = 0.0

	ELSE
	  TH_BC_YMIN(1) = 0
	  TH_BC_YMAX(1) = 0
	  TH_BC_YMIN_C1(1) = 1.0
	  TH_BC_YMAX_C1(1) = -1.0
	END IF

	KICK = 0.01
	UBULK0 = 0.0

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_GRID_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*35 FNAME
      INTEGER I,J,K

         WRITE (6,*) 'Fourier in X'
         DO I=0,NX
           GX(I)=(I*LX)/NX
           DX(I)=LX/NX
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO
         WRITE (6,*) 'Fourier in Z'
         DO K=0,NZ
           GZ(K)=(K*LZ)/NZ
           DZ(K)=LZ/NZ
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO
         WRITE (6,*) 'Finite-difference in Y'

         write(*,*) 'USE_MPI: ',USE_MPI

         IF (USE_MPI) THEN
           FNAME='./ygrid'//MPI_IO_NUM//'.txt'
           write(*,*) 'FNAME: ',FNAME
           write(*,*) 'MPI_IO_NUM: ****',MPI_IO_NUM,'*****'
         ELSE
           FNAME='./ygrid.txt'
         END IF

         OPEN (30,file=FNAME,form='formatted',status='old')
         READ (30,*) NY_T
C Check to make sure that grid file is the correct dimensions
         IF (NY_T.ne.NY) THEN
           WRITE(6,*) 'NY, NY_T',NY,NY_T
           STOP 'Error: ygrid.txt wrong dimensions'
         END IF
         DO J=1,NY+1
           READ(30,*) GY(j)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO
         DO J=1,NY
           READ(30,*) GYF(j)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GYF(',J,') = ',GYF(J)
         END DO
         CLOSE(30)

         IF (USE_MPI) THEN
           CALL GHOST_GRID_MPI
         ELSE
C Define ghost cells
           GYF(0)=2.d0*GYF(1)-GYF(2)
           GYF(NY+1)=2.d0*GYF(NY)-GYF(NYM)
           GY(0)=2.d0*GY(1)-GY(2)
         END IF

C Define grid spacing
         DO J=1,NY+1
           DY(J)=(GYF(J)-GYF(J-1))
         END DO
         DO J=1,NY
           DYF(J)=(GY(J+1)-GY(J))
         END DO
         DYF(NY+1)=DYF(NY)

         RETURN 
         END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_CHAN(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      CHARACTER*35 FNAME
      LOGICAL FINAL
      integer i,j,k,n
      real*8 uc, ubulk

      WRITE(6,*) 'Saving flow statistics.'


C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
      END IF

      if (FINAL) then
! We are done with the simulation
! Close the NetCDF file
       CALL NETCDF_CLOSE_STATS_CHAN


        IF (USE_MPI) THEN
          FNAME='stats'//MPI_IO_NUM//'.txt'
        ELSE
          FNAME='stats.txt'
        END IF

        open(20,file=FNAME,form='formatted',status='unknown')
        do j=1,NY
          write(20,201) j,GYF(j),UBAR(j),VBAR(j),WBAR(j)
        end do
201     format(I3,',',F16.9,',',F16.9,',',F16.9,',',F16.9)
        do n=1,N_TH
        do j=1,NY
          write(20,202) j,GYF(j),THBAR(j,n)
        end do
        end do
202     format(I3,',',F16.9,',',F16.9)
        close(20)

      else

! We are in the middle of a run, compile statistics
! First get the number of samples taken so far
      write(*,*) 'TIME, DELTA_T: ',TIME, DELTA_T
      NSAMPLES=NSAMPLES+1
! Get the mean velocity
      do j=1,NY
        UBAR(j)=(1./float(NSAMPLES))*dble(CU1(0,0,j))
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*UBAR(j)
        VBAR(j)=(1./float(NSAMPLES))*dble(CU2(0,0,j))
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*VBAR(j)
        WBAR(j)=(1./float(NSAMPLES))*dble(CU3(0,0,j))
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*WBAR(j)
        do n=1,N_TH
          THBAR(j,n)=(1./float(NSAMPLES))*dble(CTH(0,0,j,n))
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*THBAR(j,n)
        end do
      end do

! Compute and write out the centerline velocity
      if (int(float(NY)/2.) .eq. float(NY)/2.) then
! IF NY is even
        uc=dble(CU1(0,0,int(float(NY)/2.))) 
      else
        uc=0.5*(dble(CU1(0,0,int(float(NY)/2.)-1))
     +         +dble(CU1(0,0,int(float(NY)/2.))))
      end if
      write(*,*) 'Centerline velocity = ', uc 
! Compute and write out bulk velocity
! Integrat the instantaneous mean profile numerically at GY points
      UBULK=0.
      do j=2,NY
        UBULK=UBULK+0.5*(dble(CU1(0,0,j)+CU1(0,0,j-1)))*DY(j)
      end do
      UBULK=UBULK/LY      
! Write out UBULK
      write(*,*) 'UBULK: ',UBULK

! Save CUi
      do k=0,TNKZ
        do i=0,NKX
          do j=0,NY+1
            CR1(i,k,j)=CU1(i,k,j)
            CR2(i,k,j)=CU2(i,k,j)
            CR3(i,k,j)=CU3(i,k,j)
          end do
        end do
      end do

! Convert to physical space
      call fft_xz_to_physical(CU1,U1,0,NY+1)
      call fft_xz_to_physical(CU2,U2,0,NY+1)
      call fft_xz_to_physical(CU3,U3,0,NY+1)

! Get the turbulent kinetic energy at each level 
      do j=1,NY
        urms(j)=0.
        vrms(j)=0.
        wrms(j)=0.
      do k=0,NZM
      do i=0,NXM 
        urms(j)=urms(j)+(abs(U1(i,k,j)-dble(CR1(0,0,j))))**2.
        vrms(j)=vrms(j)+(abs(U2(i,k,j)-dble(CR2(0,0,j))))**2.
        wrms(j)=wrms(j)+(abs(U3(i,k,j)-dble(CR3(0,0,j))))**2.
      end do
      end do
        urms(j)=sqrt(urms(j)/(float(NZ)*float(NX)))
        vrms(j)=sqrt(vrms(j)/(float(NZ)*float(NX)))
        wrms(j)=sqrt(wrms(j)/(float(NZ)*float(NX)))
      end do 
! Get the bulk rms value
      urms_b=0.
      do j=2,NY
        urms_b=urms_b+0.5*(urms(j)+urms(j-1))*DY(j)
      end do
      urms_b=urms_b/LY      
      vrms_b=0.
      do j=1,NY
        vrms_b=vrms_b+vrms(j)*DY(j)
      end do
      vrms_b=vrms_b/LY
      wrms_b=0.
      do j=2,NY
        wrms_b=wrms_b+0.5*(wrms(j)+wrms(j-1))*DY(j)
      end do
      wrms_b=wrms_b/LY
! Compute the Reynolds stress and mean velocity gradient
      do j=1,NY
        uv(j)=0. 
        uw(j)=0.
        wv(j)=0.
      do k=0,NZM
      do i=0,NXM
        uv(j)=uv(j)+(U1(i,k,j)-dble(CR1(0,0,j)))
     +    *(0.5*(U2(i,k,j)+U2(i,k,j+1))
     &    -0.5*dble(CR2(0,0,j)+CR2(0,0,j+1)))
        wv(j)=wv(j)+(U3(i,k,j)-dble(CR3(0,0,j)))
     +    *(0.5*(U2(i,k,j)+U2(i,k,j+1))
     &      -0.5*dble(CR2(0,0,j)+CR2(0,0,j+1)))
        uw(j)=uw(j)+(U1(i,k,j)-dble(CR1(0,0,j)))
     +    *(U3(i,k,j)-dble(CR3(0,0,j)))
      end do
      end do
        uv(j)=uv(j)/(float(NZ)*float(NX))
        uw(j)=uw(j)/(float(NZ)*float(NX))
        wv(j)=wv(j)/(float(NZ)*float(NX))
      end do
              
! Get the y-derivative of the mean velocity at GYF points
!      do j=1,NY
!        dudy(j)=dble(CR1(0,0,j+1)-CR1(0,0,j-1))/(2.*DYF(j))
!        dwdy(j)=dble(CR3(0,0,j+1)-CR3(0,0,j-1))/(2.*DYF(j))
!      end do
! Get the y-derivative of the mean velocity at GY points
      do j=1,NY
        dudy(j)=dble(CR1(0,0,j)-CR1(0,0,j-1))/(GYF(j)-GYF(j-1))
        dwdy(j)=dble(CR3(0,0,j)-CR3(0,0,j-1))/(GYF(j)-GYF(j-1))
      end do

! Calculate the mean square shear
      do j=1,NY
        shear(j)=0.d0
        do k=0,NZM
          do i=0,NXM
            shear(j)=shear(j)
     &            +((U1(i,k,j+1)-U1(i,k,j-1))/(2.d0*DYF(j)))**2.d0
     &            +((U3(i,k,j+1)-U3(i,k,j-1))/(2.d0*DYF(j)))**2.d0
          end do
        end do
        shear(j)=shear(j)/dble(NX*NZ)
      end do

! Write out the bulk rms velocity
      write(*,*) '<U_rms>: ',urms_b
      write(*,*) '<V_rms>: ',vrms_b
      write(*,*) '<W_rms>: ',wrms_b

! Get the rms vorticity
! First, get the x-component in fourier space
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=(CR3(i,k,j+1)-CR3(i,k,j-1))/(2.d0*DYF(j))
     &            -CIKZ(K)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
! Get the rms value
      do j=1,NY
      omega_x(j)=0.d0
      do k=1,NZM
      do i=1,NXM
        omega_x(j)=omega_x(j)+S1(i,k,j)**2.d0
      end do
      end do
      omega_x(j)=sqrt(omega_x(j)/(dble(NX-1)*dble(NZ-1)))
      end do

! Now, get the y-component in fourier space
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKZ(k)*CR1(i,k,j)-CIKX(i)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
! Get the rms value
      do j=1,NY
      omega_y(j)=0.d0
      do k=0,NZM
      do i=0,NXM
        omega_y(j)=omega_y(j)+S1(i,k,j)**2.d0
      end do
      end do
      omega_y(j)=sqrt(omega_y(j)/(dble(NX)*dble(NZ)))
      end do

! Now, get the y-component in fourier space
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
     &             -(CR1(i,k,j+1)-CR1(i,k,j-1))/(2.d0*DYF(j))
      end do
      end do
        CS1(0,0,j)=CS1(0,0,j)-dudy(j)
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
! Get the rms value
      do j=1,NY
      omega_z(j)=0.d0
      do k=0,NZM
      do i=0,NXM
        omega_z(j)=omega_z(j)+S1(i,k,j)**2.d0
      end do
      end do
      omega_z(j)=sqrt(omega_z(j)/(dble(NX)*dble(NZ)))
      end do


! Write out the mean statistics at each time
      IF (USE_MPI) THEN
        FNAME='mean'//MPI_IO_NUM//'.txt'
      ELSE
        FNAME='mean.txt'
      END IF
      open(40,file=FNAME,form='formatted',status='unknown')
      write(40,*) TIME_STEP,TIME,DELTA_T,UBULK
      do j=1,NY
        write(40,401) j,GYF(J),dble(CR1(0,0,j))
     +      ,0.5*dble((CR2(0,0,j-1)+CR2(0,0,j)))
     +      ,dble(CR3(0,0,j)),urms(j),vrms(j),wrms(j)
     +      ,uv(j),uw(j),wv(j),dudy(j),dwdy(j),dble(cp(0,0,j)),shear(j)
     &      ,omega_x(j),omega_y(j),omega_z(j)
      end do

401   format(I3,' ',17(F20.9,' '))

C Call netcdf
      CALL NETCDF_WRITE_STATS_CHAN

! Do over the number of passive scalars
      do n=1,N_TH

! Save CTH
      do k=0,TNKZ
        do i=0,NKX
          do j=0,NY+1
            CRTH(i,k,j,n)=CTH(i,k,j,n)
          end do
        end do
      end do

! Compute the scalar gradient and store in CRi
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
! Store gradients of TH(:,:,:,n) (if it is used) in CRi
          CR1(i,k,j)=CIKX(i)*CTH(i,k,j,n)
          CR2(i,k,j)=(CTH(i,k,j+1,n)-CTH(i,k,j-1,n))/(GYF(j+1)-GYF(j-1))
          CR3(i,k,j)=CIKZ(k)*CTH(i,k,j,n)
          end do
        end do
      end do
! Convert gradients to physical space
      CALL FFT_XZ_TO_PHYSICAL(CR1,R1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CR2,R2,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CR3,R3,0,NY+1)


! Convert to physical space
      call fft_xz_to_physical(CTH(0,0,0,n),TH(0,0,0,n),0,NY+1)
      do j=1,NY
        thrms(j,n)=0.
      do k=0,NZM
      do i=0,NXM
        thrms(j,n)=thrms(j,n)+(abs(TH(i,k,j,n)-dble(CRTH(0,0,j,n))))**2.
      end do
      end do
        thrms(j,n)=sqrt(thrms(j,n)/(float(NZ)*float(NX)))
      end do
! Compute the Reynolds stress and mean velocity gradient
      do j=1,NY
        thv(j,n)=0.
      do k=0,NZM
      do i=0,NXM
       thv(j,n)=thv(j,n)+(TH(i,k,j,n)-dble(CRTH(0,0,j,n)))
     +    *(0.5*(U2(i,k,j)+U2(i,k,j+1))
     &      -0.5*dble(CR2(0,0,j)+CR2(0,0,j+1)))
      end do
      end do
      thv(j,n)=thv(j,n)/(float(NZ)*float(NX))
      end do

! Get the y-derivative of the mean scalar at GYF points
      do j=1,NY
        dthdy(j,n)=dble(CRTH(0,0,j+1,n)-CRTH(0,0,j-1,n))/(2.*DYF(j))
      end do

! Compute the potential energy dissipation, grad(TH) \cdot grad(TH)
      do j=1,NY
        pe_diss(j,n)=0.d0
        do k=0,NZM
          do i=0,NXM
            pe_diss(j,n)=pe_diss(j,n)
     &          +R1(i,k,j)**2.d0+R2(i,k,j)**2.d0+R3(i,k,j)**2.d0
          end do
        end do
        pe_diss(j,n)=pe_diss(j,n)/dble(NX*NZ)
      end do
 
        if (n.eq.1) then
        open(88,file='movie_th.txt',status='unknown',form='formatted')
        do I=0,NXM
        do J=1,NY
          write(88,*) TH(I,0,J,n)
        end do
        end do

! This file will contain a single plane and is used in conjunction with
! the matlab script 'realtime_movie' to visualize data during
! simulation
        open (76,file='temp.txt',status='unknown',form='formatted')
        do J=1,NY
          write(76,*) gyf(j)
        end do
        do I=0,NXM
        do J=1,NY
          write(76,*) TH(I,0,J,n)
        end do
        end do
        close (76)
 !       CALL SYSTEM('mv temp.txt ../post_process/matlab/latest_slice.txt
 !    & ')
	CALL SYSTEM('mv temp.txt latest_slice.txt')
        end if

! Convert back to Fourier space
      call FFT_XZ_TO_FOURIER(TH(0,0,0,n),CTH(0,0,0,n),0,NY+1)

! End do over number of passive scalars, n
      end do


! Write out the mean statistics at each time
      IF (USE_MPI) THEN
        FNAME='mean_th'//MPI_IO_NUM//'.txt'
      ELSE
        FNAME='mean_th.txt'
      END IF
      open(41,file=FNAME,form='formatted',status='unknown')
      write(41,*) TIME_STEP,TIME,DELTA_T,UBULK
      do n=1,N_TH 
      do j=1,NY
        write(41,402) j,GYF(J),dble(CTH(0,0,j,n))
     +      ,dthdy(j,n),thrms(j,n),thv(j,n),pe_diss(j,n)
      end do
      end do

402   format(I3,' ',6(F20.9,' '))

      write(*,*) 'VERBOSITY: ',VERBOSITY
      if (VERBOSITY.gt.4) then 
      write(*,*) 'Outputting info for gnuplot...'
      open (unit=10, file="solution")
      do i=2,NXM
        do j=2,NYM
          write (10,*) i, j, U1(i,0,j)
        end do
        write (10,*) ""
      end do
      close (10)
      call system ('gnuplot <gnuplot.in') 
      end if

      if (MOVIE) then
! This file will contain a time history over a plane
        open(85,file='movie_u.txt',status='unknown',form='formatted')
        do I=0,NXM
        do J=1,NY
          write(85,*) U1(I,0,J)
        end do
        end do
        open(86,file='movie_v.txt',status='unknown',form='formatted')
        do I=0,NXM
        do J=1,NY
          write(86,*) U2(I,0,J)
        end do
        end do
        open(87,file='movie_w.txt',status='unknown',form='formatted')
        do I=0,NXM
        do J=1,NY
          write(87,*) U3(I,0,J)
        end do
        end do

        END IF

C Convert velocity back to Fourier space
      call fft_xz_to_fourier(U1,CU1,0,NY+1)
      call fft_xz_to_fourier(U2,CU2,0,NY+1)
      call fft_xz_to_fourier(U3,CU3,0,NY+1)

      end if

      write(*,*) 'done save_stats chan' 

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine tkebudget_chan
C----*|--.---------.---------.---------.---------.---------.---------.-|--

! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps
      character*35 FNAME
      include 'header'

      integer i,j,k

! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      do j=2,NYM
        epsilon(j)=0.
      end do
! Store du/dx in CS1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dx in CS1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*(CR2(i,k,j)+CR2(i,k,j+1))/2.0
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dy at GYF gridpoints, note remove mean
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        F1(i,k,j)=((U1(i,k,j+1)-CR1(0,0,j+1))
     &      -(U1(i,k,j-1)-CR1(0,0,j-1)))/(GY(j)+GY(j+1))
      end do
      end do
      end do
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dvdx*dudy
        epsilon(j)=epsilon(j)+(S1(i,k,j)*F1(i,k,j))
      end do
      end do
      end do
! Store dw/dx in CS1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dz at GYF gridpoints, note remove mean
! Store du/dz in CS1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKZ(k)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CF1,F1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dudz*dwdx
        epsilon(j)=epsilon(j)+S1(i,k,j)*F1(i,k,j)
      end do
      end do
      end do
! Compute dv/dy at GYF gridpoints, note remove mean
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        S1(i,k,j)=((U2(i,k,j+1)-CR2(0,0,j+1))-(U2(i,k,j)-CR2(0,0,j)))
     &            /GYF(j)
      end do
      end do
      end do
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute dw/dy at GYF gridpoints, note remove mean
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        S1(i,k,j)=((U3(i,k,j+1)-CR3(0,0,j+1))
     &      -(U3(i,k,j-1)-CR3(0,0,j-1)))/(GY(j)+GY(j+1))
      end do
      end do
      end do
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dz in CF1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKZ(k)*(CR2(i,k,j)+CR2(i,k,j+1))/2.0
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CF1,F1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dvdz*dwdy
        epsilon(j)=epsilon(j)+S1(i,k,j)*F1(i,k,j)
      end do
      end do
      end do
! Store dw/dz in CS1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKZ(k)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
      do j=2,NYM
        epsilon(j)=epsilon(j)/float(NX*NZ)
      end do


! Write out the bulk rms velocity
      write(*,*) '<U_rms>: ',urms_b


! Write out the mean statistics at each time
      IF (USE_MPI) THEN
        FNAME='tke'//MPI_IO_NUM//'.txt'
      ELSE
        FNAME='tke.txt'
      END IF
      open(45,file=FNAME,form='formatted',status='unknown')
      write(45,*) TIME_STEP,TIME,DELTA_T
      do j=2,NYM
        write(45,401) j,GYF(J),epsilon(j)
      end do
401   format(I3,' ',2(F20.9,' '))


      return 
      end
 
C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_1_LOWER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
C Dirichlet
        DO I=0,NXM
          MATL(I,0)=0. 
          MATD(I,0)=1.
          MATU(I,0)=0.                   
          VEC(I,0)=0.

          MATL(I,1)=0. 
          MATD(I,1)=1.
          MATU(I,1)=0.                   
          VEC(I,1)=U_BC_YMIN_C1 
        END DO

      RETURN 
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_1_UPPER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Top wall
C Dirichlet
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=U_BC_YMAX_C1
        END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|---
      SUBROUTINE APPLY_BC_2_LOWER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
C Dirichlet
        DO I=0,NXM
          MATL(I,1)=0.d0 
          MATD(I,1)=1.d0
          MATU(I,1)=0.d0                   
          VEC(I,1)=V_BC_YMIN_C1 

          MATL(I,2)=0.d0 
          MATD(I,2)=1.d0
          MATU(I,2)=0.d0                   
          VEC(I,2)=V_BC_YMIN_C1 
        END DO

C The following is only a placeholder, this row is used for U1 and U3
      DO I=0,NXM
        MATL(I,0) = 0.
        MATD(I,0) = 1.
        MATU(I,0) = 0.
        VEC(I,0) = 0.
      END DO

      RETURN
      END

 
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_2_UPPER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I
C Top wall
C Dirichlet
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=V_BC_YMAX_C1
          
          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=V_BC_YMAX_C1
        END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_3_LOWER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
C Dirichlet
        DO I=0,NXM
          MATL(I,0)=0. 
          MATD(I,0)=1.
          MATU(I,0)=0.                   
          VEC(I,0)=0.

          MATL(I,1)=0. 
          MATD(I,1)=1.
          MATU(I,1)=0.                   
          VEC(I,1)=W_BC_YMIN_C1
        END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_3_UPPER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Top wall
C Dirichlet
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=W_BC_YMAX_C1
        END DO

      RETURN 
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine APPLY_BC_TH_LOWER(MATL,MATD,MATU,VEC,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      include 'header'
      integer i,N

! Bottom Wall:
	IF (RHO_TYPE .EQ. 1) THEN
! Neumann
! NOTE: BC enforced at GY(2)
        do i=0,NXM
          MATL(i,1)=0.
          MATD(i,1)=-1.
          MATU(i,1)=1.
          VEC(i,1)=DY(2)*TH_BC_YMIN_C1(N)
        end do
        do i=0,NXM
          MATL(i,0)=0.
          MATD(i,0)=-1.
          MATU(i,0)=1.
          VEC(i,0)=DY(1)*TH_BC_YMIN_C1(N)
        end do

	ELSE
! Dirichlet
	DO I=0,NXM
          MATL(I,0)=0. 
          MATD(I,0)=1.
          MATU(I,0)=0.                   
          VEC(I,0)=0.

          MATL(I,1)=0. 
          MATD(I,1)=1.
          MATU(I,1)=0.                   
          VEC(I,1)=TH_BC_YMIN_C1(N)
        END DO	
	
	END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine APPLY_BC_TH_UPPER(MATL,MATD,MATU,VEC,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      include 'header'
      integer i,N
! Top wall
! Neumann
! NOTE: BC enforced at GY(NY)
	IF (RHO_TYPE .EQ. 1) THEN
        do i=0,NXM
          MATL(i,NY)=-1.
          MATD(i,NY)=1.
          MATU(i,NY)=0.
          VEC(i,NY)=DY(NY)*TH_BC_YMAX_C1(N)
        end do
        do i=0,NXM
          MATL(i,NY+1)=-1.
          MATD(i,NY+1)=1.
          MATU(i,NY+1)=0.
          VEC(i,NY+1)=DY(NY+1)*TH_BC_YMAX_C1(N)
        end do

	ELSE
! Dirichlet
	DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=TH_BC_YMAX_C1(N)
        END DO	

	END IF

      return
      end

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_LOWER
C----*|--.---------.---------.---------.---------.---------.---------.-|--
C This subroutine is called after initializing the flow
C It sets the appropriate boundary conditions including ghost cell values
C  on the velocity field in Fourier space
      INCLUDE 'header'
      INTEGER I,K      

C Now, apply the boundary conditions depending on the type specified 
C U1
C Dirichlet 
C Start with zero
         DO K=0,TNKZ
           DO I=0,NKX
             CU1(I,K,1)=0.d0
           END DO
         END DO
C Now, set only the mean
         CU1(0,0,1)=U_BC_YMIN_C1
C Ghost cell not used
         CU1(0,0,0)=0.d0


C Dirichlet
C U3
C Start with zero
         DO K=0,TNKZ
           DO I=0,NKX
             CU3(I,K,1)=0.d0
           END DO
         END DO
C Now, set only the mean
         CU3(0,0,1)=W_BC_YMIN_C1
C Ghost cell not used
         CU3(0,0,0)=0.d0

C Dirichlet
C U2
C Set the vertical velocity at GYF(1) (halfway between GY(2) and GY(1))
         DO K=0,TNKZ
           DO I=0,NKX
             CU2(I,K,1)=2.d0*V_BC_YMIN_C1-CU2(I,K,2) 
           END DO
         END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_UPPER
C----*|--.---------.---------.---------.---------.---------.---------.-|--
C This subroutine is called after initializing the flow
C It sets the appropriate boundary conditions including ghost cell values
C  on the velocity field in Fourier space
      INCLUDE 'header'
      INTEGER I,K      

! Now, apply boundary conditions to the top of the domain
C Dirichlet
C U1 
C Start with zero
         DO K=0,TNKZ
           DO I=0,NKX
             CU1(I,K,NY)=0.d0
           END DO
         END DO
C Now, set only the mean
         CU1(0,0,NY)=U_BC_YMAX_C1
C Ghost cell not used
         CU1(0,0,NY+1)=0.d0

C Dirichlet
C U3
C Start with zero
         DO K=0,TNKZ
           DO I=0,NKX
             CU3(I,K,NY)=0.d0
           END DO
         END DO
C Now, set only the mean
         CU3(0,0,NY)=W_BC_YMAX_C1
C Ghost cell not used
         CU3(0,0,NY+1)=0.d0

C Dirichlet
C U2
C Set the vertical velocity at GYF(NY) (halfway between GY(NY) and GY(NY+1))
         DO K=0,TNKZ
           DO I=0,NKX
             CU2(I,K,NY+1)=2.d0*V_BC_YMAX_C1-CU2(I,K,NY)
           END DO
         END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine filter_chan
C----*|--.---------.---------.---------.---------.---------.---------.-|--

C This subroutine applies a filter to the highest wavenumbers
C It should be applied to the scalars in Fourier space
C The filter used is a sharpened raised cosine filter in the horizontal
C and a fourth order implicit compact filter in the vertical, with the
C parameter alpha determining the width of the vertical filtering window

      include 'header'

      integer I,J,K,js,je,N

! Variables for horizontal filtering
      real*8 sigma(0:NKX,0:TNKZ),sigma0

! Variables for vertical filtering
      real*8 alpha
      parameter (alpha=0.0d0)
! Parameters for a larger stencil filter
      real*8 f_a,f_b,f_c

      js=0
      je=NY+1

C Set the filtering constants for the horizontal direction
      DO i=0,NKX
       DO k=0,TNKZ
        sigma0=0.5d0*(1.d0+
     &       cos(sqrt((KX(i)*LX*1.d0/float(NX))**2.d0
     &            +(KZ(k)*LZ*1.d0/float(NZ))**2.d0)))
! Apply a sharpened raised cosine filter
        sigma(i,k)=sigma0**4.d0*(35.d0-84.d0*sigma0
     &        +70.d0*sigma0**2.d0-20.d0*sigma0**3.d0)
       END DO
      END DO

      DO N=1,N_TH
C Do the spectral filtering in the horizontal
        DO K=0,TNKZ
          DO I=0,NKX
            DO J=js+1,je-1
              CTH(I,K,J,N)=CTH(I,K,J,N)*sigma(i,k)
            END DO
          END DO
        END DO
      END DO
C Set the filtering constants
      f_a=(1.d0/8.d0)*(5.d0+6.d0*alpha)
      f_b=0.5d0*(1.d0+2.d0*alpha)
      f_c=(-1.d0/8.d0)*(1.d0-2.d0*alpha)


      DO N=1,N_TH
C First, zero the tridiagonal matrix components
      DO I=0,NKX
        DO J=0,NY+1
          MATD_C(I,J)=1.d0
          MATL_C(I,J)=0.d0
          MATU_C(I,J)=0.d0
          VEC_C(I,J)=0.d0
        END DO
      END DO


C Filter the passive scalar, TH in the vertical direction
      DO K=1,TNKZ
        DO I=1,NKX
C Construct the centered difference terms
          DO J=2,NY-1
            MATL_C(I,J)=alpha
            MATD_C(I,J)=1.d0
            MATU_C(I,J)=alpha
            VEC_C(I,J)=f_a*CTH(I,K,J,N)
     &                +(f_b/2.d0)*(CTH(I,K,J+1,N)+CTH(I,K,J-1,N))
     &                +(f_c/2.d0)*(CTH(I,K,J+2,N)+CTH(I,K,J-2,N))
          END DO
C Now, construct the equations for the boundary nodes
          J=1
            MATL_C(I,J)=0.d0
            MATD_C(I,J)=1.d0
            MATU_C(I,J)=0.d0
            VEC_C(I,J)=CTH(I,K,J,N)
          J=NY
            MATL_C(I,J)=0.d0
            MATD_C(I,J)=1.d0
            MATU_C(I,J)=0.d0
            VEC_C(I,J)=CTH(I,K,J,N)
         END DO
C Now, solve the tridiagonal system
         CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY,NKX)
         DO I=1,NKX
           DO J=js+1,je-1
             CTH(I,K,J,N)=VEC_C(I,J)
           END DO
         END DO
C END DO K  
       END DO

C END DO N 
       END DO
       return
       end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_REAL(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NX, NY
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY), G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO J=NY-1,0,-1
        DO I=0,NX
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|    
      SUBROUTINE THOMAS_COMPLEX(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
C The RHS vector and solution is complex
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NY, NX
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY)
      COMPLEX*16 G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO I=0,NX
        DO J=NY-1,0,-1
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END


