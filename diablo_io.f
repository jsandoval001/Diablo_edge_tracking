C This file contains subroutines for inputting and outputting data in
C Diablo as well as all subroutines called directly from diablo.f
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INITIALIZE
        INCLUDE 'header'

        REAL    VERSION, CURRENT_VERSION
        logical RESET_TIME
        INTEGER I, J, K, N
        INTEGER ISTAT, ENERGY_BISECTION_COUNTER

        ! Dummy variabls

        logical AUX_LOG
        integer AUX_INT
        real*8  AUX_REAL1, AUX_REAL2, AUX_REAL3

        REAL*8 DE

        WRITE(*,*) 'here'

        CURRENTFILE = INPATH(:LIP)//'input.dat'
        WRITE(*,*) 'file ',CURRENTFILE
        OPEN (11,file=CURRENTFILE,form='formatted',status='old')      
        WRITE(*,*) 'INPUT.DAT OPENED'

C Read input file.
C   (Note - if you change the following section of code, update the
C    CURRENT_VERSION number to make obsolete previous input files!)

        CURRENT_VERSION=2.0
        READ(11,*) 
        READ(11,*) 
        READ(11,*) 
        READ(11,*) 
        READ(11,*) FLAVOR,   VERSION
        IF (VERSION .NE. CURRENT_VERSION) THEN
          PRINT *, 'Wrong input data format.'
          STOP
        END IF
        READ(11,*) 
        READ(11,*) USE_MPI,    LES,    READ_HDF5,    SAVE_HDF5
        READ(11,*) 
        READ(11,*) RE, LX, LY, LZ, TIME_LIMIT
        READ(11,*) 
        READ(11,*) NUM_PER_DIR, CREATE_NEW_FLOW, RHO_TYPE, TIME_AVERAGE
        READ(11,*) 
        READ(11,*) N_TIME_STEPS, DELTA_T, RESET_TIME, VARIABLE_DT 
     &, CFL , UPDATE_DT
        READ(11,*) 
        READ(11,*) VERBOSITY, SAVE_FLOW_INT, SAVE_STATS_INT 
     & , MOVIE, INIT_E , T0
        READ(11,*) 

        ! Read in the parameters for the N_TH scalars
        DO N=1,MAX(N_TH,1)
          READ(11,*) 
          READ(11,*) AUX_LOG
          
          IF (N_TH >= 1) THEN
            CREATE_NEW_TH(N) = AUX_LOG
          END IF
          
          READ(11,*) 
          READ(11,*) AUX_LOG, AUX_INT
          
          IF (N_TH >= 1) THEN
            FILTER_TH(N)  = AUX_LOG
            FILTER_INT(N) = AUX_INT
          END IF
          
          READ(11,*) 
          READ(11,*) AUX_REAL1, AUX_REAL2, AUX_REAL3

          IF (N_TH >= 1) THEN
            RI_TAU(N)   = AUX_REAL1
            PR(N)       = AUX_REAL2
            REACTION(N) = AUX_REAL3
          END IF          

        END DO !N=1,MAX(N_TH,1)

        READ(11,*) 
        READ(11,*) 
        READ(11,*) 
        READ(11,*) 
        READ(11,*) 
        READ(11,*) WALL_STOKES , PRESSURE_STOKES
        READ(11,*) 
        READ(11,*) ET_MODE
        READ(11,*) 
        READ(11,*) RESUME_RESCALING_SIM , RESUME_BISECTION_SIM , 
     &             SAVE_ALL_BSCTN_FFS
        READ(11,*) 
        READ(11,*) ET_RESCALING , ET_BISECTION, SET_LAMBDA_INIT
        READ(11,*) 
        READ(11,*) LATEST_LAMBDA_LAM, LATEST_LAMBDA_TURB
        READ(11,*) 
        READ(11,*) MIN_T_ES_EVAL, AVG_WINDOW_ES_EVAL
        READ(11,*) 
        READ(11,*) DELTA_LAMBDA_THRS, COMP_WINDOW_SPS_EVAL , 
     &             DELTA_TRAJECTORIES_THRS, ENERGY_BISECT_THRS
        READ(11,*) 
        READ(11,*) KIN_E_LAM_THRS, KIN_E_TURB_THRS
        READ(11,*) 
        READ(11,*) LATEST_ALPHA_E_LAM, LATEST_ALPHA_E_TURB
        READ(11,*) 
        READ(11,*) TOL_ET_ITERATION , T_ET_COMPARISON
        READ(11,*) 
        READ(11,*) APPLY_SYMMETRIES
        
        ! After this line there is more info in input.dat, but it's
        ! about stratified flow, so will not be considered in this 
        ! version of the code.

        CLOSE(11)

        ! check if the flags for the edge trackin mode are consistent
        CALL CHECK_ET_INPUTS_ET( ISTAT )

        ! if the check comes back with istat/=0, then i terminate the
        ! simulation 
        IF ( ISTAT /= 0 ) CALL END_RUN_MPI_INIT()

        IF ( RESUME_BISECTION_SIM ) THEN

          ! read all the indexes and turbulent/laminar lambda values 
          ! from the log files
          CALL INIT_RESUME_BISECTION()

          ! LATEST_LAMBDA_LAM and LATEST_LAMBDA_TURB read from 
          ! lambda_shifts.ind
          LAMBDA = 0.5D0 * ( LATEST_LAMBDA_LAM + LATEST_LAMBDA_TURB )

          IF (RANK_G == 0) THEN
            ! I reset save.flow.ind
            OPEN( unit = 20, file = 'save_flow.ind',  
     &           status = 'replace', action = 'write')
          
            WRITE(20, '(I0)') 0
            CLOSE(20)

          END IF

        ELSE ! SET EVERYTHING AS A NEW BISECTION RUN

          ! shift counters and vars for bisection mode
          
          ! this counters (ENERGYCOUNT_LAM and ENERGYCOUNT_TURB) put me 
          ! in the 1st position of the arrays that append the KE every
          ! time step
          ENERGYCOUNT_LAM       = 1
          ENERGYCOUNT_TURB      = 1
          
          LAMBDA_LAM_FLAG       = .FALSE.
          LAMBDA_TURB_FLAG      = .FALSE.
  
          LAMBDA_SHIFTS_COUNT      = 0
          LATEST_LAMBDA_LAM_COUNT  = 0
          LATEST_LAMBDA_TURB_COUNT = 0
          
          SP_SHIFTS_COUNT       = 0
          
          IF (.NOT. SET_LAMBDA_INIT ) THEN
  
            LATEST_LAMBDA_LAM     = 0.D0
            LATEST_LAMBDA_TURB    = 1.D0
  
          END IF
  
          LAMBDA = 0.5D0 * ( LATEST_LAMBDA_LAM + LATEST_LAMBDA_TURB )

        END IF

        NU = 1.0d0/RE
        
        ! # of times the energy has been shifted
        SHIFTS_COUNT = 0
        ALPHA_E = 0.5D0 * ( LATEST_ALPHA_E_LAM + LATEST_ALPHA_E_TURB )


C If we are using MPI, then Initialize the MPI Variables
        WRITE(*,*) 'ABOUT TO INITIALIZE MPI'
        IF (USE_MPI) THEN
          CALL INIT_MPI
        ELSE
          RANK_G = 0
        END IF
        WRITE(*,*) 'MPI INITIALIZED'

        IF (RANK_G.EQ.0) THEN
          WRITE(*,*) 'RHO_TYPE: ',RHO_TYPE
        END IF
  
C Initialize channel package.

        CALL INPUT_CHAN
        CALL CREATE_GRID_CHAN
      
        IF (USE_MPI) THEN
          CALL INIT_CHAN_MPI
        ELSE 
          CALL INIT_CHAN
        END IF
     

C Initialize grid
        IF (FLAVOR.NE.'Ensemble') THEN
        IF (RANK_G.EQ.0) THEN

        WRITE(6,*) 'Note that this code is distributed under the '
          WRITE(6,*) 'GNU General Public License.'
          WRITE(6,*) 'No warranty is expressed or implied.'
          WRITE(6,*)
          write(*,*) 'Flavor: ',FLAVOR
          WRITE(6,*) 'Grid size: NX =',NX,', NY =',NY,', NZ =',NZ,'.'
          WRITE(6,*) 'RE: ',RE,', NU: ',NU
         
          DO N=1,N_TH
            WRITE(6,*) 'Scalar number: ',N
            WRITE(6,*) '  Richardson number: ',RI_TAU(N)
            WRITE(6,*) '  Prandlt number: ',PR(N)
          END DO
      
        END IF
        END IF
      

        NXM=NX-1
        NYM=NY-1
        NZM=NZ-1

C Initialize storage arrays.
        DO K=0,NZ+1
          DO I=0,NX+1 
            DO J=0,NY+1
              U1(I,K,J)=0.d0
              U3(I,K,J)=0.d0
              U2(I,K,J)=0.d0
              P (I,K,J)=0.d0
              R1(I,K,J)=0.d0
              R2(I,K,J)=0.d0
              R3(I,K,J)=0.d0
              F1(I,K,J)=0.d0
              F2(I,K,J)=0.d0
              F3(I,K,J)=0.d0
              BU1(I,K,J)=0.d0
              DBU1(I,K,J) = 0.D0
              BU1F(I,K,J)=0.d0
              DO N = 1,N_TH
              BTH(I,K,J,N)=0.d0
              END DO
            END DO
          END DO
        END DO

C Initialize FFT package (includes defining the wavenumber vectors).
        CALL INIT_FFT

C Initialize energy storing arrays and counters     
        KINS(:)         = 0.D0
        KINS_n(:)       = 0.D0 ! KE array from previous iteration
        ENERGIES(:)     = 0.D0
        DISSIPATIONS(:) = 0.D0
        ENERGYCOUNT     = 1
        DISSCOUNT       = 1

        KINS_LAM (:)    = 0.D0  
        KINS_TURB(:)    = 0.D0  

C Initialize RKW3 parameters.
        H_BAR(1)=DELTA_T*(8.0d0/15.0d0)
        H_BAR(2)=DELTA_T*(2.0d0/15.0d0)
        H_BAR(3)=DELTA_T*(5.0d0/15.0d0)
        BETA_BAR(1)=1.0d0
        BETA_BAR(2)=25.0d0/8.0d0
        BETA_BAR(3)=9.0d0/4.0d0
        ZETA_BAR(1)=0.0d0
        ZETA_BAR(2)=-17.0d0/8.0d0
        ZETA_BAR(3)=-5.0d0/4.0d0
        H_ALPHA(1)=DELTA_T*(0.D0)
        H_ALPHA(2)=DELTA_T*(8.D0/15.D0)
        H_ALPHA(3)=DELTA_T*(2.D0/3.D0)
        DISSIPATION = 0.D0
        NSAMPLES = 0

C Initialize values for reading of scalars
        NUM_READ_TH=0
        DO N=1,N_TH
          IF (CREATE_NEW_TH(N)) THEN
            NUM_READ_TH=NUM_READ_TH 
          ELSE
            NUM_READ_TH=NUM_READ_TH + 1
            READ_TH_INDEX(NUM_READ_TH)=N
          END IF
        END DO
   
        CALL CREATE_BACKGROUND_TH_CHAN
        CALL CREATE_TH_CHAN 
  
C initialise adjoint system (to be done)

!   CALL INIT_ADJOINT

C Create flow.
C Initialize flow.
        IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
          PREVIOUS_TIME_STEP = 0
          TIME_STEP          = 0
          TIME               = T0 ! From input.dat
        END IF

        TIME_TIME_STEP = TIME
        ! I CREATE THE STOKES BOUNDARY LAYER PROFILE
        ! If CREATE_NEW_FLOW is fals, it reads the flow
        ! from start.h5 anyway
        CALL CREATE_BACKGROUND_FLOW_CHAN
  
        IF ( CREATE_NEW_FLOW ) THEN
    
          CALL CREATE_FLOW_CHAN
     
          IF (FLAVOR.NE.'Ensemble') THEN
          IF (RANK_G.EQ.0) THEN
            write(*,*) 'A new flowfield has been created'
          END IF
          END IF
          
          IF (FLAVOR.EQ.'Basic') THEN
            CALL SAVE_STATS(.FALSE.)
!          CALL SAVE_FLOW(.FALSE.)
          END IF
      
        ELSE ! Reads the flow variables from start.h5
      
          IF (RANK_G.EQ.0) THEN
            write(*,*) 'Reading flow...'
          END IF
        
          IF ( .NOT. ET_BISECTION ) THEN
            
            CALL READ_FLOW
            
            ! Time counters
            TIME_STEP = 0
            TIME      = T0FILE

            IF (RANK_G.EQ.0) THEN
              write(*,*) 'Done reading flow'
            END IF

          END IF
!        CALL PERTURB
               
        END IF ! IF (CREATE_NEW_FLOW)
 
C Initialize flow.
        IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN

          PREVIOUS_TIME_STEP=0
          TIME_STEP=0
          TIME=T0

        END IF


        IF ( ET_MODE .AND. ET_BISECTION ) THEN
          ! I call read flow with lambda = 0 (laminar) and then 
          ! lambda = 1 (turbulent) to know the initial turbulent and 
          ! laminar energies

          KINS(:)         = 0.D0
          KINS_n(:)       = 0.D0 
          ENERGIES(:)     = 0.D0
          DISSIPATIONS(:) = 0.D0
          ENERGYCOUNT     = 1
          DISSCOUNT       = 1
          
          ! laminar flow field initial energy
          CALL READ_FLOW_BISECTION( 0.D0 )        

          CALL GET_ENERGY(.FALSE.)
          ENERGY_LAM_INI = ENERGIES(1)

          print *, 'ENERGY_LAM_INI = ', ENERGY_LAM_INI

          KINS(:)         = 0.D0
          KINS_n(:)       = 0.D0 
          ENERGIES(:)     = 0.D0
          DISSIPATIONS(:) = 0.D0
          ENERGYCOUNT     = 1
          DISSCOUNT       = 1
          
          ! turbulent flow field initial energy
          CALL READ_FLOW_BISECTION( 1.D0 )        

          CALL GET_ENERGY(.FALSE.)
          ENERGY_TURB_INI = ENERGIES(1)

          print *, 'ENERGY_TURB_INI = ', ENERGY_TURB_INI

          ! I restore the arrays as they were meant to be
          KINS(:)         = 0.D0
          KINS_n(:)       = 0.D0 
          ENERGIES(:)     = 0.D0
          DISSIPATIONS(:) = 0.D0
          ENERGYCOUNT     = 1
          DISSCOUNT       = 1

          ENERGY_AVG = 0.5D0 * ( ENERGY_LAM_INI + ENERGY_TURB_INI )
          DE         = ABS( ENERGY_TURB_INI - ENERGY_LAM_INI )

          ! I read the flow from start_turb.h5 and start_lam.h5 and
          ! using LAMBDA
          CALL READ_FLOW_BISECTION( LAMBDA )
          CALL GET_KE_ONLY( ENERGY_BISECTION )

          IF ( RANK_G == 0 ) THEN
            PRINT *, 'ENERGY( LAMBDA ) INI = ', ENERGY_BISECTION
          END IF

          ! I iterate over LAMBDA until the energy is close the midpoint
          ! between laminar and turbulent
          ENERGY_BISECTION_COUNTER = 0

          CALL WAIT()

          DO WHILE( ABS( ENERGY_BISECTION - ENERGY_AVG ) > 0.05D0 * DE
     &      .AND. ENERGY_BISECTION_COUNTER < 30  
     &      .AND. .NOT. RESUME_BISECTION_SIM )

            ! Slightly modifying LAMBDA to converge to the midpoint
            IF ( ENERGY_BISECTION > ENERGY_AVG ) LAMBDA = LAMBDA - 0.05
            IF ( ENERGY_BISECTION < ENERGY_AVG ) LAMBDA = LAMBDA + 0.05

            ! read and update the velocity field
            CALL READ_FLOW_BISECTION( LAMBDA )
            ! compute the KE without writing the files
            CALL GET_KE_ONLY( ENERGY_BISECTION )

            ENERGY_BISECTION_COUNTER = ENERGY_BISECTION_COUNTER + 1

          END DO
          
          CALL WAIT()

          IF ( RANK_G == 0 ) THEN

            PRINT *, 'ENERGY( LAMBDA ) CONVERGED = ', ENERGY_BISECTION
            PRINT *, 'AFTER = ', ENERGY_BISECTION_COUNTER, ' ITERS '

          END IF

          ! Time counters
          TIME_STEP = 0
          TIME      = T0FILE

        END IF

        ! if RESUME_BISECTION_SIM is TRUE, then KINS_n is read from 
        ! lambda_shifts.ind. Otherwise, is set to 0.
        !if ( .NOT. RESUME_BISECTION_SIM ) KINS_n = 0.D0

        ! If I'm in Edge Tracking mode, it doesn't matter the INIT_E 
        ! given in input.dat. What the model is going to consider is 
        ! the one that start.h5 has got.

        CALL GET_ENERGY(.FALSE.)
        INIT_E0 = ENERGIES(1)
      
        IF ( ET_MODE .AND. ET_RESCALING ) THEN
          ! Rescale the energy to start the first ET iteration
          INIT_E = ALPHA_E * INIT_E0
          
          IF ( RANK_G == 0 ) THEN
            WRITE(*,*) 'INIT_E0: ', INIT_E0          
          END IF

          CALL GET_ENERGY(.TRUE.)

          IF (RANK_G.EQ.0) THEN
            WRITE(*,*) 'INIT_E after rescaling : ', INIT_E
          END IF

        ELSE

          INIT_E = INIT_E0

        END IF

        CALL WAIT()
        ! I initialise the shift log files
        IF ( RANK_G == 0 ) THEN

          ! I write to energy_shift_history.dat, the alpha value and 
          ! the new energy after the last energy shift (write header
          ! set to .TRUE.)
          IF ( ET_RESCALING ) CALL WRITE_ENERGY_SHIFT_HISTORY(.TRUE.)

          ! I write to lambda_shift_history.dat the starting point
          ! information (count and physical time) and the lambda shift
          ! information (count, current lambda value and energy values
          ! associated). I also write the header of the file (that's the
          ! argument that's set to .TRUE.). Same for sp_shift_history

          IF ( ET_BISECTION .AND. .NOT. RESUME_BISECTION_SIM ) THEN
            CALL WRITE_LAMBDA_SHIFT_HISTORY(.TRUE.,.FALSE.)
            CALL WRITE_SP_SHIFT_HISTORY    (.TRUE.)
            CALL WRITE_BISECTION_LOG()
          END IF

        END IF

        CALL WAIT()

        CALL SAVE_STATS(.FALSE.)
        CALL POISSON_P_CHAN

        RETURN
   
      END SUBROUTINE INITIALIZE

      SUBROUTINE APPLY_SYMMETRY()

        INCLUDE 'header'

        INTEGER I, J, K, N

        REAL*8 G1(0:NXM,0:NZM,0:NY+1), G2(0:NXM,0:NZM,0:NY+1),
     &         G3(0:NXM,0:NZM,0:NY+1), GP(0:NXM,0:NZM,0:NY+1)
     
        CALL FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
        CALL FFT_XZ_TO_PHYSICAL(CU2,U2,0,NY+1)
        CALL FFT_XZ_TO_PHYSICAL(CU3,U3,0,NY+1)
        CALL FFT_XZ_TO_PHYSICAL(CP,P,0,NY+1)
      
        DO J=0,NY+1
        DO K=0,NZM
        DO I=0,NXM
          IF (I.LT.(NX/2)) THEN
            G1(I,K,J) = (U1(I,K,J)+U1(I+NX/2,NZM-K,J))/2.D0
            G2(I,K,J) = (U2(I,K,J)+U2(I+NX/2,NZM-K,J))/2.D0
            G3(I,K,J) = (U3(I,K,J)-U3(I+NX/2,NZM-K,J))/2.D0
            GP(I,K,J) = (P(I,K,J)+P(I+NX/2,NZM-K,J))/2.D0
!          END IF
          ELSE
            G1(I,K,J) = (U1(I,K,J)+U1(I-NX/2,NZM-K,J))/2.D0
            G2(I,K,J) = (U2(I,K,J)+U2(I-NX/2,NZM-K,J))/2.D0
            G3(I,K,J) = (U3(I,K,J)-U3(I-NX/2,NZM-K,J))/2.D0
            GP(I,K,J) = (P(I,K,J)+P(I-NX/2,NZM-K,J))/2.D0
          END IF
        END DO
        END DO
        END DO

        DO J=0,NY+1
        DO K=0,NZM
        DO I=0,NXM
          U1(I,K,J) = G1(I,K,J)
          U2(I,K,J) = G2(I,K,J)
          U3(I,K,J) = G3(I,K,J)
          P(I,K,J)  = GP(I,K,J)
        END DO
        END DO
        END DO

        CALL FFT_XZ_TO_FOURIER(U1,CU1,0,NY+1)
        CALL FFT_XZ_TO_FOURIER(U2,CU2,0,NY+1)
        CALL FFT_XZ_TO_FOURIER(U3,CU3,0,NY+1)
        CALL FFT_XZ_TO_FOURIER(P,CP,0,NY+1)      

        IF (USE_MPI) THEN
          CALL GHOST_CHAN_MPI
        END IF

        ! Apply Boundary conditions to velocity field
        IF (USE_MPI) THEN
          CALL APPLY_BC_VEL_MPI
        ELSE
          CALL APPLY_BC_VEL_LOWER
          CALL APPLY_BC_VEL_UPPER
        END IF

      END SUBROUTINE APPLY_SYMMETRY

      !=================================================================
      ! 
      ! ###### #    # ###### #####   ####  #   #                        
      ! #      ##   # #      #    # #    #  # #                         
      ! #####  # #  # #####  #    # #        #                          
      ! #      #  # # #      #####  #  ###   #                          
      ! #      #   ## #      #   #  #    #   #                          
      ! ###### #    # ###### #    #  ####    #                          
      !                                                                 
      !                                                                 
      ! ###### #    #   ##   #      #    #   ##   ##### #  ####  #    # 
      ! #      #    #  #  #  #      #    #  #  #    #   # #    # ##   # 
      ! #####  #    # #    # #      #    # #    #   #   # #    # # #  # 
      ! #      #    # ###### #      #    # ######   #   # #    # #  # # 
      ! #       #  #  #    # #      #    # #    #   #   # #    # #   ## 
      ! ######   ##   #    # ######  ####  #    #   #   #  ####  #    # 
      !
      !=================================================================
                                                                 

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ET_ES_EVAL( ENERGY_RESCALING_FLAG )
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
   
        ! --------------------------------------------------------------
        ! DESCRIPTION:
        ! ------------
        ! 
        ! This routine:
        !
        ! 1. Calculates the average KIN_AVG of the energy for the time
        !     window defined by [ TIME - AVG_WINDOW_ES_EVAL , TIME ] 
        ! 
        ! 2. Compares with the laminar and turbulend user defined energy
        !    thresholds KIN_E_LAM_THRS, KIN_E_TURB_THRS
        ! 
        ! 3. If KIN_AVG < KIN_E_LAM_THRS, then it sets
        ! 
        !        ENERGY_RESCALING_FLAG = -1
        !        SHIFTS_COUNT = SHIFTS_COUNT + 1
        !    
        !    And renames
        ! 
        !        kin.txt ---> kin_SHIFTS_COUNT_lam.txt   
        !  
        !    If KIN_AVG > KIN_E_TURB_THRS, then it sets
        !
        !        ENERGY_RESCALING_FLAG = 1
        !        SHIFTS_COUNT = SHIFTS_COUNT + 1
        !
        !    And renames
        ! 
        !        kin.txt ---> kin_SHIFTS_COUNT_turb.txt   
        !
        !    So in the main folder, we will see outfile files like
        !    
        !    kin_000_turb.txt
        !    kin_001_lam.txt
        !    kin_002_lam.txt
        !    kin_003_turb.txt
        !        .
        !        .
        !        .
        !    kin_SHIFTS_COUNT_turb.txt
        !
        !    If KIN_AVG is between the thresholds, then 
        !    ENERGY_RESCALING_FLAG = 0 (it actually remain  unchanged 
        !    as it is set to 0 within diablo.f)
        !
        ! --------------------------------------------------------------

        INCLUDE 'header'

        ! I/O
        INTEGER :: ENERGY_RESCALING_FLAG

        ! Local variables
        CHARACTER*512 :: KIN_F_NAME
        INTEGER :: N_TSTEPS_AVG
        REAL*8  :: KIN_AVG
        INTEGER :: SUM_IDX

        ! If the actual simulation time (TIME-T0FILE) is less than the
        ! necessary time for averaging, I just return. This should be
        ! controlled in diablo.f anyway, but just to be sure
        IF ( TIME-AVG_WINDOW_ES_EVAL < T0FILE ) RETURN

        ! N_TSTEPS_AVG: amount of energy values contained in an 
        ! averaging window
        N_TSTEPS_AVG = NINT( ( AVG_WINDOW_ES_EVAL / DELTA_T )
     &                       / SAVE_STATS_INT )

        ! Another caution condition. If the amount of energy values 
        ! saved are less than the one I am meant to average, then
        ! the routine just returns
        IF ( ENERGYCOUNT < N_TSTEPS_AVG ) RETURN

        KIN_AVG = 0.D0

        DO SUM_IDX = ENERGYCOUNT - N_TSTEPS_AVG, ENERGYCOUNT-1
          KIN_AVG = KIN_AVG + KINS(SUM_IDX)
        END DO

        KIN_AVG = KIN_AVG / REAL( N_TSTEPS_AVG )
   
        ! TURBULENCE
    
        IF ( KIN_AVG > KIN_E_TURB_THRS ) THEN
    
          ENERGY_RESCALING_FLAG = 1
          LATEST_ALPHA_E_TURB   = ALPHA_E
    
          ! Rename kin
    
          CALL WAIT()

          IF ( RANK_G == 0 ) THEN
              
            WRITE( KIN_F_NAME , '(a, i3.3, a)' ) 
     &           "kin_", SHIFTS_COUNT , "_turb.txt"
    
            CALL RENAME( 'kin.txt', TRIM( KIN_F_NAME ) )
           
            ! Create kin again as an empty file
    
            OPEN ( unit=10, file='kin.txt', status='replace', 
     &             action='write')
            CLOSE(unit=10)
  
          END IF
          
          SHIFTS_COUNT = SHIFTS_COUNT + 1
    
          ! ALPHA_E UPDATE
          ALPHA_E = 0.5D0*(LATEST_ALPHA_E_LAM + LATEST_ALPHA_E_TURB)

          ! Synchronise all the procs before reinitialising
          CALL WAIT()
    
        END IF
   
        ! LAMINARISATION
    
        IF ( KIN_AVG < KIN_E_LAM_THRS ) THEN
    
          ENERGY_RESCALING_FLAG = -1
          LATEST_ALPHA_E_LAM = ALPHA_E
    
          ! Rename kin
          
          CALL WAIT()

          IF ( RANK_G == 0 ) THEN
    
            WRITE( KIN_F_NAME, '(a, i3.3, a)')  
     &      'kin_' , SHIFTS_COUNT , '_lam.txt'
        
            CALL RENAME( 'kin.txt', TRIM( KIN_F_NAME ) )
  
            ! Create kin again as an empty file
        
            OPEN ( unit=10, file='kin.txt', status='replace',  
     &             action='write')
            CLOSE(unit=10)

          END IF
  
          SHIFTS_COUNT = SHIFTS_COUNT + 1
    
          ! ALPHA_E UPDATE
          ALPHA_E = 0.5D0 * ( LATEST_ALPHA_E_LAM + LATEST_ALPHA_E_TURB )

          ! Synchronise all the procs before reinitialising
          CALL WAIT()
    
        END IF
   
      END SUBROUTINE ET_ES_EVAL


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ET_LAMBDA_EVAL( LAMBDA_RESCALING_FLAG )
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
   
      ! This subroutine checks if a lambda shift has to be done. If so,
      ! then it updates lamda, lambda counters, and renames kin.txt as 
      ! kin_lam or kin_turb, depending on the evaluated regime of the 
      ! current lambda value.

        INCLUDE 'header'

        ! I/O
        INTEGER :: LAMBDA_RESCALING_FLAG !  1 --> DOWNSCALING (TURB)
                                         ! -1 --> UPSCALING (LAM)
        ! Local variables
        CHARACTER*512 :: KIN_F_NAME
        INTEGER :: N_TSTEPS_AVG
        REAL*8  :: KIN_AVG
        INTEGER :: SUM_IDX

        CALL WAIT()

        ! Initialise before evaluation
        LAMBDA_RESCALING_FLAG = 0        

        ! If the simulation time (TIME-T0FILE) is less than the
        ! necessary time for averaging, it just returns. This should
        ! be controlled in diablo.f anyway, but just to be sure
        IF ( TIME - T0FILE < AVG_WINDOW_ES_EVAL ) RETURN

        ! # of elements to average to evaluate a lambda shift
        N_TSTEPS_AVG = NINT( ( AVG_WINDOW_ES_EVAL / DELTA_T )
     &                       / SAVE_STATS_INT )

        ! Another caution condition. If the amount of energy values 
        ! saved are less than the one I am meant to average, then
        ! the routine just returns
        IF ( ENERGYCOUNT < N_TSTEPS_AVG ) RETURN

        KIN_AVG = 0.D0

        DO SUM_IDX = ENERGYCOUNT - N_TSTEPS_AVG, ENERGYCOUNT-1
          KIN_AVG = KIN_AVG + KINS( SUM_IDX ) 
        END DO

        KIN_AVG = KIN_AVG / REAL( N_TSTEPS_AVG )
   
        ! TURBULENCE
        IF ( KIN_AVG > KIN_E_TURB_THRS ) THEN
          
          LAMBDA_RESCALING_FLAG = 1

          LAMBDA_TURB_FLAG    = .TRUE.
          LATEST_LAMBDA_TURB  = LAMBDA

          ! Here I save the amount of energy values stored the last
          ! time a turbulent case was run
          ENERGYCOUNT_TURB = ENERGYCOUNT-1

          ! I update the latest turbulent KE time series
          KINS_TURB(:) = KINS(:)
    
          ! Rename kin.txt and move all the solutions to the folder
          ! "latest_turbulent_solutions"

          !CALL WAIT()

          IF ( RANK_G == 0 ) THEN
              
            ! For example: kin.txt --> kin_002_005_turb.txt
            WRITE( KIN_F_NAME, '(a, i3.3, a, i3.3, a)') 
     &           "kin_" , SP_SHIFTS_COUNT , "_" , LAMBDA_SHIFTS_COUNT , 
     &           "_turb.txt"
    
            CALL RENAME( 'kin.txt', TRIM( KIN_F_NAME ) )
           
            ! Create kin again as an empty file
            OPEN(unit=10,file='kin.txt',status='replace',action='write')
            CLOSE(unit=10)
            
            ! Remove al the previous files from the destination folder
            call execute_command_line(
     &      'rm -f latest_turbulent_solutions/*.h5')
             
            ! Move all the solutions to latest_turbulent_solutions
            call execute_command_line ( 
     &      'mv *out.h5 latest_turbulent_solutions/' )

          END IF

          ! lambda shift counter update
          LATEST_LAMBDA_TURB_COUNT = LAMBDA_SHIFTS_COUNT 
          LAMBDA_SHIFTS_COUNT      = LAMBDA_SHIFTS_COUNT + 1
    
          ! lambda update
          LAMBDA = 0.5D0 * ( LATEST_LAMBDA_LAM + LATEST_LAMBDA_TURB )
          
          ! write the bisection log files sp_shifts.ind and 
          ! lambda_shifts.ind
          IF ( RANK_G == 0 ) CALL WRITE_BISECTION_LOG()

          ! Synchronise all the procs before reinitialising
          !CALL WAIT()

        END IF ! KIN_AVG > KIN_E_TURB_THRS
   
        ! --------------------------------------------------------------
        !
        ! LAMINARISATION
    
        IF ( KIN_AVG < KIN_E_LAM_THRS ) THEN
    
          LAMBDA_RESCALING_FLAG = -1

          LAMBDA_LAM_FLAG    = .TRUE.
          LATEST_LAMBDA_LAM  = LAMBDA
    
          ! Here I save the amount of energy values stored the last
          ! time an laminar case was run
          ENERGYCOUNT_LAM = ENERGYCOUNT-1

          ! I update the latest turbulent KE time series
          KINS_LAM(:) = KINS(:)

          ! Rename kin.txt and move all the solutions to a the folder
          ! latest_laminar_solutions

          !CALL WAIT()

          IF ( RANK_G == 0 ) THEN
    
            ! For example: kin.txt --> kin_003_002_lam.txt
            WRITE(KIN_F_NAME, '(a, i3.3, a, i3.3, a)') 
     &           "kin_" , SP_SHIFTS_COUNT , "_" , LAMBDA_SHIFTS_COUNT , 
     &           "_lam.txt"
        
            CALL RENAME( 'kin.txt', TRIM( KIN_F_NAME ) )
  
            ! Create kin again as an empty file
            OPEN(unit=10,file='kin.txt',status='replace',action='write')
            CLOSE(unit=10)

            ! Remove al the previous files from the destination folder
            call execute_command_line (
     &      'rm -f latest_laminar_solutions/*.h5' )
             
            ! Move all the solutions to latest_laminar_solutions
            call execute_command_line( 
     &      'mv *out.h5 latest_laminar_solutions/')

          END IF

          LATEST_LAMBDA_LAM_COUNT  = LAMBDA_SHIFTS_COUNT 
          LAMBDA_SHIFTS_COUNT      = LAMBDA_SHIFTS_COUNT + 1
    
          ! LAMBDA UPDATE
          LAMBDA = 0.5D0 * ( LATEST_LAMBDA_LAM + LATEST_LAMBDA_TURB )

          ! write the bisection log files sp_shifts.ind and 
          ! lambda_shifts.ind
          IF( RANK_G == 0 ) CALL WRITE_BISECTION_LOG()
    
        END IF ! KIN_AVG < KIN_E_LAM_THRS
   
        CALL WAIT()

      END SUBROUTINE ET_LAMBDA_EVAL


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ET_SP_EVAL( SP_SHIFT_FLAG )
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
   
      ! This subroutine evals whether to rescale lambda or not

        INCLUDE 'header'

        ! I/O
        LOGICAL :: SP_SHIFT_FLAG

        ! Local variables
        CHARACTER*512 :: KIN_F_NAME
        INTEGER :: N_TSTEPS_AVG
        REAL*8  :: KIN_AVG
        INTEGER :: SUM_IDX
        INTEGER :: N_TSTEPS_COMP

        REAL*8  :: DELTA_LAMBDA
        REAL*8  :: MAX_DELTA_KIN

        CALL WAIT ()

        SP_SHIFT_FLAG = .FALSE.
        DELTA_LAMBDA  = ABS ( LATEST_LAMBDA_TURB - LATEST_LAMBDA_LAM )

        ! If I haven't generated two new series of solution to find a
        ! a new starting point (SP), or the difference between the 
        ! turbulent lambda and the laminar one is larger than the  
        ! threshold, then it just returns
        IF (.NOT. LAMBDA_TURB_FLAG .OR. .NOT. LAMBDA_LAM_FLAG .OR. 
     &      DELTA_LAMBDA > DELTA_LAMBDA_THRS ) RETURN
          
        ! If the simulation time (TIME-T0FILE) is less than the
        ! necessary time for averaging, I just return. This should
        ! be controlled in diablo.f anyway, but just to be sure
        IF ( TIME - T0FILE < AVG_WINDOW_ES_EVAL ) RETURN

        N_TSTEPS_COMP = NINT( ( COMP_WINDOW_SPS_EVAL / DELTA_T )
     &                       / SAVE_STATS_INT )

        ! Another caution condition. If the amount of energy values 
        ! saved are less than the one I am meant to compare, then
        ! the routine just returns
        IF ( ENERGYCOUNT_TURB < N_TSTEPS_COMP .OR.
     &       ENERGYCOUNT_LAM  < N_TSTEPS_COMP ) RETURN

        MAX_DELTA_KIN = MAXVAL(ABS( KINS_LAM ( 2:N_TSTEPS_COMP-1 ) -
     &                              KINS_TURB( 2:N_TSTEPS_COMP-1 ) ) )

        ! This is for us to check if the trajectories has gotten close 
        ! enough
        IF ( MAX_DELTA_KIN > DELTA_TRAJECTORIES_THRS ) RETURN

        SP_SHIFT_FLAG = .TRUE.

        CALL WAIT ()

      END SUBROUTINE ET_SP_EVAL

      !=================================================================
      !
      ! #####  ###### # #    # # ##### #   ##       
      ! #    # #      # ##   # #   #   #  #  #      
      ! #    # #####  # # #  # #   #   # #    #     
      ! #####  #      # #  # # #   #   # ######     
      ! #   #  #      # #   ## #   #   # #    #     
      ! #    # ###### # #    # #   #   # #    #     
      !                                             
      !                                             
      ! #      #  ####    ##   ##### #  ####  #    #
      ! #      # #       #  #    #   # #    # ##   #
      ! #      #  ####  #    #   #   # #    # # #  #
      ! #      #      # ######   #   # #    # #  # #
      ! #      # #    # #    #   #   # #    # #   ##
      ! ###### #  ####  #    #   #   #  ####  #    #
      !
      !=================================================================


      SUBROUTINE ET_ENERGY_REINITIALIZATION( ENERGY_RESCALING_FLAG )

        ! --------------------------------------------------------------
        ! DESCRIPTION:
        ! ------------
        ! 
        ! This routine performs the reinitilisation of the flow to start
        ! a new iteration in the RESCALING mode. It also makes a copy
        ! of the *out.h5 files of the last iteration
        ! 
        ! 
        ! The steps it does are:
        !
        ! 1. Set all the flow variables and time series to zero
        !
        ! 2. Reads again the flow variables from start.h5
        !
        ! 3. Redefines INIT_E as INIT_E = ALPHA_E * INIT_E0, where:
        !
        !         INIT_E0: Energy read from start.h5 when INITIALIZE()
        !                  was called at the beginning of the simualtion
        ! 
        !         ALPHAE : Rescaling parameter calculated in ET_ES_EVAL
        !                  as  
        !
        !                  ALPHA_E = 0.5D0*( LATEST_ALPHA_E_LAM + 
        !                                    LATEST_ALPHA_E_TURB  )
        !
        ! 4. Rescale the flow variables with this new INIT_E by calling
        !        
        !       CALL GET_ENERGY(.TRUE.)
        !
        !    where the '.TRUE.' flag indicates that the flow has to be 
        !    rescaled using the value of INIT_E      
        !
        ! 5. WRITING LOG FILES
        !    
        !    it calls WRITE_ENERGY_SHIFT_HISTORY(.FALSE.), where the 
        !    .FALSE. flag is to tell the subroutine not to write the
        !    header (which is written only in INITIALIZE). This routine
        !    updates the file energy_shift_history.dat with the
        !    REINITIALISED values of
        !
        !    SHIFTS_COUNT , ALPHA_E , INIT_E  
        !
        !
        ! 6. FILE MANAGEMENT
        !    
        !    6.1 It removes all the files from the folder
        !        'last_ET_it_flow_field' and copy all the *out.h5 files
        !        into it
        ! 
        !    6.2 It renames the file 000_out.h5 either as 
        !        000_out_lam.h5 or 000_out_turb.h5 depending on the 
        !        value of ENERGY_RESCALING_FLAG (either -1 or 1).
        !
        !    6.3 Depending on the value ENERGY_RESCALING_FLAG, it  
        !        either cleans up the folder latest_laminar_solutions 
        !        or latest_turbulent_solutions and moves all the *out.h5
        !        files into it
        !                          
        ! --------------------------------------------------------------
        
        INCLUDE 'header'

        CHARACTER*512 FNAME
        CHARACTER*512 INIT_SOL_FNAME
        CHARACTER*512 LATEST_INITIAL_SOL_STR

        INTEGER , INTENT(IN) :: ENERGY_RESCALING_FLAG
        REAL    VERSION, CURRENT_VERSION
        logical RESET_TIME
        INTEGER I, J, K, N

        WRITE(*,*) 'Reinitialising Flow for Edge Tracking'
  
        CALL WAIT()

        ! Initialize storage arrays: set all of them to zero.
        DO K=0,NZ+1
          DO I=0,NX+1 
            DO J=0,NY+1
              U1(I,K,J)=0.d0
              U3(I,K,J)=0.d0
              U2(I,K,J)=0.d0
              P (I,K,J)=0.d0
              R1(I,K,J)=0.d0
              R2(I,K,J)=0.d0
              R3(I,K,J)=0.d0
              F1(I,K,J)=0.d0
              F2(I,K,J)=0.d0
              F3(I,K,J)=0.d0
              BU1(I,K,J)=0.d0
              DBU1(I,K,J) = 0.D0
              BU1F(I,K,J)=0.d0
              DO N = 1,N_TH
              BTH(I,K,J,N)=0.d0
              END DO
            END DO
          END DO
        END DO

        ! Initialize FFT package (includes defining the wavenumber 
        ! vectors)
        CALL INIT_FFT
        
        ! Store the kinetic energy time series of the latest
        ! ET iteration
        KINS_n(:) = KINS(:)
        
        ! Reinitialisation of the energy time series
        KINS(:)         = 0.D0
        ENERGIES(:)     = 0.D0
        DISSIPATIONS(:) = 0.D0
        
        ! Amount of energy values stored
        ENERGYCOUNT     = 1
        DISSCOUNT       = 1
      
        ! Initialize RKW3 parameters.
        H_BAR(1)=DELTA_T*(8.0d0/15.0d0)
        H_BAR(2)=DELTA_T*(2.0d0/15.0d0)
        H_BAR(3)=DELTA_T*(5.0d0/15.0d0)
        BETA_BAR(1)=1.0d0
        BETA_BAR(2)=25.0d0/8.0d0
        BETA_BAR(3)=9.0d0/4.0d0
        ZETA_BAR(1)=0.0d0
        ZETA_BAR(2)=-17.0d0/8.0d0
        ZETA_BAR(3)=-5.0d0/4.0d0
        H_ALPHA(1)=DELTA_T*(0.D0)
        H_ALPHA(2)=DELTA_T*(8.D0/15.D0)
        H_ALPHA(3)=DELTA_T*(2.D0/3.D0)
        DISSIPATION = 0.D0
        NSAMPLES = 0

        ! Create flow.
        ! Initialize flow.
        
        IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
          PREVIOUS_TIME_STEP=0
          TIME_STEP=0
          TIME=T0
        END IF
        
        IF (CREATE_NEW_FLOW) THEN
        
           CALL CREATE_FLOW_CHAN
      
           IF (FLAVOR.NE.'Ensemble') THEN
           IF (RANK_G.EQ.0) THEN
              write(*,*) 'A new flowfield has been created'
           END IF
           END IF
           
           IF (FLAVOR.EQ.'Basic') THEN
              CALL SAVE_STATS(.FALSE.)
           END IF
        
        ELSE
        
           IF (RANK_G.EQ.0) THEN
              write(*,*) 'Reading flow...'
           END IF
           
           CALL READ_FLOW
           
           IF (RANK_G.EQ.0) THEN
              write(*,*) 'Done reading flow'
           END IF
        
        END IF
        
        !Initialize flow.
        IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
          PREVIOUS_TIME_STEP=0
          TIME_STEP=0
          TIME=T0
        END IF
  
        ! Time counters
        TIME_STEP = 0
        TIME      = T0FILE

        ! This variable is used to create a new background
        ! flow with the right phase
        
        TIME_TIME_STEP = TIME

        CALL CREATE_BACKGROUND_FLOW_CHAN
  
        ! Iniitial energy reinitialisation. This ALPHA_E comes
        ! recalculated from ET_ES_EVAL, where it was defined if it was
        ! either scaled down or up. INIT_E0 is obtained from the 
        ! start.h5 file read in READ_FLOW
        INIT_E = ALPHA_E * INIT_E0
  
        IF ( RANK_G == 0 ) THEN
          WRITE(*,*) 'INIT_E0               : ', INIT_E0
          WRITE(*,*) 'ALPHA_E Reinitialised : ', ALPHA_E          
          WRITE(*,*) 'INIT_E  Reinitialised : ', INIT_E
        END IF
  
        ! I rescale the flow with the NEW initial energy value INIT_E
        CALL GET_ENERGY(.TRUE.)
        ! Save the new energy value
        CALL GET_ENERGY(.FALSE.)
        CALL SAVE_STATS(.FALSE.)
  
        IF ( RANK_G == 0 ) CALL WRITE_ENERGY_SHIFT_HISTORY(.FALSE.)
        
        ! All the processors wait until proc 0 has reinitialised 
        ! save_flow.ind
        CALL WAIT()

        ! Save a copy of the velocity field files corresponding to the 
        ! last time the energy was shifted. I move all the out.h5 files 
        ! to a folder called "last_ET_iteration_flow_field". It is 
        ! created with the execution of initialise_simulation.sh 

        ! I copy the initial condition to tag it with the shifts counter
        
        ! --------------------------------------------------------------
        ! FILE MANAGEMENT
        ! --------------------------------------------------------------

        IF ( RANK_G == 0 ) THEN

        ! Clear last_ET_it_flow_field and copy all the *out.h5 files 
        ! into it so I have a copy of the latest rescaling run no 
        ! no matters it's laminar or turbulent

        call execute_command_line('rm -f last_ET_it_flow_field/*')
        call execute_command_line('cp *out.h5 last_ET_it_flow_field/')

          SELECT CASE ( ENERGY_RESCALING_FLAG )
  
            CASE(-1) ! Relaminarisation

              ! Name of the initial condition of the latest energy
              ! rescaling iteration
              WRITE( LATEST_INITIAL_SOL_STR, '(i3.3, a)') 
     &               SHIFTS_COUNT-1, '_0000_out.h5'

              ! Name to store the latest initial condition
              WRITE( INIT_SOL_FNAME, '(a, i3.3, a)') 
     &              '0000_out_', SHIFTS_COUNT-1, '_lam.h5'

              call execute_command_line('cp "'// 
     &                      trim( LATEST_INITIAL_SOL_STR )// 
     &             '" "' // trim(INIT_SOL_FNAME)// '"')

              ! Remove al the previous files from the destination folder
              call execute_command_line( 
     &            'rm -f latest_laminar_solutions/*')


              ! Move all the solutions to latest_laminar_solutions
              call execute_command_line( 
     &             'mv *out.h5 latest_laminar_solutions/')

            !-----------------------------------------------------------  
            
            CASE (1) ! Turbulence

              ! Name of the initial condition of the latest energy
              ! rescaling iteration
              WRITE( LATEST_INITIAL_SOL_STR, '(i3.3, a)') 
     &               SHIFTS_COUNT-1, '_0000_out.h5'

              ! Name to store the latest initial condition
              WRITE(INIT_SOL_FNAME, '(a, i3.3, a)') 
     &        '0000_out_', SHIFTS_COUNT-1, '_turb.h5'
  
              call execute_command_line('cp "'// 
     &                      trim( LATEST_INITIAL_SOL_STR )// 
     &             '" "' // trim(INIT_SOL_FNAME)// '"')

             ! Remove al the previous files from the destination folder
             call execute_command_line( 
     &           'rm -f latest_turbulent_solutions/*')


             ! Move all the solutions to latest_turbulent_solutions
             call execute_command_line( 
     &            'mv *out.h5 latest_turbulent_solutions/')

          END SELECT
  
        END IF

        ! All the processors wait until proc 0 has moved the files
        CALL WAIT()

        ! Save the new initial condition
        CALL SAVE_FLOW(.FALSE.)
  
        CALL POISSON_P_CHAN
  
        RETURN
   
      END SUBROUTINE ET_ENERGY_REINITIALIZATION


      SUBROUTINE ET_LAMBDA_REINITIALIZATION( LAMBDA_RESCALING_FLAG )
        
        INCLUDE 'header'

        CHARACTER*512 FNAME
        CHARACTER*512 INIT_SOL_FNAME

        ! This flag tells me where the simulation went to last iteration
        !  1 --> it went to turbulence
        ! -1 --> it decayed to laminar

        INTEGER , INTENT(IN) :: LAMBDA_RESCALING_FLAG
        REAL    VERSION, CURRENT_VERSION
        logical RESET_TIME
        INTEGER I, J, K, N

        WRITE(*,*) 'Reinitialising LAMBDA for Edge Tracking'
  
        CALL WAIT()

        ! Initialize storage arrays.
        DO K=0,NZ+1
          DO I=0,NX+1 
            DO J=0,NY+1
              U1(I,K,J)=0.d0
              U3(I,K,J)=0.d0
              U2(I,K,J)=0.d0
              P (I,K,J)=0.d0
              R1(I,K,J)=0.d0
              R2(I,K,J)=0.d0
              R3(I,K,J)=0.d0
              F1(I,K,J)=0.d0
              F2(I,K,J)=0.d0
              F3(I,K,J)=0.d0
              BU1(I,K,J)=0.d0
              DBU1(I,K,J) = 0.D0
              BU1F(I,K,J)=0.d0
              DO N = 1,N_TH
              BTH(I,K,J,N)=0.d0
              END DO
            END DO
          END DO
        END DO

        ! Initialize FFT package (includes defining the wavenumber 
        ! vectors).
        CALL INIT_FFT
        
        ! Store the kinetic energy time series of the latest
        ! ET iteration
        KINS_n(:) = KINS(:)
        
        ! Reinitialisation of the energy time series
        KINS(:)         = 0.D0
        ENERGIES(:)     = 0.D0
        DISSIPATIONS(:) = 0.D0
        
        ENERGYCOUNT     = 1
        DISSCOUNT       = 1
      
        ! Initialize RKW3 parameters.
        H_BAR(1)=DELTA_T*(8.0d0/15.0d0)
        H_BAR(2)=DELTA_T*(2.0d0/15.0d0)
        H_BAR(3)=DELTA_T*(5.0d0/15.0d0)
        BETA_BAR(1)=1.0d0
        BETA_BAR(2)=25.0d0/8.0d0
        BETA_BAR(3)=9.0d0/4.0d0
        ZETA_BAR(1)=0.0d0
        ZETA_BAR(2)=-17.0d0/8.0d0
        ZETA_BAR(3)=-5.0d0/4.0d0
        H_ALPHA(1)=DELTA_T*(0.D0)
        H_ALPHA(2)=DELTA_T*(8.D0/15.D0)
        H_ALPHA(3)=DELTA_T*(2.D0/3.D0)
        DISSIPATION = 0.D0
        NSAMPLES = 0

        ! I create a new flow field using the lambda value updated in 
        ! ET_LAMBDA_EVAL
        CALL READ_FLOW_BISECTION( LAMBDA )        
  
        ! Time counters
        TIME_STEP = 0
        TIME      = T0FILE
  
        ! This variable is used to create a new flow
        ! with the right phase
        TIME_TIME_STEP = TIME

        CALL CREATE_BACKGROUND_FLOW_CHAN
      
        ! Save the new energy value
        CALL GET_ENERGY(.FALSE.)
        CALL SAVE_STATS(.FALSE.)
  
        ! I write to lambda_shift_history.dat the starting point
        ! information (count and physical time) and the lambda shift
        ! information (count, current lambda value and energy values
        ! associated)

        CALL WAIT()

        IF( RANK_G==0 ) CALL WRITE_LAMBDA_SHIFT_HISTORY(.FALSE.,.FALSE.)
          
        ! All the processors wait until proc 0 has reinitialised 
        ! save_flow.ind
        CALL WAIT()

        ! Save the new initial condition
        ! I commented this, because at the beginning of the while loop, 
        ! there is a CALL SAVE_FLOW, which should be the 0000_out.h5 
        ! instead of this one
        
        ! CALL SAVE_FLOW(.FALSE.)
  
        CALL POISSON_P_CHAN
  
        CALL WAIT()

        RETURN
   
      END SUBROUTINE ET_LAMBDA_REINITIALIZATION

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ET_SP_REINITIALIZATION( )
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

        ! --------------------------------------------------------------
        ! DESCRIPTION:
        ! ------------
        ! 
        ! This routine performs the reinitilisation of the flow to start
        ! a new iteration in the RESCALING mode. It also makes a copy
        ! of the *out.h5 files of the last iteration
        ! 
        ! 
        ! The steps it does are:
        !
        ! 1. Set all the flow variables and time series to zero
        !
        ! 2. Reads again the flow variables from start.h5
        !
        ! 3. Redefines INIT_E as INIT_E = ALPHA_E * INIT_E0, where:
        !
        !         INIT_E0: Energy read from start.h5 when INITIALIZE()
        !                  was called at the beginning of the simualtion
        ! 
        !         ALPHAE : Rescaling parameter calculated in ET_ES_EVAL
        !                  as  
        !
        !                  ALPHA_E = 0.5D0*( LATEST_ALPHA_E_LAM + 
        !                                    LATEST_ALPHA_E_TURB  )
        !
        ! 4. Rescale the flow variables with this new INIT_E by calling
        !        
        !       CALL GET_ENERGY(.TRUE.)
        !
        !    where the '.TRUE.' flag indicates that the flow has to be 
        !    rescaled using the value of INIT_E      
        !
        ! 5. WRITING LOG FILES
        !    
        !    it calls WRITE_ENERGY_SHIFT_HISTORY(.FALSE.), where the 
        !    .FALSE. flag is to tell the subroutine not to write the
        !    header (which is written only in INITIALIZE). This routine
        !    updates the file energy_shift_history.dat with the
        !    REINITIALISED values of
        !
        !    SHIFTS_COUNT , ALPHA_E , INIT_E  
        !
        !
        ! 6. FILE MANAGEMENT
        !    
        !    6.1 It removes all the files from the folder
        !        'last_ET_it_flow_field' and copy all the *out.h5 files
        !        into it
        ! 
        !    6.2 It renames the file 000_out.h5 either as 
        !        000_out_lam.h5 or 000_out_turb.h5 depending on the 
        !        value of ENERGY_RESCALING_FLAG (either -1 or 1).
        !
        !    6.3 Depending on the value ENERGY_RESCALING_FLAG, it  
        !        either cleans up the folder latest_laminar_solutions 
        !        or latest_turbulent_solutions and moves all the *out.h5
        !        files into it
        !                          
        ! --------------------------------------------------------------

        INCLUDE 'header'

        ! This routines does:
        ! 1. Find a new starting point for the edge tracking simulation
        ! 2. Reinitialise the flow field
        ! 3. Reinitialise the lambda values and flags 

        logical RESET_TIME
        INTEGER J, K, N
        REAL*8 DE

        INTEGER :: MIN_ECOUNT_COMMON, ENERGY_BISECTION_COUNTER
        INTEGER :: I
        INTEGER :: FLOW_SP_IDX
        CHARACTER*512 COUNT_STR
        CHARACTER*512 FNAME
        
        CHARACTER*512 T0FILE_STR
        CHARACTER*512 LAMBDA_TURB_STR
        CHARACTER*512 LAMBDA_LAM_STR
        CHARACTER*512 SP_SHIFTS_COUNT_STR
        CHARACTER*512 LTEST_LMBD_LAM_COUNT_STR
        CHARACTER*512 LTEST_LMBD_TURB_COUNT_STR

        CALL WAIT()

      ! I save these values in strings to copy the corresponding files
      ! afterwards
      WRITE ( T0FILE_STR            ,       *  ) T0FILE
      WRITE ( LAMBDA_TURB_STR       ,       *  ) LATEST_LAMBDA_TURB
      WRITE ( LAMBDA_LAM_STR        ,       *  ) LATEST_LAMBDA_LAM
      WRITE ( SP_SHIFTS_COUNT_STR   , '(I3.3)' ) SP_SHIFTS_COUNT
      
      WRITE ( LTEST_LMBD_LAM_COUNT_STR  , '(I3.3)' ) 
     &        LATEST_LAMBDA_LAM_COUNT
      
      WRITE ( LTEST_LMBD_TURB_COUNT_STR , '(I3.3)' ) 
     &        LATEST_LAMBDA_TURB_COUNT

      ! Flags and indexes reinitialisation

      LAMBDA_SHIFTS_COUNT      = 0
      LATEST_LAMBDA_LAM_COUNT  = 0
      LATEST_LAMBDA_TURB_COUNT = 0

      SP_SHIFTS_COUNT     = SP_SHIFTS_COUNT + 1
      
      LAMBDA_LAM_FLAG  = .FALSE.
      LAMBDA_TURB_FLAG = .FALSE.
      
      LATEST_LAMBDA_TURB = 1.0D0
      LATEST_LAMBDA_LAM  = 0.0D0
      
      LAMBDA = 0.5D0 * ( LATEST_LAMBDA_TURB + LATEST_LAMBDA_LAM )

      MIN_ECOUNT_COMMON = MIN( ENERGYCOUNT_LAM-1 , ENERGYCOUNT_TURB-1 )

      ENERGYCOUNT_LAM  = 1
      ENERGYCOUNT_TURB = 1

      ! Find the new starting point

      ! I need to find the index. I move backwards to perform less
      ! iterations

      I = MIN_ECOUNT_COMMON

      DO WHILE ( I > 3 .AND. ABS(  KINS_TURB(I) 
     &                           - KINS_LAM (I) ) > ENERGY_BISECT_THRS )
        I = I-1

      END DO

      ! With I and assuming the *_out.h5 files are stored starting from 
      ! 0, I can find which file I need to retrieve to reinitialise the
      ! starting point

      ! Also, I need to assume that  SAVE_FLOW_INT is multiple of
      ! SAVE_STAT_INT

      FLOW_SP_IDX = NINT( ( REAL ( I , 8 )                / 
     &                      REAL ( SAVE_STATS_INT , 8 ) ) /
     &                      REAL ( SAVE_FLOW_INT  , 8 )         ) - 1

      WRITE( COUNT_STR , '(1I0.4)' ) FLOW_SP_IDX

      ! I Copy the new start_turb.h5 and start_lam.h5 files from the 
      ! last couple of simulations and then I call the file reading
      ! routines

      CALL WAIT()

      IF ( RANK_G == 0 ) THEN

        !---------------------------------------------------------------
        !
        ! I make a copy of start_turb.h5 with the values of T0FILE and 
        ! LATEST_LAMBDA_TURB and LATEST_LAMBDA_LAM
        !
        !              start_turb.h5 
        !                   |
        !                   |
        !                   v 
        !   start_turb_SPSHIFTSCOUNTS_T0FILE_LATESTLAMBDATURB.h5
        !
        !
        !              start_lam.h5 
        !                   |
        !                   |
        !                   v 
        !   start_lam_SPSHIFTSCOUNTS_T0FILE_LATESTLAMBDALAM.h5
        !
        !
        !
        ! Example: 
        ! ---------
        !              start_turb.h5 
        !                   |
        !                   |
        !                   v 
        ! start_turb_0008_0.11800000071525574_0.75000000000000000.h5
        !
        !              start_lam.h5 
        !                   |
        !                   |
        !                   v 
        ! start_lam_0008_0.11800000071525574_0.50000000000000000.h5
        !

        call execute_command_line('cp start_turb.h5 start_turb_'  
     &            // TRIM ( ADJUSTL ( SP_SHIFTS_COUNT_STR ) ) // '_'
     &            // TRIM ( ADJUSTL ( T0FILE_STR          ) ) // '_'
     &            // TRIM ( ADJUSTL ( LAMBDA_TURB_STR     ) ) // '.h5')

        call execute_command_line('cp start_lam.h5 start_lam_'  
     &            // TRIM ( ADJUSTL ( SP_SHIFTS_COUNT_STR ) ) // '_'
     &            // TRIM ( ADJUSTL ( T0FILE_STR          ) ) // '_'
     &            // TRIM ( ADJUSTL ( LAMBDA_LAM_STR      ) ) // '.h5')

        !---------------------------------------------------------------
        !
        ! Getting new start.h5 files from the latest saved laminar and 
        ! turbulent solutions
        ! 
        ! latest_turbulent_solutions/COUNT_STR_out.h5
        !                                   |
        !                                   |
        !                                   v
        !                              start_turb.h5
        !
        ! latest_laminar_solutions/COUNT_STR_out.h5
        !                                   |
        !                                   |
        !                                   v
        !                              start_lam.h5
        !

        ! copy the turbulent solution for the new starting point
        FNAME = SAVPATH(:LSP)//'latest_turbulent_solutions/*'//
     &          TRIM( COUNT_STR )//'_out.h5'

        call execute_command_line('cp '//trim(FNAME)//' start_turb.h5') 

        ! copy the laminar solution for the new starting point
        FNAME = SAVPATH(:LSP)//'latest_laminar_solutions/*'//
     &          TRIM( COUNT_STR )//'_out.h5'

        call execute_command_line('cp '//trim(FNAME)//' start_lam.h5') 

        !---------------------------------------------------------------
        ! Saving the h5 files before moving to a new starting point
        IF ( SAVE_ALL_BSCTN_FFS ) THEN

          print *, ' '
          print *, 'Compressing h5 files'
          print *, ' '

          ! taring turbulent h5 files
          FNAME = './latest_turbulent_solutions/turbulent_flow_field_'  
     &            // TRIM ( ADJUSTL ( SP_SHIFTS_COUNT_STR ) )  // '_'
     &            // TRIM ( ADJUSTL ( T0FILE_STR          ) )  // '_'
     &            // TRIM ( ADJUSTL ( LAMBDA_TURB_STR     ) )  // '.tar'

          call execute_command_line('find ./latest_turbulent_solutions/'
     &    //' -name "*.h5" -exec tar -cvzf '// TRIM( FNAME ) // ' {} +')


          ! taring laminar h5 files
          FNAME = './latest_laminar_solutions/laminar_flow_field_'  
     &            // TRIM ( ADJUSTL ( SP_SHIFTS_COUNT_STR ) )  // '_'
     &            // TRIM ( ADJUSTL ( T0FILE_STR          ) )  // '_'
     &            // TRIM ( ADJUSTL ( LAMBDA_LAM_STR      ) )  // '.tar'

          call execute_command_line('find ./latest_laminar_solutions/'
     &    //' -name "*.h5" -exec tar -cvzf '// TRIM( FNAME ) // ' {} +')

       END IF ! SAVE_ALL_BSCTN_FFS

      END IF ! ( RANK_G == 0 )

      ! Wait until proc 0 has finished with the compression
      CALL WAIT()

      WRITE(*,*) 'Reinitialising Flow for Edge Tracking by BISECTION'
  
      ! Initialize storage arrays.
      DO K=0,NZ+1
        DO I=0,NX+1 
          DO J=0,NY+1
          
            U1(I,K,J)=0.d0
            U3(I,K,J)=0.d0
            U2(I,K,J)=0.d0
            P (I,K,J)=0.d0
            R1(I,K,J)=0.d0
            R2(I,K,J)=0.d0
            R3(I,K,J)=0.d0
            F1(I,K,J)=0.d0
            F2(I,K,J)=0.d0
            F3(I,K,J)=0.d0
            BU1(I,K,J)=0.d0
            DBU1(I,K,J) = 0.D0
            BU1F(I,K,J)=0.d0
            
            DO N = 1,N_TH
              BTH(I,K,J,N)=0.d0
            END DO
          
          END DO
        END DO
      END DO
      
      ! Initialize FFT packg (includes defining the wavenumber vectors).
      CALL INIT_FFT

      KINS_TURB = 0.D0
      KINS_LAM  = 0.D0

      ! Reinitialisation of the energy time series
      KINS(:)         = 0.D0
      ENERGIES(:)     = 0.D0
      DISSIPATIONS(:) = 0.D0

      ENERGYCOUNT     = 1
      DISSCOUNT       = 1
    
      ! laminar flow field initial energy
      CALL READ_FLOW_BISECTION( 0.D0 )        

      CALL GET_ENERGY(.FALSE.)
      ENERGY_LAM_INI = ENERGIES(1)

      print *, 'ENERGY_LAM_INI = ', ENERGY_LAM_INI

      KINS(:)         = 0.D0
      KINS_n(:)       = 0.D0 
      ENERGIES(:)     = 0.D0
      DISSIPATIONS(:) = 0.D0

      ENERGYCOUNT     = 1
      DISSCOUNT       = 1
          
      ! turbulent flow field initial energy
      CALL READ_FLOW_BISECTION( 1.D0 )        

      CALL GET_ENERGY(.FALSE.)
      ENERGY_TURB_INI = ENERGIES(1)

      print *, 'ENERGY_TURB_INI = ', ENERGY_TURB_INI

      ENERGY_AVG = 0.5D0 * ( ENERGY_LAM_INI + ENERGY_TURB_INI )
      DE         = ABS( ENERGY_TURB_INI - ENERGY_LAM_INI )
      
      ! I restore the arrays as they were meant to be
      KINS(:)         = 0.D0
      KINS_n(:)       = 0.D0 
      ENERGIES(:)     = 0.D0
      DISSIPATIONS(:) = 0.D0

      ENERGYCOUNT     = 1
      DISSCOUNT       = 1

      ! Initialize RKW3 parameters.
      H_BAR(1)=DELTA_T*(8.0d0/15.0d0)
      H_BAR(2)=DELTA_T*(2.0d0/15.0d0)
      H_BAR(3)=DELTA_T*(5.0d0/15.0d0)
      BETA_BAR(1)=1.0d0
      BETA_BAR(2)=25.0d0/8.0d0
      BETA_BAR(3)=9.0d0/4.0d0
      ZETA_BAR(1)=0.0d0
      ZETA_BAR(2)=-17.0d0/8.0d0
      ZETA_BAR(3)=-5.0d0/4.0d0
      H_ALPHA(1)=DELTA_T*(0.D0)
      H_ALPHA(2)=DELTA_T*(8.D0/15.D0)
      H_ALPHA(3)=DELTA_T*(2.D0/3.D0)
      DISSIPATION = 0.D0
      NSAMPLES = 0

      ! reads start_lam.h5 and start_turb.h5 and computes the new flow 
      ! field as uλ = ulam + λ * ( uturb - ulam )
      CALL READ_FLOW_BISECTION( LAMBDA )
      ! Save the new energy value
      CALL GET_KE_ONLY( ENERGY_BISECTION )

      IF ( RANK_G == 0 ) THEN
        PRINT *, 'ENERGY( LAMBDA ) INI = ', ENERGY_BISECTION
      END IF

      ! I iterate over LAMBDA until the energy is close the midpoint
      ! between laminar and turbulent
      ENERGY_BISECTION_COUNTER = 0

      CALL WAIT()

      DO WHILE( ABS( ENERGY_BISECTION - ENERGY_AVG ) > 0.05D0 * DE
     &      .AND. ENERGY_BISECTION_COUNTER < 30 )

        ! Slightly modifying LAMBDA to converge to the midpoint
        IF ( ENERGY_BISECTION > ENERGY_AVG ) LAMBDA = LAMBDA - 0.05
        IF ( ENERGY_BISECTION < ENERGY_AVG ) LAMBDA = LAMBDA + 0.05

        ! read and update the velocity field
        CALL READ_FLOW_BISECTION( LAMBDA )
        ! compute the KE without writing the files
        CALL GET_KE_ONLY( ENERGY_BISECTION )

        ENERGY_BISECTION_COUNTER = ENERGY_BISECTION_COUNTER + 1

      END DO

      CALL WAIT()
          
      IF ( RANK_G == 0 ) THEN
        PRINT *, 'ENERGY( LAMBDA ) CONVERGED = ', ENERGY_BISECTION
        PRINT *, 'AFTER = ', ENERGY_BISECTION_COUNTER, ' ITERS '
      END IF

      CALL GET_ENERGY(.FALSE.)

      ! Time counters
      TIME_STEP = 0
      TIME      = T0FILE
  
      ! I write the shift log files with the latest update of the values
      IF ( RANK_G == 0 ) THEN
        
        ! Write lambda history file with no header but with separator
        CALL WRITE_LAMBDA_SHIFT_HISTORY(.FALSE.,.TRUE.)
        CALL WRITE_SP_SHIFT_HISTORY(.FALSE.)        

      END IF
        
      CALL POISSON_P_CHAN
  
      ! All the processors wait until proc 0 has reinitialised 
      ! save_flow.ind
      CALL WAIT()

      RETURN  

      END SUBROUTINE ET_SP_REINITIALIZATION


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ET_CHECK_STOP_COND( ET_COND )
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

        INCLUDE 'header'
        LOGICAL :: ET_COND
        REAL*8  :: MAX_DIFF
        INTEGER :: N_TSTEPS_DIFFERENCE
        INTEGER :: MAX_LOOP

        N_TSTEPS_DIFFERENCE = INT(   NINT( T_ET_COMPARISON/DELTA_T )
     &                             / SAVE_STATS_INT ) + 1
        
        MAX_DIFF = 0.D0

        DO MAX_LOOP = 1, N_TSTEPS_DIFFERENCE

          ! I update the maximum value if a new one comes through
          MAX_DIFF = MAX( MAX_DIFF , ABS( KINS  (MAX_LOOP) - 
     &                                    KINS_n(MAX_LOOP) ) )

        END DO

        IF ( MAX_DIFF < TOL_ET_ITERATION ) THEN
          
          IF ( RANK_G == 0 ) THEN

            PRINT *, ' '
            PRINT *, ' **************************************'
            PRINT *, ' '
            PRINT *, '    Edge-Tracking stop condition met!'
            PRINT *, ' '
            PRINT *, ' **************************************'
            PRINT *, ' '
          END IF

          ET_COND = .FALSE.
    
        END IF 
  
      END SUBROUTINE ET_CHECK_STOP_COND


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      LOGICAL FINAL

        CALL SAVE_STATS_CHAN(FINAL)          
   
  
         if (flavor.eq.'Basic') then
         IF (RANK_G.EQ.0) THEN
                write(*,*) 'done save_stats diablo'
              END IF
         end if

        IF (FINAL) THEN
      
            CALL VIS_FLOW_CHAN         
  
        END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE READ_FLOW
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*512 FNAME
      CHARACTER*512 FNAME_TH(N_TH)
      INTEGER I, J, K, N, NUM_PER_DIR_T

      !=================================================================
      ! DEFINING FNAME
      !=================================================================

#ifdef HDF5
      IF (READ_HDF5) THEN
        FNAME = INPATH(:LIP)//'start.h5'
      ELSE IF (USE_MPI) THEN
        FNAME = INPATH(:LIP)//'diablo_'//MPI_IO_NUM//'.start'
      ELSE
        FNAME=INPATH(:LIP)//'diablo.start'
      END IF
#else
      IF (USE_MPI) THEN
        FNAME = INPATH(:LIP)//'diablo_'//MPI_IO_NUM//'.start'
      ELSE
        FNAME=INPATH(:LIP)//'diablo.start'
      END IF
#endif

      IF (RANK_G.EQ.0) THEN
        WRITE(6,*) 'Reading flow from ',FNAME
      END IF

      !=================================================================
      ! READING INITIAL SOLUTION
      !=================================================================

      IF (FNAME(LEN_TRIM(FNAME)-2:LEN_TRIM(FNAME)).EQ.".h5") THEN
#ifdef HDF5      

         CALL READHDF5(FNAME)

         ! For edge tracking time thresholds (minimum time to evaluate
         ! an energy shift and the minimum size for the average window)
         T0FILE = TIME
#else
         write(*,*) ' **** ERROR ******************************'
         write(*,*) ' Program not compiled with HDF5 libraries.'
         stop       
#endif      
      
      ELSE

         DO N=1,N_TH
            IF (USE_MPI) THEN
                FNAME_TH(N)=INPATH(:LIP)//'diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) //'_'//MPI_IO_NUM// '.start'
            ELSE
               FNAME_TH(N)=INPATH(:LIP)//'diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) // '.start'
            END IF
         END DO
         WRITE(6,*)   'Reading flow from ',FNAME

         OPEN(UNIT=10,FILE=FNAME,STATUS="OLD",FORM="UNFORMATTED")
         READ (10) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP


         write(*,*) 'NX_T, NY_T, NZ_T: ',NX_T,NY_T,NZ_T

         ! TEMPORARY, UNCOMMENT!!!
         !      IF ((NX .NE. NX_T) .OR. (NY .NE. NY_T) .OR. (NZ .NE. NZ_T))
         !     *     STOP 'Error: old flowfield wrong dimensions. '
         !      IF (NUM_PER_DIR .NE. NUM_PER_DIR_T)
         !     *     STOP 'Error: old flowfield wrong NUM_PER_DIR. '

         write(*,*) 'READING FLOW'

         READ (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *             (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *             (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *             (((CP(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1)

         DO N=1,NUM_READ_TH
! Specify in input.dat which scalars are to be read
            OPEN(UNIT=11,FILE=FNAME_TH(READ_TH_INDEX(N)),STATUS="OLD"
     &              ,FORM="UNFORMATTED")
            READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP
            READ (11) (((CTH(I,K,J,READ_TH_INDEX(N))
     &           ,I=0,NKX),K=0,TNKZ),J=0,NY+1)
            CLOSE(11)
         END DO
 
         CLOSE(10)
         CLOSE(11)
      
      END IF
      
      !=================================================================

C Apply initial boundary conditions, set ghost cells
      IF (USE_MPI) THEN
        call APPLY_BC_VEL_MPI
      ELSE
        call APPLY_BC_VEL_LOWER
        call APPLY_BC_VEL_UPPER
      END IF

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE READ_FLOW_BISECTION( LAMBD )
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      CHARACTER*512 FNAME_LAM
      CHARACTER*512 FNAME_TURB
      REAL*8 LAMBD

      FNAME_TURB = INPATH(:LIP)//'start_turb.h5'
      FNAME_LAM  = INPATH(:LIP)//'start_lam.h5'


      !IF (FNAME(LEN_TRIM(FNAME)-2:LEN_TRIM(FNAME)).EQ.".h5") THEN
#ifdef HDF5      
         CALL READHDF5_BISECTION( FNAME_TURB , FNAME_LAM , LAMBD )

         ! For edge tracking time thresholds (minimum time to evaluate
         ! an energy shift and the minimum size for the average window)
         T0FILE = TIME

#else
         write(*,*) ' **** ERROR ******************************'
         write(*,*) ' Program not compiled with HDF5 libraries.'
         stop       
#endif      
      
      !END IF

C Apply initial boundary conditions, set ghost cells
      IF (USE_MPI) THEN
        call APPLY_BC_VEL_MPI
      ELSE
        call APPLY_BC_VEL_LOWER
        call APPLY_BC_VEL_UPPER
      END IF

      RETURN
      
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_FLOW(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*512 FNAME
      CHARACTER*512 FNAME_TH(N_TH)
      INTEGER      I, J, K, N
      LOGICAL      FINAL, FLAG
      INTEGER ID
      CHARACTER*512 STR

      CHARACTER*512 ENERGY_SHIFTS_COUNT_STR
      CHARACTER*512 SP_SHIFTS_COUNT_STR
      CHARACTER*512 LAMBDA_SHIFTS_COUNT_STR

      CHARACTER*512 counter_str

      ! Shift indexes as strings
      WRITE ( ENERGY_SHIFTS_COUNT_STR , '(I3.3)' ) SHIFTS_COUNT
      WRITE ( SP_SHIFTS_COUNT_STR     , '(I3.3)' ) SP_SHIFTS_COUNT
      WRITE ( LAMBDA_SHIFTS_COUNT_STR , '(I3.3)' ) LAMBDA_SHIFTS_COUNT

      IF (FINAL) THEN
#ifdef HDF5
            IF(SAVE_HDF5) THEN
              FNAME = SAVPATH(:LSP)//'end.h5'
            ELSE IF (USE_MPI) THEN
              FNAME = SAVPATH(:LSP)//'diablo_'//MPI_IO_NUM//'.res'
            ELSE
              FNAME=SAVPATH(:LSP)//'diablo.res'
            END IF
#else      
            IF (USE_MPI) THEN
               FNAME = SAVPATH(:LSP)//
     &              'diablo_'//MPI_IO_NUM//'.res'
            ELSE
               FNAME=SAVPATH(:LSP)//'diablo.res'
            END IF
#endif            
         
      ELSE

#ifdef HDF5

        IF (SAVE_HDF5) THEN
          FNAME = SAVPATH(:LSP)//'out.h5'                
        ELSE IF (USE_MPI) THEN
          FNAME = SAVPATH(:LSP)//
     &          'diablo_'//MPI_IO_NUM//'.saved'
        ELSE
          FNAME=SAVPATH(:LSP)//'diablo.saved'
        END IF
#else
        IF (USE_MPI) THEN
          FNAME = SAVPATH(:LSP)//
     &          'diablo_'//MPI_IO_NUM//'.saved'
        ELSE
          FNAME=SAVPATH(:LSP)//'diablo.saved'
        END IF            
#endif        
         
      END IF
      
      ! check that I'm working with h5 outfiles
      IF ( FNAME( LEN_TRIM(FNAME)-2 : LEN_TRIM(FNAME) ).EQ.".h5") THEN
      
#ifdef HDF5
         
        ! check if save_flow.ind exists
        inquire( file=SAVPATH(:LSP)//'save_flow.ind',exist=flag )
         
        if (flag) then
         
          open(unit=500,file=SAVPATH(:LSP)//
     &         'save_flow.ind',status='old'
     &         ,form='formatted')
         
          read(500,'(1I10)') id
         
          close(unit=500) 
         
          write( counter_str ,'(1I0.4)' ) id

          if ( ET_MODE ) then

            if ( ET_RESCALING ) then
            
              str = TRIM( ENERGY_SHIFTS_COUNT_STR )//'_'//
     &              TRIM( counter_str             )
            
            end if

            if ( ET_BISECTION ) then

              str =  TRIM( SP_SHIFTS_COUNT_STR     ) //'_'//
     &               TRIM( LAMBDA_SHIFTS_COUNT_STR ) //'_'//
     &               TRIM( counter_str             )
            
            end if

          else ! not ET_MODE

            write( str,'(1I0.4)' ) id

          end if         

          FNAME = SAVPATH(:LSP) // TRIM(str) // '_' // 'out.h5'
         
        end if
         
        call WriteHDF5(FNAME)
         
        if (flag) then
           
           ! update save_flow.ind 
           open(unit=500 , file = SAVPATH(:LSP)//'save_flow.ind',
     &          form='formatted')
           write(500,'(1I10)') id+1
           close(unit=500)
        
        end if
#else
        write(*,*) ' **** ERROR ******************************'
        write(*,*) ' Program not compiled with HDF5 libraries.'
        stop
#endif      
      
      ELSE ! not h5 output

        IF (FINAL) THEN
          DO N=1,N_TH
            IF (USE_MPI) THEN
              FNAME_TH(N) = SAVPATH(:LSP)//'diablo_th'
     &         //CHAR(MOD(N,100)/10+48)
     &         //CHAR(MOD(N,10)+48) //'_'//MPI_IO_NUM// '.res'             
            ELSE
              FNAME_TH(N)=SAVPATH(:LSP)//'diablo_th'
     &         //CHAR(MOD(N,100)/10+48)
     &         //CHAR(MOD(N,10)+48) // '.res'
            END IF
 
          END DO      
        ELSE
       
          DO N=1,N_TH
            IF (USE_MPI) THEN
              FNAME_TH(N) = SAVPATH(:LSP)//'diablo_th'
     &         //CHAR(MOD(N,100)/10+48)
     &         //CHAR(MOD(N,10)+48) //'_'//MPI_IO_NUM// '.saved'
            ELSE
              FNAME_TH(N)=SAVPATH(:LSP)//'diablo_th'
     &         //CHAR(MOD(N,100)/10+48)
     &         //CHAR(MOD(N,10)+48) // '.saved'
            END IF
 
 
          END DO      
       
        END IF ! IF (FINAL)
      
      
        IF ( RANK_G == 0 ) WRITE(6,*) 'Writing flow to ',FNAME
  
        OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")
        WRITE(10) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP



        WRITE (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CP(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1)
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN"
     &       ,FORM="UNFORMATTED")
          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=0,NY+1)
          CLOSE(11)
        END DO

        CLOSE(10)
        CLOSE(11)
      
      END IF ! not h5 output

      RETURN
      
      END SUBROUTINE SAVE_FLOW
      
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE GET_ENERGY(SET_ENERGY)
C----*|--.---------.---------.---------.---------.---------.---------.-|--      

      INCLUDE 'header'
      REAL*8 TEMP1(0:NY+1)
      REAL*8 TEMP2(0:NY+1)
      REAL*8 ENERGY,ANS,KINENERGY,POTENERGY
      LOGICAL, INTENT(IN) :: SET_ENERGY            
      INTEGER I,J,K,N
      COMPLEX*16 CTHTEMP(0:NKX,0:TNKZ,1:NY,1:MAX(1,N_TH))
      CHARACTER*512 FNAMET,FNAMEK,FNAMEP

      FNAMET = SAVPATH(:LSP)//'./tot.txt'     
      FNAMEK = SAVPATH(:LSP)//'./kin.txt'     
      FNAMEP = SAVPATH(:LSP)//'./pot.txt'     

      OPEN(UNIT=20,FILE=FNAMET,ACCESS='APPEND',STATUS='OLD')
      OPEN(UNIT=40,FILE=FNAMEK,ACCESS='APPEND',STATUS='OLD')
      OPEN(UNIT=60,FILE=FNAMEP,ACCESS='APPEND',STATUS='OLD')

      TEMP1(:) = 0.d0
      TEMP2(:) = 0.D0
      ENERGY = 0.d0
      ANS = 0.d0

      DO J = 1,NY
        DO K = 0,TNKZ
          TEMP1(J) = TEMP1(J) + LX*LZ
     &                  * CDABS(CU1(0,K,J))**2.0d0
          TEMP1(J) = TEMP1(J) + LX*LZ 
     &                  * CDABS(CU3(0,K,J))**2.0d0
          DO N = 1,N_TH
              TEMP2(J) = TEMP2(J) + RI_TAU(N) * LX*LZ 
     &         * ABS(CTHTEMP(0,K,J,N))**2.0d0
          END DO
          DO I = 1,NKX
            TEMP1(J) = TEMP1(J) + 2.0d0*LX*LZ
     &                  * CDABS(CU1(I,K,J))**2.0d0
            TEMP1(J) = TEMP1(J) + 2.0d0*LX*LZ
     &                  * CDABS(CU3(I,K,J))**2.0d0
            DO N = 1,N_TH
              TEMP2(J) = TEMP2(J) + RI_TAU(N) * 2.0d0*LX*LZ
     &         * CDABS(CTH(I,K,J,N))**2.0d0
            END DO
          END DO
        END DO
      END DO 
      CALL TRAPZ(TEMP1,KINENERGY,1)
      CALL TRAPZ(TEMP2,POTENERGY,1)

      TEMP1(:) = 0.d0
      
      DO J = 1,NY
        DO K = 0,TNKZ
          TEMP1(J) = TEMP1(J) + LX*LZ
     &                  * CDABS(CU2(0,K,J))**2.0d0
          DO I = 1,NKX
            TEMP1(J) = TEMP1(J) + 2.0d0 * LX*LZ 
     &                  * CDABS(CU2(I,K,J))**2.0d0
          END DO
        END DO
      END DO
      CALL TRAPZ(TEMP1,ANS,2)
      
      KINENERGY = KINENERGY + ANS
      
      KINENERGY = KINENERGY/2.0d0   
      POTENERGY = POTENERGY/2.0D0
      
      KINENERGY = KINENERGY/(LX*LZ*LY)                                      
      POTENERGY = POTENERGY/(LX*LZ*LY)
             
      CALL GATHER_SCALAR_MPI(KINENERGY,KINENERGY)
      CALL GATHER_SCALAR_MPI(POTENERGY,POTENERGY)
      ENERGY = KINENERGY + POTENERGY


      ! If I'm in Edge Tracking mode, it doesn't matter the INIT_E given
      ! in input.dat. What the model is going to consider is the one that
      ! start.h5 has got.

      !IF ( ET_MODE ) THEN
      !  INIT_E = ENERGY
      !END IF
      
      IF (RANK_G.EQ.0) THEN
        !WRITE(*,*) 'ENERGY:', ENERGY
        WRITE(20,*) ENERGY
        WRITE(40,*) KINENERGY
        WRITE(60,*) POTENERGY                       
      END IF

      CLOSE(20)
      CLOSE(40)
      CLOSE(60)
      
      IF (SET_ENERGY) THEN      
        
        ! Rescaling of the flow field 
        CU1(:,:,:)   = CU1(:,:,:)   * DSQRT( INIT_E / ENERGY )
        CU2(:,:,:)   = CU2(:,:,:)   * DSQRT( INIT_E / ENERGY )
        CU3(:,:,:)   = CU3(:,:,:)   * DSQRT( INIT_E / ENERGY )
        CTH(:,:,:,:) = CTH(:,:,:,:) * DSQRT( INIT_E / ENERGY )      
                         
        ! Update ghost cells with the new flow field                         
        IF (USE_MPI) THEN
          CALL GHOST_CHAN_MPI
        END IF                           
                         
        ! Get the pressure from the poisson equation                         
        CALL POISSON_P_CHAN
      
        ! Update ghost cells again                         
        IF (USE_MPI) THEN
          CALL GHOST_CHAN_MPI
        END IF                         
        
      ELSE               
        
        POTS     ( ENERGYCOUNT ) = POTENERGY
        KINS     ( ENERGYCOUNT ) = KINENERGY
        ENERGIES ( ENERGYCOUNT ) = ENERGY

        ENERGYCOUNT = ENERGYCOUNT + 1
                         
      END IF
                               
      RETURN
      END 


C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE GET_KE_ONLY( KE )
C----*|--.---------.---------.---------.---------.---------.---------.-|--      

      INCLUDE 'header'
      REAL*8 TEMP1(0:NY+1)
      REAL*8 TEMP2(0:NY+1)
      REAL*8 ENERGY,ANS,KINENERGY,POTENERGY, KE
      INTEGER I,J,K,N
      COMPLEX*16 CTHTEMP(0:NKX,0:TNKZ,1:NY,1:MAX(1,N_TH))
      CHARACTER*512 FNAMET,FNAMEK,FNAMEP

      TEMP1(:) = 0.d0
      TEMP2(:) = 0.D0
      ENERGY = 0.d0
      ANS = 0.d0

      DO J = 1,NY
        DO K = 0,TNKZ
          TEMP1(J) = TEMP1(J) + LX*LZ
     &                  * CDABS(CU1(0,K,J))**2.0d0
          TEMP1(J) = TEMP1(J) + LX*LZ 
     &                  * CDABS(CU3(0,K,J))**2.0d0
          DO N = 1,N_TH
              TEMP2(J) = TEMP2(J) + RI_TAU(N) * LX*LZ 
     &         * ABS(CTHTEMP(0,K,J,N))**2.0d0
          END DO
          DO I = 1,NKX
            TEMP1(J) = TEMP1(J) + 2.0d0*LX*LZ
     &                  * CDABS(CU1(I,K,J))**2.0d0
            TEMP1(J) = TEMP1(J) + 2.0d0*LX*LZ
     &                  * CDABS(CU3(I,K,J))**2.0d0
            DO N = 1,N_TH
              TEMP2(J) = TEMP2(J) + RI_TAU(N) * 2.0d0*LX*LZ
     &         * CDABS(CTH(I,K,J,N))**2.0d0
            END DO
          END DO
        END DO
      END DO 
      CALL TRAPZ(TEMP1,KINENERGY,1)
      CALL TRAPZ(TEMP2,POTENERGY,1)

      TEMP1(:) = 0.d0
      
      DO J = 1,NY
        DO K = 0,TNKZ
          TEMP1(J) = TEMP1(J) + LX*LZ
     &                  * CDABS(CU2(0,K,J))**2.0d0
          DO I = 1,NKX
            TEMP1(J) = TEMP1(J) + 2.0d0 * LX*LZ 
     &                  * CDABS(CU2(I,K,J))**2.0d0
          END DO
        END DO
      END DO
      CALL TRAPZ(TEMP1,ANS,2)
      
      KINENERGY = KINENERGY + ANS
      
      KINENERGY = KINENERGY/2.0d0   
      
      KINENERGY = KINENERGY/(LX*LZ*LY)                                      
      POTENERGY = POTENERGY/(LX*LZ*LY)
             
      CALL GATHER_SCALAR_MPI(KINENERGY,KINENERGY)
      
      KE = KINENERGY
                               
      RETURN
      END       
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE GET_DISSIPATION(LAST_TIME)
C----*|--.---------.---------.---------.---------.---------.---------.-|--

! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps

      include 'header'
      character*512 FNAME
      LOGICAL LAST_TIME      
      REAL*8 NEW_DISS
      REAL*8 ANS
      REAL*8 VECTOR(0:NY+1)

      integer I,J,K,N

      IF (LAST_TIME) THEN                                          
                       
      IF (RANK_G.EQ.0) THEN     
        WRITE(*,*) 'Dissipation = ',DISSIPATION        
      END IF
            
      ELSE

      FNAME = SAVPATH(:LSP)//'diss.txt'

      OPEN(UNIT=80,FILE=FNAME,ACCESS='APPEND',STATUS='OLD')

      NEW_DISS = 0.d0
      ANS = 0.d0
      VECTOR(:) = 0.d0
      NSAMPLES = NSAMPLES + 1

! DU1DX
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J)+NU*LX*LZ*KX2(0)*CDABS(CU1(0,K,J))**2.0d0
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KX2(I)*CDABS(CU1(I,K,J))**2.0d0
          END DO
        END DO
      END DO
      CALL TRAPZ(VECTOR,ANS,1)
      NEW_DISS = NEW_DISS + ANS
      
! DU1DY
      VECTOR(:) = 0.d0      
      
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ
     &            *(CDABS(CU1(0,K,J)-CU1(0,K,J-1))**2.0d0)/(DY(J)*DY(J))
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 2.0d0*NU*LX*LZ
     &            *(CDABS(CU1(I,K,J)-CU1(I,K,J-1))**2.0d0)/(DY(J)*DY(J))                    
          END DO
        END DO            
      END DO       
      CALL TRAPZ(VECTOR,ANS,2)    
      NEW_DISS = NEW_DISS + ANS      
      
! DU1DZ
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ*KZ2(K)*CDABS(CU1(0,K,J))
     &          **2.0d0
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KZ2(K)*CDABS(CU1(I,K,J))**2.0d0
          END DO         
        END DO
      END DO 
      CALL TRAPZ(VECTOR,ANS,1)
      NEW_DISS = NEW_DISS + ANS
      
! DU3DX
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ*KX2(0)*CDABS(CU3(0,K,J))
     &          **2.0d0
          DO I = 1,NKX    
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KX2(I)*CDABS(CU3(I,K,J))**2.0d0
          END DO          
        END DO
      END DO
      CALL TRAPZ(VECTOR,ANS,1)
      NEW_DISS = NEW_DISS + ANS
      
! DU3DY
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ
     &            *(CDABS(CU3(0,K,J)-CU3(0,K,J-1))**2.0d0)/(DY(J)*DY(J))
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 2.0d0*NU*LX*LZ
     &            *(CDABS(CU3(I,K,J)-CU3(I,K,J-1))**2.0d0)/(DY(J)*DY(J))
          END DO         
        END DO
      END DO  
      CALL TRAPZ(VECTOR,ANS,2)      
      NEW_DISS = NEW_DISS + ANS    
      
! DU3DZ
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ*KZ2(K)*CDABS(CU3(0,K,J))
     &          **2.0d0
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KZ2(K)*CDABS(CU3(I,K,J))**2.0d0
          END DO          
        END DO
      END DO  
      CALL TRAPZ(VECTOR,ANS,1)
      NEW_DISS = NEW_DISS + ANS
      
! DU2DX
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + 
     &             NU*LX*LZ*KX2(0)*CDABS(CU2(0,K,J))**2.0d0
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KX2(I)*CDABS(CU2(I,K,J))**2.0d0
          END DO                    
        END DO
      END DO
      CALL TRAPZ(VECTOR,ANS,2)
      NEW_DISS = NEW_DISS + ANS
      
! DU2DY
      DO J = 1,NY
      VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ
     &            *(CDABS(CU2(0,K,J+1)-CU2(0,K,J))**2.0d0)/
     &     (DYF(J)*DYF(J))
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 2.0d0*NU*LX*LZ
     &            *(CDABS(CU2(I,K,J+1)-CU2(I,K,J))**2.0d0)/
     &     (DYF(J)*DYF(J))
          END DO         
        END DO
      END DO           
      CALL TRAPZ(VECTOR,ANS,1) 
      NEW_DISS = NEW_DISS + ANS
      
! DU2DZ
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + 
     &             NU*LX*LZ*KZ2(K)*CDABS(CU2(0,K,J))**2.0d0
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KZ2(K)*CDABS(CU2(I,K,J))**2.0d0
          END DO          
        END DO
      END DO  
      CALL TRAPZ(VECTOR,ANS,2)   
      NEW_DISS = NEW_DISS + ANS                            
      
! DTHDX
!      DO N = 1,N_TH
!      DO J = 1,NY
!        VECTOR(J) = 0.d0
!        DO K = 0,TNKZ
!          VECTOR(J) = VECTOR(J) + 
!     &       NU*LX*LZ*KX2(0)*CDABS(CTH(0,K,J,N))**2.0d0
!     &            *RI_TAU(N)/PR(N)
!          DO I = 1,NKX
!            VECTOR(J) = VECTOR(J) + 
!     &       2.0d0*NU*LX*LZ*KX2(I)*CDABS(CTH(I,K,J,N))**2.0d0
!     &            *RI_TAU(N)/PR(N)
!          END DO         
!        END DO
!      END DO   
!      CALL TRAPZ(VECTOR,ANS,1) 
!      NEW_DISS = NEW_DISS + ANS
!      END DO
      
! DTHDY
!      DO N = 1,N_TH
!      DO J = 1,NY
!        VECTOR(J) = 0.d0
!        DO K = 0,TNKZ
!          VECTOR(J) = VECTOR(J) + 
!     &             NU*LX*LZ
!     &          *(CDABS(CTH(0,K,J,N)-CTH(0,K,J-1,N))**2.0d0)/
!     &           (DY(J)*DY(J))
!     &            *RI_TAU(N)/PR(N)
!          DO I = 1,NKX
!            VECTOR(J) = VECTOR(J) + 
!     &             2.0d0*NU*LX*LZ
!     &          *(CDABS(CTH(I,K,J,N)-CTH(I,K,J-1,N))**2.0d0)/
!     &           (DY(J)*DY(J))
!     &            *RI_TAU(N)/PR(N)
!          END DO          
!        END DO
!      END DO  
!      CALL TRAPZ(VECTOR,ANS,2)   
!      NEW_DISS = NEW_DISS + ANS 
!      END DO
      
! DTHDZ
!      DO N = 1,N_TH
!      DO J = 1,NY
!        VECTOR(J) = 0.d0
!        DO K = 0,TNKZ
!          VECTOR(J) = VECTOR(J) + 
!     &             NU*LX*LZ*KZ2(K)*CDABS(CTH(0,K,J,N))**2.0d0
!     &            *RI_TAU(N)/PR(N)
!          DO I = 1,NKX
!            VECTOR(J) = VECTOR(J) + 
!     &             2.0d0*NU*LX*LZ*KZ2(K)*CDABS(CTH(I,K,J,N))**2.0d0
!     &            *RI_TAU(N)/PR(N)
!          END DO          
!        END DO
!      END DO    
!      CALL TRAPZ(VECTOR,ANS,1)
!      NEW_DISS = NEW_DISS + ANS
!      END DO      
      
      NEW_DISS = NEW_DISS / (LX*LY*LZ)
      
      IF (USE_MPI) THEN
        CALL GATHER_SCALAR_MPI(NEW_DISS,NEW_DISS)
      END IF      
      
      IF (RANK_G.EQ.0) THEN
        WRITE(80,*) NEW_DISS
      END IF
      DISSIPATIONS(DISSCOUNT) = NEW_DISS
      DISSCOUNT = DISSCOUNT + 1
      
      CLOSE(80)

! TRAPEZOIDAL RULE IN TIME      
      IF (TIME_AVERAGE) THEN
      IF (NSAMPLES.EQ.1) THEN
      DISSIPATION = DISSIPATION + 0.5d0 * NEW_DISS / FLOAT(N_TIME_STEPS)
      ELSEIF (NSAMPLES.EQ.N_TIME_STEPS) THEN
      DISSIPATION = DISSIPATION + 0.5d0 * NEW_DISS / FLOAT(N_TIME_STEPS)
      ELSE
        DISSIPATION = DISSIPATION + NEW_DISS / FLOAT(N_TIME_STEPS)
      END IF
      ELSE
        DISSIPATION = NEW_DISS
      END IF                  
      
! TIME AVERAGED
!      DISSIPATION = (1.d0/FLOAT(NSAMPLES)*NEW_DISS)
!     &      +((FLOAT(NSAMPLES)-1.d0)/FLOAT(NSAMPLES))*DISSIPATION
      
      END IF
      
      RETURN
      END   
      
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|      
      SUBROUTINE END_RUN(FLAG)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
      INCLUDE 'header'

      LOGICAL FLAG,FILE_EXISTS
      
      FLAG=.FALSE.
      ! Check for the time
      call WALL_TIME(END_TIME)
      if (END_TIME-START_TIME.gt.TIME_LIMIT) THEN
         write(*,*) ' STOP beacuse of wall-time hit!'
         FLAG=.TRUE.
      END IF
      
      INQUIRE(FILE="stop.now", EXIST=FILE_EXISTS)
      IF ( FILE_EXISTS ) THEN
         write(*,*) ' STOP beacuse of stop.now file!'
         FLAG=.TRUE.
      END IF
      
      RETURN
      END 


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|      
      SUBROUTINE CHECK_ET_INPUTS_ET(ISTAT)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
        INCLUDE 'header'

        INTEGER ISTAT
      
        ISTAT = 0

        ! CASES WHERE THE PROGRAM WOULDN'T WORK
        IF ( ET_MODE .AND. .NOT. ET_RESCALING .AND. 
     &                     .NOT. ET_BISECTION ) THEN

          ISTAT = 1
          print *, 'CHECK ET FLAGS IN input.dat'

        END IF

        ! CASES WHERE THE PROGRAM WOULDN'T WORK
        IF ( ET_RESCALING .AND. ET_BISECTION ) THEN

          ISTAT = 2
          print *, 'ET_RESCALING and ET_BISECTION cannot be both TRUE'

        END IF

        IF ( WALL_STOKES .AND. PRESSURE_STOKES ) THEN

          ISTAT = 3
          print *, 'WALL_STOKES and PRESSURE_STOKES cannot be both TRUE'

        END IF

      END 

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|      
      SUBROUTINE INIT_RESUME_BISECTION()
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
        INCLUDE 'header'

        real*8 aux
        integer i

        ! starting point shifts
        open( unit=500 , file=SAVPATH(:LSP)//'sp_shifts.ind', 
     &        status='old',form='formatted')

          read(500,*) SP_SHIFTS_COUNT

        close(unit=500) 

        ! lambda shifts 
        open( unit=501 , file=SAVPATH(:LSP)//'lambda_shifts.ind', 
     &        status='old',form='formatted')

          read(501,*) LAMBDA_SHIFTS_COUNT ! integer
          read(501,*) LATEST_LAMBDA_LAM   ! real
          read(501,*) LATEST_LAMBDA_TURB  ! real
          read(501,*) LAMBDA_LAM_FLAG     ! logical
          read(501,*) LAMBDA_TURB_FLAG    ! logical

          ! DO i = 1 , size( KINS_n )
          !   read(501,*) KINS_n(i)
          ! END DO

        close(unit=501) 

        ! If the values were entered in the wrong order, I correct
        ! that here
        if ( LATEST_LAMBDA_LAM > LATEST_LAMBDA_TURB ) then

          print *, 'Lambda values were entered in the wrong order'

          aux                = LATEST_LAMBDA_LAM
          LATEST_LAMBDA_LAM  = LATEST_LAMBDA_TURB
          LATEST_LAMBDA_TURB = aux

        end if

        ! this counters put me in the first position of the arrays that
        ! append the kinetic energy
        ENERGYCOUNT_LAM   = 1
        ENERGYCOUNT_TURB  = 1


      END 

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|      
      SUBROUTINE WRITE_ENERGY_SHIFT_HISTORY( WRITE_HEADER )
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
        INCLUDE 'header'

        LOGICAL WRITE_HEADER

        IF( WRITE_HEADER ) THEN

          print *, ' '
          print *, 'UPDATING HEADER energy_shift_history.dat ...'
          print *, ' '

          ! Open file
          OPEN( UNIT=19,FILE=SAVPATH(:LSP)//'energy_shift_history.dat', 
     &          status='replace')
  
          ! Write header
          WRITE(19, '(A,3x,A,16x,A)') 'SHIFTS_COUNT','ALPHA_E','INIT_E'
          
          CLOSE(19)

          print *, ' '
          print *, 'UPDATING HEADER energy_shift_history.dat DONE'
          print *, ' '

        END IF

        ! I write to energy_shift_history.dat, the alpha value and 
        ! the new energy after the last energy shift

        print *, ' '
        print *, 'UPDATING energy_shift_history.dat ...'
        print *, ' '

        OPEN(unit=19,file=SAVPATH(:LSP)// 'energy_shift_history.dat',  
     &        status = 'old', position = 'append', action = 'write')
        
        WRITE(19, '(I10,5X,E21.15E2,2X,E21.15E2)') 
     &        SHIFTS_COUNT, ALPHA_E, INIT_E
        
        CLOSE(19)

        print *, ' '
        print *, 'UPDATING energy_shift_history.dat DONE'
        print *, ' '

        ! I reset save.flow.ind
        OPEN( unit = 20, file = SAVPATH(:LSP)//'save_flow.ind',  
     &        status = 'replace', action = 'write')
        WRITE(20, '(I0)') 0
        CLOSE(20)

      END 


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|      
      SUBROUTINE WRITE_LAMBDA_SHIFT_HISTORY( WRITE_HEADER, WRITE_SEP )
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
        INCLUDE 'header'

        LOGICAL WRITE_HEADER , WRITE_SEP
        CHARACTER*512 FMAT1, FMAT2, FMAT
        CHARACTER*154 SEPARATOR
        INTEGER I

        ! Initialize the line with dashes
        DO I = 1, 154
         SEPARATOR(I:I) = '-'
        END DO

        IF( WRITE_HEADER ) THEN

          ! Initialize the line with dashes
          DO I = 1, 154
           SEPARATOR(I:I) = '='
          END DO

          print *, ' '
          print *, 'UPDATING HEADER lambda_shift_history.dat ...'
          print *, ' '

          ! Open file
          OPEN( UNIT=19,FILE=SAVPATH(:LSP)//'lambda_shift_history.dat', 
     &          status='replace')
  
          ! Write header
          WRITE(19, '(A,3x,A,3x,A,16x,A,17x,A,8x,A,8x,A)') 
     &              'SP_SHIFTS_COUNT'  , 'LAMBDA_SHIFTS_COUNT',  
     &              'T0FILE '          ,  'LAMBDA'            , 
     &              'ENERGY_TURB_INI'  , 'ENERGY_LAM_INI '    , 
     &              'ENERGY(LAMBDA)'
          
          WRITE(19, '(A)') SEPARATOR

          CLOSE(19)

          print *, ' '
          print *, 'UPDATING HEADER lambda_shift_history.dat DONE'
          print *, ' '

        END IF

        ! The format string is too long in this case, so I needed to
        ! put it into two strings
        FMAT1 = TRIM('(I10,12X,I10,8X,E21.15E2,2X,E21.15E2,2X,')
        FMAT2 = TRIM( 'E21.15E2,2X,E21.15E2,2X,E21.15E2)')
 
        WRITE(FMAT,*) TRIM(FMAT1) // TRIM(FMAT2)

        print *, ' '
        print *, 'UPDATING lambda_shift_history.dat ...'
        print *, ' '

        OPEN(unit=19,file=SAVPATH(:LSP)//'lambda_shift_history.dat',  
     &        status = 'old', position = 'append', action = 'write')
        
        IF ( WRITE_SEP ) WRITE(19, '(A)') SEPARATOR

        WRITE(19, FMAT) SP_SHIFTS_COUNT, LAMBDA_SHIFTS_COUNT, T0FILE , 
     &                  LAMBDA , ENERGY_TURB_INI, ENERGY_LAM_INI      ,  
     &                  ENERGIES(1)    

        CLOSE(19)

        print *, ' '
        print *, 'UPDATING lambda_shift_history.dat DONE'
        print *, ' '

        ! I reset save.flow.ind
        OPEN( unit = 20, file = 'save_flow.ind',  
     &        status = 'replace', action = 'write')
        WRITE(20, '(I0)') 0
        CLOSE(20)

      END 

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|      
      SUBROUTINE WRITE_SP_SHIFT_HISTORY( WRITE_HEADER )
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
        INCLUDE 'header'

        LOGICAL WRITE_HEADER
        CHARACTER*154 SEPARATOR
        INTEGER I

        ! Initialize the line with dashes
        DO I = 1, 154
         SEPARATOR(I:I) = '-'
        END DO

        IF( WRITE_HEADER ) THEN
          ! Initialize the line with dashes
          DO I = 1, 154
           SEPARATOR(I:I) = '='
          END DO

          print *, ' '
          print *, 'UPDATING HEADER sp_shifts_history.dat ...'
          print *, ' '

          ! Open file
          OPEN( UNIT=19,FILE=SAVPATH(:LSP)//'sp_shifts_history.dat', 
     &          status='replace')
  
          ! Write header
          WRITE(19, '(A,3x,A,16x,A)') 
     &    'SP_SHIFTS_COUNT','T0FILE','ENERGY'
          
          WRITE(19, '(A)') SEPARATOR

          CLOSE(19)

          print *, ' '
          print *, 'UPDATING HEADER sp_shifts_history.dat DONE'
          print *, ' '

        END IF

        print *, ' '
        print *, 'UPDATING sp_shifts_history.dat ...'
        print *, ' '

        ! I write to energy_shift_history.dat, the alpha value and 
        ! the new energy after the last energy shift
        OPEN(unit=19,file=SAVPATH(:LSP)//'sp_shifts_history.dat',  
     &        status = 'old', position = 'append', action = 'write')
        
        WRITE(19, '(I10,5X,E21.15E2,2X,E21.15E2)') 
     &        SP_SHIFTS_COUNT, T0FILE, ENERGIES(1)
        
        CLOSE(19)

        print *, ' '
        print *, 'UPDATING sp_shifts_history.dat DONE'
        print *, ' '

        ! I reset save.flow.ind
        OPEN( unit = 20, file = 'save_flow.ind',  
     &        status = 'replace', action = 'write')
        WRITE(20, '(I0)') 0
        CLOSE(20)

      END 


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|      
      SUBROUTINE WRITE_BISECTION_LOG()
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
        INCLUDE 'header'

        integer i

        print *, ' '
        print *, 'UPDATING sp_shifts.ind ...'
        print *, ' '

        ! starting point shifts
        open( unit=500 , file=SAVPATH(:LSP)//'sp_shifts.ind', 
     &        form='formatted')

          write(500,*) SP_SHIFTS_COUNT

        close(unit=500) 

        print *, ' '
        print *, 'UPDATING sp_shifts.ind DONE'
        print *, ' '

        ! lambda shifts 

        print *, ' '
        print *, 'UPDATING lambda_shifts.ind ...'
        print *, ' '

        open( unit=501 , file=SAVPATH(:LSP)//'lambda_shifts.ind', 
     &        status='replace',form='formatted')

          write(501,*) LAMBDA_SHIFTS_COUNT ! integer
          write(501,*) LATEST_LAMBDA_LAM   ! real
          write(501,*) LATEST_LAMBDA_TURB  ! real
          write(501,*) LAMBDA_LAM_FLAG     ! logical
          write(501,*) LAMBDA_TURB_FLAG    ! logical

          !DO i = 1 , size(KINS)
          !  write(501,*) KINS(i)
          !END DO

        close(unit=501) 

        print *, ' '
        print *, 'UPDATING lambda_shifts.ind DONE'
        print *, ' '

      END 


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine wall_time(wt)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c
c     Return wall-clock time as seconds after Jan. 1, 2014.
c     Support for leap year is not included anymore.
c
c     Next leap year is 2016!
c
c     By using a 'save' statement, the wall-time after the first
c     call to the subroutine could be computed, but that is not
c     intended with the present subroutine (e.g. the history file)
c
      implicit none

      real*8 wt
      integer val(8),i,shift,day

      integer mon(12,2)
      data mon /
     &     31,28,31,30,31,30,31,31,30,31,30,31,
     &     31,29,31,30,31,30,31,31,30,31,30,31/
c
c     Get current date and time
c     val(1) : year
c     val(2) : month
c     val(3) : day
c     val(4) : difference to GMT
c     val(5) : hour
c     val(6) : minute
c     val(7) : second
c     val(8) : 1/1000 second
c
      call date_and_time(values=val)
c
c     Determine leap year
c
      if (mod(val(1),4).eq.0) then
         if (mod(val(1),100).eq.0) then
            if (mod(val(1),400).eq.0) then
               shift=2
            else
               shift=1
            end if
         else
            shift=2
         end if
      else
         shift = 1
      end if
c
c     Construct day of the year
c
      day = val(3)-1
      do i=1,val(2)-1
         day=day+mon(i,shift)
      end do
c
c     And compute wall-clock time
c
      wt = (val(1)-2014)*365*86400+
     &     day*86400+val(5)*3600+val(6)*60+val(7)+val(8)/1000.

      end       
