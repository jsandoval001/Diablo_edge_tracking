C******************************************************************************|
C diablo.f -> DNS In A Box, Laptop Optimized                       VERSION 0.9
C
C This Fortran 77 code computes incompressible flow in a box.
C
C Primative variables (u,v,w,p) are used, and continuity is enforced with a
C fractional step algorithm.
C
C SPATIAL DERIVATIVES:
C   0, 1, 2, or 3 directions are taken to be periodic and handled spectrally
C   (these cases are referred to as the "periodic", "channel", "duct", and
C    "cavity" cases respectively).
C   The remaining directions are taken to be bounded by walls and handled with
C   momentum- and energy-conserving second-order central finite differences.
C
C TIME ADVANCEMENT
C   Two main approaches are implemented:
C     1. RKW3 on nonlinear terms and CN on viscous terms over each RK substep.
C     2. RKW3 on y-derivative terms and CN on other terms over each RK substep.
C
C The emphasis in this introductory code is on code simplicity:
C   -> All variables are in core.
C   -> The code is not explicitly designed for use with either MPI or SMP.
C   -> Overindexing is not used.
C A few simple high-performance programming constructs are used:
C   -> The inner 2 loops are broken out in such a way as to enable out-of-order
C      execution of these loops as much as possible, thereby leveraging
C      vector and superscalar CPU architectures.
C   -> The outer loops are fairly long (including as many operations as
C      possible inside on a single J plane of data) in order to make effective
C      use of cache.
C Multiple time advancement algorithms are implemented for the periodic,
C channel, duct, and cavity cases in order to compare their efficiency for
C various flows on different computational architectures.  In a few of the
C algorithms, the code is broken out fully into a PASS1/PASS2 architecture
C to maximize the efficient use of cache.
C
C This code was developed as a joint project in MAE 223 (CFD), taught by
C Thomas Bewley, at UC San Diego (spring 2001, spring 2005).
C Primary contributions follow:
C Thomas Bewley was the chief software architect
C John R. Taylor wrote the channel flow solvers
C******************************************************************************|
C
C This code is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the
C Free Software Foundation; either version 2 of the License, or (at your
C option) any later version. This code is distributed in the hope that it
C will be useful, but WITHOUT ANY WARRANTY; without even the implied
C warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C GNU General Public License for more details. You should have received a
C copy of the GNU General Public License along with this code; if not,
C write to the Free Software Foundation, Inc., 59 Temple Place - Suite
C 330, Boston, MA 02111-1307, USA.
C
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      PROGRAM DIABLO
       !   implicit none
      INCLUDE 'header'
      INTEGER N
      LOGICAL FLAG

      LOGICAL MIN_TIME_AVG_COND , ET_LOOP , ET_COND  
      INTEGER ENERGY_RESCALING_FLAG
      INTEGER LAMBDA_RESCALING_FLAG
      LOGICAL SP_SHIFT_FLAG

      WRITE(6,*) 
      WRITE(6,*) '             ****** WELCOME TO DIABLO ******'
      WRITE(6,*) '                      VERTICAL CODE'
      WRITE(6,*)

!      OPEN(11,FILE='/media/data-2/stokes_edge_tracking/wall/Re_1000/'//
!     &     'first_guests_2023_11_03/005_1e-08/path.dat'
!     &     ,FORM='FORMATTED',
!     &     STATUS='OLD')
      OPEN(11,FILE='path.dat',FORM='FORMATTED', STATUS='OLD')
      READ(11,'(A)') INPATH
      READ(11,'(A)') SAVPATH

      CLOSE(11)
      WRITE(*,*) 'INPATH: ',INPATH
      WRITE(*,*) 'OUTPATH: ',SAVPATH
   
      LIP = LEN_TRIM(INPATH)
      LSP = LEN_TRIM(SAVPATH)

      CALL INITIALIZE

      ! COMPUTE TIME TAKEN
      CALL WALL_TIME(START_TIME)

      ! A flag to determine if we are considering the first time-step
      FIRST_TIME=.TRUE.

      ! CALL SAVE_FLOW(.FALSE.)

      ! EDGE TRACKING LOOP. IF EDGE_TRACKING_MODE == .FALSE., THEN IT 
      ! WORKS AS THE ORIGINAL VERSION

      MIN_TIME_AVG_COND = .FALSE.

      ET_LOOP = .TRUE. ! SET TO .FALSE. AFTER THE FIRST ITERATION
                       ! IF EDGE_TRACKING_MODE == .FALSE, SO IT
                       ! WORKS AS THE ORIGINAL VERSION

      ET_COND = .TRUE. ! WHENEVER THE STOP CONDITION FOR THE 
                       ! EDGE TRACKING ALGORITHM IS MET, THE
                       ! LOOP ENDS.

      IF ( APPLY_SYMMETRIES ) CALL APPLY_SYMMETRY()

      ! EDGE-TRACKING (ET) LOOP
      ET: DO WHILE ( ET_LOOP .AND. ET_COND ) 

        ! If the simulation is not in Edge Tracking mode, it just run over
        ! the main loop (fixed initial energy)
        IF ( .NOT. ET_MODE ) THEN
          ET_LOOP = .FALSE.
        END IF

        CALL SAVE_FLOW(.FALSE.)

        ENERGY_RESCALING_FLAG = 0
        LAMBDA_RESCALING_FLAG = 0

        ! MAIN LOOP ( TIME-MARCHING (TM) )
        TM: DO TIME_STEP = TIME_STEP+1, TIME_STEP+N_TIME_STEPS
      
          CALL WAIT()

          ! Runge-Kutta step to advance the solution in time
          DO RK_STEP=1,3                            
            CALL RK_CHAN_1         
          END DO
          
          IF ( APPLY_SYMMETRIES ) CALL APPLY_SYMMETRY()

          ! Getting the energy and dissipation and writing them onto
          ! the corresponding log files (kin.txt and dissipation.txt)
          IF ( MOD( TIME_STEP , SAVE_STATS_INT ) == 0 ) THEN
            CALL GET_DISSIPATION (.FALSE.)      
            CALL GET_ENERGY      (.FALSE.)
          END IF
                  
          TIME = TIME + DELTA_T
          FIRST_TIME=.FALSE.
   ! Save statistics to an output file
   !        IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
   !            CALL SAVE_STATS(.FALSE.)
   !        END IF
           
   ! Save the flow to a restart file
          IF ( MOD( TIME_STEP , SAVE_FLOW_INT ) == 0 ) THEN

            IF (RANK_G == 0) THEN
              
              PRINT *, ' '
              WRITE(6,*) 'BEginning TIME_STEP = ', TIME_STEP
              WRITE(6,*) 'TIME                = ', TIME + DELTA_T       
              WRITE(6,*) 'TIME - T0FILE       = ', TIME + DELTA_T-T0FILE       
              WRITE(6,*) 'ENERGY              = ', KINS(ENERGYCOUNT-1)       
          
            END IF

            CALL SAVE_FLOW(.FALSE.)
   !        CALL GET_ENERGY(.FALSE.)
          END IF
           
   ! Save the full three-dimensional fields in NETCDF format to vis.nc
   ! (optional)
          IF (MOD(TIME_STEP,SAVE_FLOW_INT).EQ.0) THEN
            CALL NETCDF_OPEN_VIS
            CALL NETCDF_WRITE_VIS
            CALL NETCDF_CLOSE_VIS
          END IF
   
   ! Filter the scalar field
   !      DO N = 1,N_TH
   !        IF (FILTER_TH(N).AND.(MOD(TIME_STEP,FILTER_INT(N)).EQ.0)) THEN
   !          IF (RANK_G.EQ.0) THEN    
   !          write(*,*) 'Filtering...'
   !          END IF
   !          CALL FILTER_CHAN
   !        END IF 
   !      END DO
   
          ! Check wether to stop the simulation or continue it
          IF (USE_MPI) THEN
            CALL END_RUN_MPI(FLAG)
          ELSE
            CALL END_RUN(FLAG)
          END IF
         
          IF (FLAG) THEN
            EXIT ET
          END IF
   
         ! -------------------------------------------------------------
         ! ENERGY SHIFT EVALUATION 
         ! 
         ! If the simulation has run for longer than 
         ! MIN_T_ENERGY_SHIFT_EVALUATION, and for a period of time 
         ! longer than AVERAGE_WINDOW_ENERGY_SHIFT_EVALUATION after 
         ! MIN_T_ENERGY_SHIFT_EVALUATION, then we evaluate if an energy 
         ! shift is required by calling 
         ! EDGE_TRACKING_ENERGY_SHIFT_EVALUATION.
         !
         ! If the conditions are met, ENERGY_RESCALING_FLAG is set to 1 
         ! and the flow is reinitialised with a different initial energy 
         ! value.
         ! 
         !  Energy
         !   .
         !  /|\           TIME MOVING WINDOW FOR AVERAGING THE ENERGY 
         !   |                     |                        |    
         !   |                ---> |<- AVG_WINDOW_ES_EVAL ->| --->
         !   |                     |                        |
         !   |            .    /\      /\      /\      /\   .  
         !   |            .   /  \    /  \    /  \    /  \  . /\
         !   |            .  /    \  /    \  /    \  /    \ ./  \
         !   |        /\/\/\/      \/      \/      \/      \/    \
         !   |       |    .                                 . 
         !   |_______|____|_________________________________|__\ t    
         !           |    |                                 |  /     
         !       T0FILE   |                                 |
         !                |                                 |
         !   T0FILE + MIN_T_ES_EVAL                       TIME-->
         !
         ! T_WINDOW_AVG ≥ AVERAGE_WINDOW_ENERGY_SHIFT_EVALUATION
         ! -------------------------------------------------------------
   
          IF ( ET_MODE ) THEN

            MIN_TIME_AVG_COND = TIME - (   T0FILE + MIN_T_ES_EVAL 
     &                                   + AVG_WINDOW_ES_EVAL ) > 0.D0

            IF ( MIN_TIME_AVG_COND .AND. 
     &           MOD( TIME_STEP , SAVE_FLOW_INT ) == 0 ) THEN      
              
              ! I check the energy shift or lambda shift condition
              ! depending on the Edge Tracking Mode defined in input.dat

              IF ( ET_RESCALING ) THEN
      
                IF(RANK_G == 0)  PRINT *, 'CHECKING ENERGY SHIFT CONDS'
                CALL ET_ES_EVAL ( ENERGY_RESCALING_FLAG )

              END IF

              IF ( ET_BISECTION ) THEN

                IF(RANK_G == 0) PRINT * , 'CHECKING LAMBDA SHIFT CONDS'

                CALL ET_LAMBDA_EVAL ( LAMBDA_RESCALING_FLAG )

                IF(RANK_G == 0) PRINT * , 'LAMBDA_RESCALING_FLAG = ' ,
     &                                     LAMBDA_RESCALING_FLAG         

              END IF

            END IF

          END IF

          IF ( ET_RESCALING .AND. ENERGY_RESCALING_FLAG /= 0 ) THEN
            
            ! CHECK IF I HAVE TO STOP THE EDGE TRACKING LOOP
            ! WHILE ET_COND = .TRUE.,  IT KEEPS ITERATING
            ! OTHERWISE IT EXTIS THE LOOP TO TERMINATES THE\
            ! EXECUTION (GOTO 20)

            !IF ( SHIFTS_COUNT > 1 ) THEN

              !IF (USE_MPI) THEN
              !  CALL ET_CHECK_STOP_COND_MPI( ET_COND )
              !ELSE
              !  CALL ET_CHECK_STOP_COND( ET_COND )
              !END IF
            
            !END IF

            IF (ET_COND) THEN
              
              IF(RANK_G == 0) THEN
                PRINT *, 'STARTING ENERGY_SHIFT ...'
                PRINT *, 'ENERGY_RESCALING_FLAG = ',
     &                    ENERGY_RESCALING_FLAG
              END IF

              CALL ET_ENERGY_REINITIALIZATION( ENERGY_RESCALING_FLAG )
              
              ! REINITIALISE FLAG VARIABLES
              MIN_TIME_AVG_COND = .FALSE.
              ENERGY_RESCALING_FLAG = 0

              IF (RANK_G == 0) THEN
                PRINT *, 'NEW INIT_E = ', INIT_E
              END IF

              EXIT TM ! TO SIMULATE ANOTHER INITIAL ENERGY VALUE
            
            ELSE ! STOP THE EDGE TRACKING
    
              EXIT ET
    
            END IF
   
          END IF ! ENERGY_RESCALING_FLAG /= 0


          IF ( ET_BISECTION .AND. LAMBDA_RESCALING_FLAG /= 0 ) THEN

            IF (RANK_G == 0) PRINT *,
     &      'CALLING ET_SP_EVAL ...'

            CALL ET_SP_EVAL( SP_SHIFT_FLAG )
            
            IF (RANK_G == 0) PRINT *,
     &      'SP_SHIFT_FLAG = ', SP_SHIFT_FLAG

            ! SP_SHIFT_FLAG, shifts starting point forward
            IF ( SP_SHIFT_FLAG ) THEN

              IF (RANK_G == 0) THEN 
                
                PRINT * , 'CALLING ET_SP_REINITIALIZATION ...'
                PRINT * , '         |                        '
                PRINT * , '         |                        '
                PRINT * , '         V                        '
              
              END IF

              CALL ET_SP_REINITIALIZATION()
              
              IF (RANK_G == 0) PRINT * ,
     &        'ET_SP_REINITIALIZATION DONE!'

            ELSE

              IF (RANK_G == 0) THEN 
              
                PRINT *, 'CALLING ET_LAMBDA_REINITIALIZATION ...'
                PRINT * , '         |                           '
                PRINT * , '         |                           '
                PRINT * , '         V                           '
              
              END IF

              CALL ET_LAMBDA_REINITIALIZATION( LAMBDA_RESCALING_FLAG )

              IF (RANK_G == 0) PRINT * ,
     &        'ET_LAMBDA_REINITIALIZATION DONE!'

            END IF

            ! Reinitialise flag variables
            SP_SHIFT_FLAG         = .FALSE.
            MIN_TIME_AVG_COND     = .FALSE.
            LAMBDA_RESCALING_FLAG = 0

            ! Sychronise all the processors before exiting the TM loop
            CALL WAIT ()

            IF (RANK_G == 0) THEN

              PRINT *, ' '
              PRINT *, '============================================== '
              PRINT *, 'STARTING A NEW TM LOOP  ---------------->      '
              PRINT *, '============================================== '
              PRINT *, ' '

            END IF

            EXIT TM

          END IF ! LAMBDA_RESCALING_FLAG /= 0

        END DO TM ! END TIME MARCHING DO

      END DO ET ! DO WHILE ( ET_LOOP .AND. ET_COND )
 
      IF( RANK_G == 0 ) PRINT *, 'END OF EDGE TRACKING LOOP' 

      ! COMPUTE TIME TAKEN
      CALL GET_DISSIPATION(.TRUE.)
      CALL GET_ENERGY(.FALSE.)
      CALL WALL_TIME(END_TIME)
     
      IF (RANK_G.EQ.0) THEN
      WRITE(*,*) 'Elapsed Time (sec): ',END_TIME-START_TIME
      WRITE(*,*) 'Seconds per Iteration: '
     &     ,(END_TIME-START_TIME)/N_TIME_STEPS
      END IF
      
      IF (RANK_G.EQ.0) THEN
!        WRITE(*,*) 'DISSCOUNT ',DISSCOUNT
!        DO N = 1,DISSCOUNT
!          WRITE(*,*) 'DISSIPATIONS ',DISSIPATIONS(N)
!        END DO
!        WRITE(*,*) 'ENERGYCOUNT ',ENERGYCOUNT
!        DO N = 1,ENERGYCOUNT
!          WRITE(*,*) 'POTS ',POTS(N)
!          WRITE(*,*) 'KINS ',KINS(N)
!        END DO
      END IF
      
      IF (USE_MPI) THEN
        CALL WAIT
      END IF
      
      TIME_STEP=TIME_STEP-1
      CALL SAVE_FLOW(.TRUE.)
      CALL SAVE_STATS(.TRUE.)
      WRITE(6,*)
      WRITE(6,*) '        ****** Hello world!  Have a nice day! ******'
      WRITE(6,*)
      END


