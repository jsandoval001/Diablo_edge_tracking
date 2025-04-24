C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE GHOST_CHAN_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'
      INCLUDE 'mpif.h'
      include 'header_mpi'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer arrys is 0:NKX,0:TNKZ,# of variables
      COMPLEX*16 OCPACK(0:NKX,0:TNKZ,4+N_TH)
      COMPLEX*16 ICPACK(0:NKX,0:TNKZ,4+N_TH)

! If we are using more than one processor, then we need to pass data
      IF (NPROCS.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANK.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

! End if RANK=0
      ELSE IF (RANK.lt.NPROCS-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

! At this point we have passed data up the chain

! Now, we need to pass data down the chain of processes
      IF (RANK.eq.NPROCS-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data down
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
      ELSE IF (RANK.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO

        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE GHOST_GRID_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'

      integer i,j,k,N

      real*8 OCPACK(3),ICPACK(3)

! First, Pass data up the chain to higher ranked processes

      IF (RANK.eq.0) THEN

! Set the lower ghost cells
        GYF(0)=2.d0*GYF(1)-GYF(2)
        GY(0)=2.d0*GY(1)-GY(2)

        OCPACK(1)=GYF(NY-1)
        OCPACK(2)=GYF(NY)
        OCPACK(3)=GY(NY)

        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

! End if RANK=0
      ELSE IF (RANK.lt.NPROCS-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        OCPACK(1)=GYF(NY-1)
        OCPACK(2)=GYF(NY)
        OCPACK(3)=GY(NY)

        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(0)=ICPACK(1)
        GYF(1)=ICPACK(2)
        GY(1)=ICPACK(3)

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Set the top ghost cell
        GYF(NY+1)=2.d0*GYF(NY)-GYF(NYM)

! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(0)=ICPACK(1)
        GYF(1)=ICPACK(2)
        GY(1)=ICPACK(3)

      END IF
! At this point we have passed data up the chain

      IF (RANK.eq.NPROCS-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information

        OCPACK(1)=GYF(2)
        OCPACK(2)=GY(2)


! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
      ELSE IF (RANK.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us

        OCPACK(1)=GYF(2)
        OCPACK(2)=GY(2)
        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(NY+1)=ICPACK(1)
        GY(NY+1)=ICPACK(2)

      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        GYF(NY+1)=ICPACK(1)
        GY(NY+1)=ICPACK(2)

      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This version solves for one full row, then passes the data
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...
      implicit none
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX-1,0:NY+1), B(0:NX-1,0:NY+1), C(0:NX-1,0:NY+1)
     &         , G(0:NX-1,0:NY+1)

      REAL*8 OCPACK(4),ICPACK(4)

      DO I = 0,NX-1

      IF (RANK.NE.0) THEN
! If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4,MPI_DOUBLE_PRECISION,RANK-1,12
     &               ,MPI_COMM_WORLD,status,ierror)
! Unpack the data
        A(I,1)=OCPACK(1)
        B(I,1)=OCPACK(2)
        C(I,1)=OCPACK(3)
        G(I,1)=OCPACK(4)
      ELSE
! We are the lowest process, so solve at J=1
        J=1
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
       END IF

! For all processes, solve from 2..NY 
        DO J=2,NY
            A(I,J)=-A(I,J)/B(I,J-1)
            B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
            G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
        END DO

      IF (RANK.NE.NPROCS-1) THEN
        ICPACK(1)=A(I,NY)
        ICPACK(2)=B(I,NY)
        ICPACK(3)=C(I,NY)
        ICPACK(4)=G(I,NY)
        CALL MPI_SEND(ICPACK,4,MPI_DOUBLE_PRECISION,RANK+1,12
     &               ,MPI_COMM_WORLD,ierror)
      ELSE
! UPPER PROCESS, SOLVE ONE MORE ROW
        J = NY+1
        A(I,J) = -A(I,J)/B(I,J-1)
        B(I,J) = B(I,J) + A(I,J)*C(I,J-1)
        G(I,J) = G(I,J) + A(I,J)*G(I,J-1)      
      END IF

      END DO

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_FORWARD_COMPLEX_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...
      implicit none
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX,0:NY+1), B(0:NX,0:NY+1), C(0:NX,0:NY+1)
      COMPLEX*16 G(0:NX,0:NY+1)

      COMPLEX*16 OCPACK(4),ICPACK(4)

      DO I=0,NX

      IF (RANK.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4,MPI_DOUBLE_COMPLEX,RANK-1,13
     &               ,MPI_COMM_WORLD,status,ierror)
C Unpack the data
        A(I,1)=REAL(OCPACK(1))
        B(I,1)=REAL(OCPACK(2))
        C(I,1)=REAL(OCPACK(3))
        G(I,1)=OCPACK(4)
      END IF

      DO J=2,NY
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
      END DO

      IF (RANK.NE.NPROCS-1) THEN
        ICPACK(1)=A(I,NY)
        ICPACK(2)=B(I,NY)
        ICPACK(3)=C(I,NY)
        ICPACK(4)=G(I,NY)
        CALL MPI_SEND(ICPACK,4,MPI_DOUBLE_COMPLEX,RANK+1,13
     &               ,MPI_COMM_WORLD,ierror)
      END IF


      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_REAL_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|------
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...
      implicit none
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX-1,0:NY+1), B(0:NX-1,0:NY+1), C(0:NX-1,0:NY+1)
     &     , G(0:NX-1,0:NY+1)
      REAL*8 ICPACK(1),OCPACK(1)

      DO I = 0,NX-1

       IF (RANK.NE.NPROCS-1) THEN
! If we aren't the highest process, then wait for data
        CALL MPI_RECV(OCPACK,1,MPI_DOUBLE_PRECISION,RANK+1,10
     &               ,MPI_COMM_WORLD,status,ierror)
          G(I,NY+1)=OCPACK(1)
      ELSE
! Else, if we are the highest process, compute the solution at j=NY 
         G(I,NY+1)=G(I,NY+1)/B(I,NY+1)
      END IF

! All processes solve from NY..1
      DO J=NY,0,-1
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      END DO

      IF (RANK.NE.0) THEN
          ICPACK(1)=G(I,2)
        CALL MPI_SEND(ICPACK,1,MPI_DOUBLE_PRECISION,RANK-1,10
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      END DO

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_COMPLEX_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...
      implicit none
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX,0:NY+1), B(0:NX,0:NY+1), C(0:NX,0:NY+1)
      COMPLEX*16 G(0:NX,0:NY+1)

      DO I=0,NX

      IF (RANK.NE.NPROCS-1) THEN
C If we aren't the highest process, then wait for data
        CALL MPI_RECV(G(I,NY+1),1,MPI_DOUBLE_COMPLEX,RANK+1,11
     &               ,MPI_COMM_WORLD,status,ierror)
        J=NY
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      ELSE
C Else, if we are the highest process, then compute the solution at j=NY
        G(I,NY)=G(I,NY)/B(I,NY)
      END IF

      DO J=NY-1,1,-1
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      END DO

      IF (RANK.NE.0) THEN
        CALL MPI_SEND(G(I,2),1,MPI_DOUBLE_COMPLEX,RANK-1,11
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      END DO

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      CALL MPI_INIT(IERROR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERROR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,IERROR)

      RANK_G = RANK

      IF (RANK_G.EQ.0) THEN

      write(*,*) 'MPI Initialized, NPROCS,RANK: ',NPROCS,RANK

      END IF

C Set a string to determine which input/output files to use
C When MPI is used, each process will read/write to files with the
C number of their rank (+1) appended to the end of the file.
C The string MPI_IO_NUM will be used to define the RANK+1 for each process
        IF (NPROCS.le.10) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCS.le.100) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,100)/10+48)
     &             //CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCS.le.1000) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,1000)/100+48)
     &             //CHAR(MOD(RANK+1,100)/10+48)
     &             //CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCS.le.10000) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,10000)/1000+48)
     &             //CHAR(MOD(RANK+1,1000)/100+48)
     &             //CHAR(MOD(RANK+1,100)/10+48)
     &             //CHAR(MOD(RANK+1,10)+48)
        ELSE
           WRITE(6,*) 'ERROR, NPROCS>10,000, Unsupported problem size'
        END IF



      RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_CHAN_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C Initialize any constants here
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER J, N

      IF (RANK_G.EQ.0) THEN

      write(*,*) '*******IN INIT_CHAN_MPI *********'

      END IF

      PI=4.D0*ATAN(1.D0)

! Set the upper and lower bounds for timestepping
      IF (RANK.eq.0) THEN
        write(*,*) 'U_BC_YMIN: ',U_BC_YMIN
        JEND=NY
!        JSTART=2
        JSTART=1
! Now, set the indexing for the scalar equations
        DO N=1,N_TH
          JEND_TH(N)=NY       
!          JSTART_TH(N)=2
          JSTART_TH(N)=1
        END DO
      ELSE IF (RANK.eq.NPROCS-1) THEN
        write(*,*) 'U_BC_YMAX: ',U_BC_YMAX
        JSTART=2
!        JEND=NY-1
        JEND = NY

! Set the upper and lower limits of timestepping of the scalar equations
        DO N=1,N_TH
        JSTART_TH(N)=2
!        JEND_TH(N)=NY-1
        JEND_TH(N) = NY
        END DO

      ELSE
! Here, we are on a middle process
        JSTART=2
        JEND=NY
        DO N = 1,N_TH
          JSTART_TH(N) = 2
          JEND_TH(N) = NY
        END DO
      END IF

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_TH_MPI(MATL,MATD,MATU,VEC,N)
!----*|--.---------.---------.---------.---------.---------.---------.-|------
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi' 
      INTEGER N
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_TH_LOWER(MATL,MATD,MATU,VEC,N)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the upper plane, apply the boundary conditions
          CALL APPLY_BC_TH_UPPER(MATL,MATD,MATU,VEC,N)
        END IF
      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_U2_MPI(MATL,MATD,MATU,VEC)
!----*|--.---------.---------.---------.---------.---------.---------.-|------
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_2_LOWER(MATL,MATD,MATU,VEC)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_2_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_U1_MPI(MATL,MATD,MATU,VEC)
!----*|--.---------.---------.---------.---------.---------.---------.-|------
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_1_LOWER(MATL,MATD,MATU,VEC)        
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_1_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_U3_MPI(MATL,MATD,MATU,VEC)
!----*|--.---------.---------.---------.---------.---------.---------.-|------
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_3_LOWER(MATL,MATD,MATU,VEC)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_3_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_REM_DIV_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
!----*|--.---------.---------.---------.---------.---------.---------.-|------
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      INTEGER I,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
C Apply the boundary conditions
        IF (RANK.EQ.0) THEN
        DO I=0,NKX
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
          IF ((K.EQ.0).AND.(I.EQ.0)) THEN
C Otherwise the matrix will be singular
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
            MATL_C(I,1)=0.d0
            MATD_C(I,1)=1.d0
            MATU_C(I,1)=0.d0
            VEC_C(I,1)=(0.d0,0.d0)
          ELSE
C Use Dirichlet boundary conditions, dp/dz=0 at walls
            MATL_C(I,1)=0.d0
            MATD_C(I,1)=1.d0
            MATU_C(I,1)=-1.d0
            VEC_C(I,1)=(0.d0,0.d0)
          END IF
        END DO
        END IF
C Apply the boundary conditions
        IF (RANK.EQ.NPROCS-1) THEN
          DO I=0,NKX
            MATL_C(I,NY)=1.d0
            MATD_C(I,NY)=-1.d0
            MATU_C(I,NY)=0.d0
            VEC_C(I,NY)=(0.d0,0.d0)
          END DO
        END IF

        RETURN
        END

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_POISSON_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
!----*|--.---------.---------.---------.---------.---------.---------.-|------
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      INTEGER I,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANK.eq.0) THEN
        DO I=0,NKX
          IF ((I.EQ.0).AND.(K.EQ.0)) THEN
            MATD_C(I,1)=1.d0
            MATU_C(I,1)=0.d0
            VEC_C(I,1)=(0.d0,0.d0)
          ELSE
! Here, apply Neumann boundary conditions (dp/dz=0) at the walls
            MATD_C(I,1)=1.d0
            MATU_C(I,1)=-1.d0
            VEC_C(I,1) = (0.d0,0.d0)
              VEC_C(I,1)=VEC_C(I,1)
     &              -DY(2)*((CU2(I,K,3)-CU2(I,K,2))/DYF(2)
     &            -(CU2(I,K,2)-CU2(I,K,1))/DYF(1))
     &             /DY(2)
          END IF
        END DO
        END IF
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANK.eq.NPROCS-1) THEN
        DO I=0,NKX
            MATD_C(I,NY)=-1.d0
            MATL_C(I,NY)=1.d0
            VEC_C(I,NY) = (0.d0,0.d0)
              VEC_C(I,NY)=VEC_C(I,NY)
     &              -NU*DY(NY-1)*((CU2(I,K,NY)-CU2(I,K,NY-1))/DYF(NY-1)
     &            -(CU2(I,K,NY-1)-CU2(I,K,NY-2))/DYF(NY-2))
     &             /DY(NY-1)
        END DO
        END IF

      RETURN
      END
      
!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_POISSON_ADJ_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
!----*|--.---------.---------.---------.---------.---------.---------.-|------
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      INTEGER I,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANK.eq.0) THEN
        DO I=0,NKX
          IF ((I.EQ.0).AND.(K.EQ.0)) THEN
            MATD_C(I,1)=1.d0
            MATU_C(I,1)=0.d0
            VEC_C(I,1)=(0.d0,0.d0)
          ELSE
! Here, apply Neumann boundary conditions (dp/dz= +-d2y u2/ReT) at the walls
            MATD_C(I,1)=1.d0
            MATU_C(I,1)=-1.d0
            VEC_C(I,1) = (0.d0,0.d0)
!            VEC_C(I,1)=NU*DY(2)*(CU2(I,K,3)-2.d0*CU2(I,K,2)+CU2(I,K,1))
!     &             /(T*DYF(1)**2)
          END IF
        END DO
        END IF
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANK.eq.NPROCS-1) THEN
        DO I=0,NKX
            MATD_C(I,NY)=-1.d0
            MATL_C(I,NY)=1.d0
            VEC_C(I,NY) = (0.d0,0.d0)
!            VEC_C(I,NY)=-NU*DY(N)*(CU2(I,K,N)-2.d0*CU2(I,K,N-1)
!     &       +CU2(I,K,N-2))/(T*DYF(N-1)**2)
        END DO
        END IF

      RETURN
      END      

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_VEL_MPI
!----*|--.---------.---------.---------.---------.---------.---------.-|------
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

! Apply Boundary conditions to velocity field
      IF (RANK.EQ.0) THEN
        CALL APPLY_BC_VEL_LOWER
      END IF
      IF (RANK.EQ.NPROCS-1) THEN
        CALL APPLY_BC_VEL_UPPER
      END IF
 
      RETURN
      END
 
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TRANSPOSE_MPI_XZ_TO_XY(A)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This subroutine starts with all arrays decomposed in x-z slabs
! and transposes the data so that it is decomposed in x-y slabs
! x-y slabs.
      include 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'

      integer i,j,k,N

      real*8 A(0:NX+1,0:NZ+1,0:NY+1) 
      real*8 buffer(0:NX+1,0:NY+1,0:NZ+1)

      real*8 test1(1:NX,1:NY,1:NZ)
      real*8 test2(1:NX,1:NY,1:NZ)

       do i=0,NXM
       do j=1,NY
         write(100+RANK,*) 'i,j,U1: ',i,j,A(I,0,J)
       end do
       end do

      do k=0,NZ-1
        do j=1,NY
          do i=0,NX-1
            buffer(i,j,k)=A(i,k,j)
            test1(i+1,j,k+1)=A(i,k,j)
          end do
        end do
      end do

      CALL MPI_AllToAll(test1,(NX)*(NY)*(NZ)/NPROCS
     &   ,MPI_DOUBLE_PRECISION,test2,(NX)*(NY)*(NZ)/NPROCS
     &       ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
      CALL MPI_Barrier(MPI_COMM_WORLD)

! A should now be indexed by A(NX,NZ/nprocs,NY*nprocs)

         n=1
         do k=0,NZ/NPROCS-1
         do i=0,NX-1 
         do j=1,NY
           A(i,k,j)=test2(i+1,j,k+1)
         end do
         end do
         end do
  
         do n=2,NPROCS
         do k=0,NZ/NPROCS-1
         do i=0,NX-1 
         do j=2,NY
           A(i,k,NY+j-1)=test2(i+1,NY+j,k+1)
         end do
         end do
         end do
         end do

         NX_t=NX
         NZ_t=NZ/NPROCS
         NY_t=NY*NPROCS-(NPROCS-1)
 
! A should now be indexed from A(0:NX_t-1,0:NZ_t/NPROCS-1,1:NY*NPROCS-(NPROCS-1))
! (The new wall locations are 1 and NY-(NPROCS-1) )

          do i=0,NXM
          do j=1,NY_T
            write(200+RANK,*) 'i,j,U1: ',i,j,A(i,0,J)
          end do
          end do


   
      RETURN
      END         
      
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TRAPZ(VECTOR,ANS,GRID)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'
      INCLUDE 'mpif.h'
      INCLUDE 'header_mpi'
      REAL*8, INTENT(IN) :: VECTOR(0:NY+1)
      REAL*8, INTENT(OUT) :: ANS
      INTEGER, INTENT(IN) :: GRID
      INTEGER J

      ANS = 0.d0  
      
      IF (GRID.EQ.1) THEN
        IF (RANK_G.EQ.0) THEN
          DO J = JSTART-1,JEND-1
            ANS = ANS + (GYF(J+1)-GYF(J))*(VECTOR(J+1)+VECTOR(J))/2.0d0
          END DO
        ELSE IF (RANK_G.EQ.NPROCS-1) THEN
          DO J = JSTART-1,JEND
            ANS = ANS + (GYF(J+1)-GYF(J))*(VECTOR(J+1)+VECTOR(J))/2.0d0
          END DO
        ELSE
          DO J = JSTART-1,JEND-1
            ANS = ANS + (GYF(J+1)-GYF(J))*(VECTOR(J+1)+VECTOR(J))/2.0d0
          END DO        
        END IF
      ELSE 
        IF (RANK_G.EQ.0) THEN
! HALF OF VECTOR(2), SINCE ALL VARS ARE ZERO AT WALLS - TRIANGLE AREA        
          ANS = ANS + (GY(2)-GY(1))*VECTOR(2)/4.0d0
          DO J = 2,NY-1
            ANS = ANS + (GY(J+1)-GY(J))*(VECTOR(J+1)+VECTOR(J))/2.0d0
          END DO
        ELSE IF (RANK_G.EQ.NPROCS-1) THEN
! HALF OF VECTOR(NY), SINCE ALL VARS ARE ZERO AT WALLS - TRIANGLE AREA        
          ANS = ANS + (GY(NY+1)-GY(NY))*VECTOR(NY)/4.0d0
          DO J = 1,NY-1
            ANS = ANS + (GY(J+1)-GY(J))*(VECTOR(J+1)+VECTOR(J))/2.0d0
          END DO
        ELSE
          DO J = 1,NY-1
            ANS = ANS + (GY(J+1)-GY(J))*(VECTOR(J+1)+VECTOR(J))/2.0d0
          END DO
        END IF      
      ENDIF                
            
      RETURN
      END                  
      
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE GATHER_SCALAR_MPI(IN,OUT)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'
      INCLUDE 'mpif.h'
      INCLUDE 'header_mpi'
      REAL*8 ICPACK(NPROCS)
      REAL*8, INTENT(IN) :: IN
      REAL*8, INTENT(OUT) :: OUT
      INTEGER J

      CALL MPI_ALLGATHER(IN,1,MPI_DOUBLE_PRECISION,
     &      ICPACK,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
      
      OUT = 0.d0
     
      DO J = 1,NPROCS
        OUT = ICPACK(J) + OUT
      END DO                
            
      RETURN
      END               
      
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE WAIT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'
      INCLUDE 'mpif.h'
      INCLUDE 'header_mpi'

      !IF (RANK_G.EQ.0) THEN
      !WRITE(*,*) 'Waiting...'
      !END IF

      CALL MPI_BARRIER( MPI_COMM_WORLD , IERROR )
      
      !IF (RANK_G.EQ.0) THEN
      !WRITE(*,*) 'Done!'
      !END IF
      
      RETURN
      END   
      
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|      
      SUBROUTINE END_RUN_MPI(FLAG)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
      INCLUDE 'header'
      INCLUDE 'mpif.h'
      INCLUDE 'header_mpi'

      LOGICAL FLAG
      INTEGER ERROR_CODE

      FLAG = .FALSE.

      IF (RANK.EQ.0) THEN
         CALL END_RUN(FLAG)
      END IF

      ! If stop.now exists, then FLAG comes back set to .TRUE.
      IF (FLAG) THEN
         print *, 'FLAG came back true in END_RUN_MPI'
         CALL MPI_ABORT( MPI_COMM_WORLD , ERROR_CODE , IERROR )
      END IF

      !CALL MPI_BCAST(FLAG,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERROR)
      
      END           

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|      
      SUBROUTINE END_RUN_MPI_INIT()
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
      INCLUDE 'header'
      INCLUDE 'mpif.h'
      INCLUDE 'header_mpi'

      LOGICAL FLAG
      INTEGER ERROR_CODE

      IF (RANK.EQ.0) THEN
         CALL END_RUN(FLAG)

         IF (FLAG) THEN
            CALL MPI_ABORT( MPI_COMM_WORLD , ERROR_CODE , IERROR )
         END IF

      END IF
      
      END           


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ET_CHECK_STOP_COND_MPI( ET_COND )
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

        INCLUDE 'header'
        INCLUDE 'mpif.h'
        INCLUDE 'header_mpi'

        LOGICAL :: ET_COND

        IF (RANK == 0) THEN
          CALL ET_CHECK_STOP_COND( ET_COND )
        END IF

        CALL MPI_BCAST(ET_COND,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERROR)
  
      END SUBROUTINE ET_CHECK_STOP_COND_MPI

