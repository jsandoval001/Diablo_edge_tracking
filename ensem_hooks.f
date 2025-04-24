C******************************************************************************|
C these are the hooks required by the 'basic' version of diablo, even
C when not running the 'ensemble' flavor    
C
C******************************************************************************|
	subroutine init_ensem
	include 'header_ensem'
	real version, current_version
	real*8 var,dLX,dLZ,zbqlu01
	integer i,j,k,grid
	
	open(12,file='input_ensem.dat',form='formatted',status='old')  
	
	current_version = 1.0
	read(12,*)
	read(12,*)
	read(12,*)
	read(12,*)
	read(12,*) version
	if (version.ne.current_version) stop 'Wrong input data format.'
	read(12,*)
	read(12,*) scaling
	read(12,*)
	read(12,*) meas_freq, meas_ave, var
	
	close(12)
	
	root = 0	
C Initialize measurements
	call zbqlini(rank)
	meas_freq = 1/(delta_t*meas_freq)
	if (meas_freq.eq.0) then
	 meas_freq = 1
	end if
	
	do i=1,NM
	  meas_var(i) = var
	end do
	
	do k=1,THL
	  do j=1,NV
	    do i=1,3
	      veh_pos(i,j,k) = 0.d0
	    end do
	  end do
	end do
	
	if (rank.eq.root) then
	grid = 4
	do j=0,grid-1
	  do i=0,grid-1
	  veh_pos(1,grid*j+i+1,0) = width*LX/(grid-1)*i
     &				   +(1.-width)*LX/2.
	  veh_pos(3,grid*j+i+1,0) = width*LX/(grid-1)*j
     &				   +(1.-width)*LZ/2.
	  end do
	end do
	end if
	sze = 3*NV*THL
	call mpi_bcast(veh_pos,sze,mpi_real8,root,mpi_comm_world,ierr)
	
	
	meas_type(1) = 'TH'
	meas_type(2) = 'U1'
	meas_type(3) = 'U3'
	
C	v0=0.
	mag_var = 0.
	theta = 2.*pi*zbqlu01(pi)  

	
	return
	end 
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine create_flow_ensem
	include 'header_ensem'
	integer i,k,j
	real*8 zbqlu01
	
	if (num_per_dir.eq.3) then
	  if (rank.eq.root) then
	  call fft_xzy_to_physical(CU3,U3)
	  call fft_xzy_to_physical(CU1,U1)
	  do j=0,nym
	    do k=0,nzm
	      do i=0,nxm
		  U3(i,k,j) = 0.10
		  U1(i,k,j) = 0.20
		end do
	    end do
	  end do
	  call fft_xzy_to_fourier(U3,CU3)
	  call fft_xzy_to_fourier(U1,CU1)
	  else
	  call fft_xzy_to_physical(CU3,U3)
	  call fft_xzy_to_physical(CU1,U1)
	  do j=0,nym
	    do k=0,nzm
	      do i=0,nxm
		  U3(i,k,j) = 0.40*zbqlu01(pi)-0.20
		  U1(i,k,j) = 0.40*zbqlu01(pi)-0.20
		end do
	    end do
	  end do
	  call fft_xzy_to_fourier(U3,CU3)
	  call fft_xzy_to_fourier(U1,CU1)
	  end if
	  
	else
	  write(6,*) 'Non-periodic cases not supported'
	end if
	
	
	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine save_stats_ensem(final)
	include 'header_ensem'
	logical final
	integer i,j,k,n
	
	if (num_per_dir.eq.3) then
	  call average
	  call variance
	  if (rank.eq.root) then
	  do j=0,tnky
	  do k=0,tnkz
	    do i=0,nkx
	      CF2(i,k,j) = ci*(KX(i)*U3(i,k,j)-KZ(k)*U1(i,k,j)) 
		CR2(i,k,j) = ci*(KX(i)*R3(i,k,j)-KZ(k)*R1(i,k,j))
	    end do
	  end do
	end do
	  do n=1,N_TH
	    call fft_xzy_to_physical(CTH(0,0,0,n),TH(0,0,0,n))
	    call fft_xzy_to_physical(CRTH(0,0,0,n),RTH(0,0,0,n))
	  end do
	  call fft_xzy_to_physical(CU1,U1)
	  call fft_xzy_to_physical(CU3,U3)
	  call fft_xzy_to_physical(CR1,R1)
	  call fft_xzy_to_physical(CR3,R3)
	  call fft_xzy_to_physical(CF2,F2)
	  call fft_xzy_to_physical(CR2,R2)
	  
	  open(unit=23,file="truth.dat",status="unknown")
	  open(unit=24,file="mean.dat",status="unknown")
	  open(unit=25,file="phivar.dat",status="unknown")
	  open(unit=26,file="u1var.dat",status="unknown")
	  open(unit=27,file="u3var.dat",status="unknown")
	  open(unit=28,file="u1.dat",status="unknown")
	  open(unit=29,file="u3.dat",status="unknown")
	  open(unit=30,file="u1mean.dat",status="unknown")
	  open(unit=31,file="u3mean.dat",status="unknown")
	  open(unit=32,file="tvort.dat",status="unknown")
	  open(unit=33,file="evort.dat",status="unknown")
	
	  if (final) then
	    write(23,*) 1
	  else
	    write(23,*) 0
	  end if
	  write(23,*) LX
	  write(23,*) LZ
	  write(23,*) width
	  write(23,*) nx
	  write(23,*) nz
	  write(23,*) nv
	  
	  do i=1,nv
	    write(23,*) veh_pos(1,i,0)
	    write(23,*) veh_pos(3,i,0)
	  end do
	  
	  do i=0,nxm
	    write(23,*) gx(i)
	  end do
	  do k=0,nzm
	    write(23,*) gz(k)
	  end do
	  
	  write(23,*) pos(1)
	  write(23,*) pos(3)
	  write(24,*) rpos(1)
	  write(24,*) rpos(3)
	  do n=1,N_TH
	  do k=0,nzm
	    do i=0,nxm
	      write(23,*) TH(i,k,0,n)		
		write(24,*) RTH(i,k,0,n)
		write(25,*) FTH(i,k,0,n)
		write(26,*) F1(i,k,0)
		write(27,*) F3(i,k,0)
		write(28,*) U1(i,k,0)
		write(29,*) U3(i,k,0)
		write(30,*) R1(i,k,0)
		write(31,*) R3(i,k,0)
		write(32,*) F2(i,k,0)
		write(33,*) R2(i,k,0)
	    end do
	  end do
	  end do
	  
	  close(23)
	  close(24)
	  close(25)
	  close(26)
	  close(27)
	  close(28)
	  close(29)
	  close(30)
	  close(31)
	  close(32)
	  close(33)
	  
	  do n=1,N_TH
	    call fft_xzy_to_fourier(TH(0,0,0,n),CTH(0,0,0,n))
	  end do
	  call fft_xzy_to_fourier(U1,CU1)
	  call fft_xzy_to_fourier(U3,CU3)

C	  call system('mv truth.dat ../post_process/matlab/truth.dat',i)
C	  call system('mv mean.dat ../post_process/matlab/mean.dat',i)
        call system('scp *.dat 
     &  dynamics:~/Desktop/diablo/post_process/.',i)
	  call system('ssh dynamics
     &  "mv ~/Desktop/diablo/post_process/*.dat 
     &  ~/Desktop/diablo/post_process/matlab/."',i)
	  end if
	else
	  write(6,*) 'Non-periodic cases not supported'
	end if
	
	
	return
	end