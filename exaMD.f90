!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   exaMD.f90
!!!!!
!!!!!   This is an example driving program for using the exaMD proto-app 
!!!!!   library.
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program exaMD

  Use exaMD_mod
  Use exaMD_partitioning
  Use exaMD_comms
  Use exaMD_forces_mod
  Use exaMD_input_formats
  Use exaMD_utils_mod
  Use user_func_mod

  Implicit None

  Real(kind=PRC) :: cutoff, block_length      ! Problem parameters
  Real(kind=PRC) :: rnd(3), ftot(3), ftot2(3) ! Arrays for position, forces sum and forces^2 sum
  Real(kind=Kind(0.0d0)) :: t0, t1, t2, t3    ! MPI timer values
  Real(kind=Kind(0.0d0)) :: total_time_and_rank(2), max_time_and_rank(2)
  Integer :: natms_tot                        ! global number of atoms
  Integer, Allocatable :: ntot(:)             ! Number of atoms per block
  Integer, Allocatable :: first_particle(:)   ! Indices of the first particle in each block
  Integer, Allocatable :: next_particle(:)    ! Indices of the next particle in the same block
  Real(kind=PRC), Allocatable :: particle_positions(:,:)
  Integer :: i, j, k, ii, jj, b               ! Loop indices
  Integer :: nblocks                          ! Total number of blocks
  Integer :: nblocks_xyz(3)                   ! Numbers of blocks in each dimension
  Integer :: comm=MPI_COMM_WORLD              ! MPI communicator
  Integer :: noprocs                          ! Number of MPI processes
  Integer :: pid                              ! MPI rank
  Integer :: argumentcount                    ! Number of command line arguments
  Integer :: AllocateStatus                   ! Status returned by Allocate
  Integer, Allocatable :: part(:,:,:)
  Character(len=32) :: arg                    ! String for each input argument
  Logical :: inputcase, overlapcase, writeout ! Variables for input arguments

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Initialization
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Call exaMD_Init(pid, noprocs, comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Default problem specification
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  natms_tot = 500000 ! global number of particles (default value)
  globalbox = (/Real(100.0, KIND=PRC), Real(100.0, KIND=PRC), Real(100.0, KIND=PRC)/) ! global box size
  cutoff = Real(8.0, KIND=PRC) ! interaction cut-off distance, R_c

  block_length = half*cutoff

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Process command line arguments
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  inputcase = .False.
  overlapcase = .False.
  writeout = .False.

  argumentcount = command_argument_count()
  i = 1
  Do While (i.Le.argumentcount)
     Call get_command_argument(i, arg)
     If (arg .Eq. "-overlap") overlapcase = .True.
     If (arg .Eq. "-readfile") inputcase = .True.
     If (arg .Eq. "-writefile") writeout = .True.
     If (arg .Eq. "-atoms") Then
        ! read the number of atoms
        i = i + 1
        Call get_command_argument(i, arg)
        Read( arg, '(i32)' ) natms_tot
     End If
     i = i + 1
  End Do

  ! If we are reading the particle data from file, 
  ! then read in natms_tot and globalbox now.
  If (inputcase) Then
     If (pid==0) Call exaMD_read_params(natms_tot)
     ! Note that only process 0 requires natms_tot but all processes need globalbox
     Call MPI_Bcast(globalbox, 3, MYFLOAT, 0, COMM_EXAMD, ierr)
  End If

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Domain creation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Call exaMD_domain_create(block_length, nblocks, nblocks_xyz)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Output the main simulation parameters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  If(pid==0) Then
     Write (*,*) '-----------------------------------------------------------------------'
     Write (*,*) 'The exaMD Proto-app Example Program'
     Write (*,*) '-----------------------------------'
     Write (*,*)
     Write (*,'(A,I0)') ' Number of atoms: ', natms_tot
     Write (*,'(A,3F7.1)') ' Simulation box size:', globalbox
     Write (*,'(A,F6.2)') ' Cut-off distance:', cutoff
     Write (*,'(A,3(X,I0))') ' Number of blocks in each dimension:', nblocks_xyz
     Write (*,'(A,I0)') ' Number of MPI processes: ', noprocs 
     Write (*,'(A,I0)') ' Floating-point precision: ', PRC
     If(overlapcase) Then
        Write (*,*) 'Overlap option: selected'
     Else
        Write (*,*) 'Overlap option: not selected'
     End If
     Write (*,*) '-----------------------------------------------------------------------'
  End If

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Generate particle positions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Allocate ntot() on all processes.
  ! The array holds the number of particles in each block.
  Allocate(ntot(nblocks), stat=AllocateStatus)
  If (AllocateStatus /= 0) Then
     Write(*,*) "Insufficient memory to allocate array ntot"
     Call MPI_Abort(COMM_EXAMD, 108, ierr)
  End If

  ! Generate particles on root processor
  If(pid==0) Then
     ! Allocate temporary memory for holding particle position information
     Allocate(particle_positions(3,natms_tot), stat=AllocateStatus)
     If (AllocateStatus /= 0) Then
        Write(*,*) "Insufficient memory to allocate array particle_positions"
        Call MPI_Abort(COMM_EXAMD, 103, ierr)
     End If
     Allocate(first_particle(nblocks), stat=AllocateStatus)
     If (AllocateStatus /= 0) Then
        Write(*,*) "Insufficient memory to allocate array first_particle"
        Call MPI_Abort(COMM_EXAMD, 104, ierr)
     End If
     Allocate(next_particle(natms_tot), stat=AllocateStatus)
     If (AllocateStatus /= 0) Then
        Write(*,*) "Insufficient memory to allocate array next_particle"
        Call MPI_Abort(COMM_EXAMD, 105, ierr)
     End If

     If (inputcase) Then
        Call exaMD_read_particles(particle_positions, natms_tot, ntot, first_particle, &
                                  next_particle, nblocks, nblocks_xyz, block_length)
     Else
        Call exaMD_particles_generate(particle_positions, natms_tot, ntot, first_particle, &
                                      next_particle, nblocks, nblocks_xyz, block_length)
     End If

     If (writeout) Then
        Call exaMD_write_particles(particle_positions, natms_tot)
     End If
     
  End If

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Parallel partitioning
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Allocate memory to hold the partition information
  Allocate(part(nblocks_xyz(1),nblocks_xyz(2),nblocks_xyz(3)), stat=AllocateStatus)
  If (AllocateStatus /= 0) Then
     Write(*,*) "Insufficient memory to allocate array part"
     Call MPI_Abort(COMM_EXAMD, 106, ierr)
  End If

  ! Do parallel partitioning on process 0
  If (pid.Eq.0) Then
     Call exaMD_partition_domain(nblocks, nblocks_xyz, ntot, pid, noprocs, part)
  Endif

  ! Broadcast the partitions to all processes
  Call MPI_Bcast(part, nblocks, MPI_INTEGER, 0, COMM_EXAMD, ierr)

  ! Set up the lists of interior, boundary and exterior blocks
  Call exaMD_setup_partitioned_grid(part, nblocks, nblocks_xyz, pid)

  Deallocate(part)

  ! Distribute work to processes
  Call exaMD_partition_work(particle_positions, ntot, first_particle, next_particle, nblocks, nblocks_xyz, pid, noprocs)

  If (pid==0) Deallocate(particle_positions, first_particle, next_particle)
  Deallocate(ntot)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  User set-up section
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! chemical types and force tables
  Allocate(kktype(1,1), ltpvdw(1), prmvdw(1,5))
  kktype = 1
  ltpvdw = 2
  prmvdw(1,1:5) = (/one, one, zero, zero, cutoff/) ! force parameters

  ! Everything from this point on is done on "own" block set

  Do b = 1, nblocks_loc
     i = locblocks(1,b)
     j = locblocks(2,b)
     k = locblocks(3,b)
     grid(i,j,k)%ltype(:) = 1
  End Do
  Do b = 1, n_exterior
     i = exterior_blks(1,b)
     j = exterior_blks(2,b)
     k = exterior_blks(3,b)
     If(grid(i,j,k)%natms.Gt.0) Then
        grid(i,j,k)%ltype(:) = 1
     End If
  End Do

!!! Structure of the rest:
!!! 
!!! Start communications
!!! Do work when relevant comms arrived
!!! Send results back and combine
!!! Final results
!!! 

  ! Set-up pairwise block interactions for Midpoint method
  Call exaMD_setup_midpoint_interactions()

  ! Synchronise processes before starting timing
  Call MPI_Barrier(COMM_EXAMD, ierr)
  t0 = MPI_WTIME()

  If (noprocs>1) Then
     ! Initiate the communication of particle positions to neighbouring domains
     Call exaMD_start_comms_posn()
  Endif

  t1 = MPI_WTIME()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Compute interactions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  If (overlapcase) Then
     Call exaMD_forces_midpoint_overlap(lennardjones_12, cutoff, pid, noprocs)
  Else
     Call exaMD_forces_midpoint_all(lennardjones_12, cutoff, noprocs)
  Endif
  
  t2 = MPI_WTIME()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Collect forces back
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  If (noprocs.Gt.1) Then
     Call exaMD_wait_comms_forces(noprocs)
  End If

  t3 = MPI_WTIME()
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Check and report on validity of calculated forces 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  jj = 0
  Do b=1, nblocks_loc
     i = locblocks(1,b)
     j = locblocks(2,b)
     k = locblocks(3,b)
     Do ii = 1, grid(i,j,k)%natms
        If(grid(i,j,k)%forces(1,ii)/=grid(i,j,k)%forces(1,ii)) &
             Print *, 'NaN for fxx of atom',ii,'block',b,'on task',pid
        If(grid(i,j,k)%forces(2,ii)/=grid(i,j,k)%forces(2,ii)) &
             Print *, 'NaN for fyy of atom',ii,'block',b,'on task',pid
        If(grid(i,j,k)%forces(3,ii)/=grid(i,j,k)%forces(3,ii)) Then
           Print *, 'NaN for fzz of atom',ii,'block',b,'on task',pid
           Print *, grid(i,j,k)%posn(1:3,:)
           jj = jj + 1
        End If
     End Do
     If (jj>10) Exit
  End Do

  ! Check results
  rnd = zero
  Do b=1, nblocks_loc
     i = locblocks(1,b)
     j = locblocks(2,b)
     k = locblocks(3,b)
     rnd = rnd + (/ Sum(grid(i,j,k)%forces(1,1:grid(i,j,k)%natms)), &
                    Sum(grid(i,j,k)%forces(2,1:grid(i,j,k)%natms)), &
                    Sum(grid(i,j,k)%forces(3,1:grid(i,j,k)%natms)) /)
  End Do
  Call MPI_REDUCE(rnd,ftot,3,MYFLOAT,MPI_SUM,0,COMM_EXAMD,ierr)

  rnd = zero
  Do b=1, nblocks_loc
     i = locblocks(1,b)
     j = locblocks(2,b)
     k = locblocks(3,b)
     rnd = rnd + (/ Sum(grid(i,j,k)%forces(1,1:grid(i,j,k)%natms)**2), &
                    Sum(grid(i,j,k)%forces(2,1:grid(i,j,k)%natms)**2), &
                    Sum(grid(i,j,k)%forces(3,1:grid(i,j,k)%natms)**2) /)
  End Do
  Call MPI_REDUCE(rnd,ftot2,3,MYFLOAT,MPI_SUM,0,COMM_EXAMD,ierr)

  If(pid == 0) Then
     Write (*,*) 'Checks on force calculations'
     If (Abs(ftot(1)).Gt.1e-10) Write(*,*) 'WARNING! Suspected error'
     Write (*,'(A,3E10.1)') '  SUM(fxx, fyy, fzz):', ftot
     Write (*,'(A,3E10.1)') '  SUM(fxx^2, fyy^2, fzz^2):', ftot2
     Write(*,*)
  End If

  ! Find the process that took the maximum total time
  total_time_and_rank(1) = t3-t0
  total_time_and_rank(2) = Real(pid, KIND=Kind(0.0d0))
  Call MPI_ALLREDUCE(total_time_and_rank, max_time_and_rank, 1, MPI_2DOUBLE_PRECISION, &
       MPI_MAXLOC, COMM_EXAMD, ierr)

  ! Print the times out from this process
  If(pid == Int(max_time_and_rank(2))) Then
     Write(*,'(A,I0)') ' Times from process ', Int(max_time_and_rank(2))
     Write(*,'(2(A,F10.4))') '  Communication times (secs):  Pre', t1-t0,';  Post ', t3-t2
     Write(*,'(A,F10.4)') '  Computation time (secs):', t2-t1
     Write(*,'(A,F10.4)') '  Total time (secs):', t3-t0
  End If

  ! Barrier needed to make sure the above appears before any subsequent output from process 0
  Call MPI_BARRIER(COMM_EXAMD, ierr)

!  write(*,'(I3,2(A,F10.4))') pid,' Communication times (secs):  Pre', t1-t0,';  Post ', t3-t2
!  write(*,'(I3,A,F10.4)') pid,' Computation time (secs):', t2-t1
!  write(*,'(I3,A,F10.4)') pid,' Total time (secs):', t3-t0

#ifdef PRINTALLFORCES
  Do ii=0, noprocs-1
     If (ii.Eq.pid) Then
        Do b=1, nblocks_loc
           i = locblocks(1,b)
           j = locblocks(2,b)
           k = locblocks(3,b)
           rnd = (/ Sum(grid(i,j,k)%forces(1,1:grid(i,j,k)%natms)**2), &
                    Sum(grid(i,j,k)%forces(2,1:grid(i,j,k)%natms)**2), &
                    Sum(grid(i,j,k)%forces(3,1:grid(i,j,k)%natms)**2) /)
           If (Sum(rnd(1:3)) .Gt. zero) &
                Write(*,'(3I4.3,3(X,E20.12))') i,j,k, rnd(1:3)
        End Do
     End If
     Call FLUSH(6)
     Call MPI_Barrier(COMM_EXAMD, ierr)
  End Do
#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Program clean-up
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Free user memory used
  Deallocate(kktype, ltpvdw, prmvdw)

  Call exaMD_Finalize('Done')

End Program exaMD
