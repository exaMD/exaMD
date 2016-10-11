!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   exaMD_utils.f90
!!!!!
!!!!!   This file provides the utility functions for the exaMD proto-app library 
!!!!!   for solving molecular dynamics problems using the midpoint method.
!!!!!
!!!!!   The routines contained in this module are:
!!!!!   * exaMD_Init : MPI initialization including communicator duplication
!!!!!   * exaMD_Finalize : does MPI finalization and memory freeing
!!!!!   * exaMD_domain_create : creates appropriate grid structure of blocks 
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module exaMD_utils_mod
  
  Implicit None
  
  Contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_Init()
!!!!!
!!!!!  * Sets precision to be used throughout the library
!!!!!  * Initialises MPI and sets up a local communicator, COMM_EXAMD, 
!!!!!    for use within exaMD.
!!!!!
!!!!!  Returns:
!!!!!  * pid - Process ID in MPI communicator COMM_EXAMD
!!!!!  * noprocs - Number of processes in COMM_EXAMD
!!!!!  
!!!!!  Global structure data set
!!!!!  * MYFLOAT - Precision to use
!!!!!  * COMM_EXAMD - MPI communicator to use within library
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_Init(pid, noprocs, input_comm)

  Use exaMD_mod
  
  Implicit None
  
  Integer, Intent(OUT) :: pid, noprocs
  Integer, Intent(IN) :: input_comm
  Logical :: flag

  ! Precision of calculations
  
  If(PRC == 8) Then
     MYFLOAT = MPI_DOUBLE_PRECISION
  Else If (PRC == 4) Then
     MYFLOAT = MPI_REAL
  Else
     Call exaMD_Finalize('wrong float type')
  End If


  ! Parallel initialisation

  ! First check whether MPI has already been initialised. If not, call MPI_Init.  
  Call MPI_Initialized(flag, ierr)
  If (.Not. flag) Then
     Call MPI_Init(ierr)
  Endif
  
  ! Create a new communicator for use within exaMD
  Call MPI_Comm_dup(input_comm, COMM_EXAMD, ierr)
  
  ! Get the rank of this process and the total number of processes
  Call MPI_COMM_RANK(COMM_EXAMD, pid, ierr)
  Call MPI_COMM_SIZE(COMM_EXAMD, noprocs, ierr)

End Subroutine exaMD_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_Finalize()
!!!!!
!!!!!  * Finalizes MPI and frees memory
!!!!!
!!!!!  Returns:
!!!!!  (nothing)
!!!!!
!!!!!  Global structure data set
!!!!!  * memory freed in many global exaMD arrays
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_Finalize(text)

  Use exaMD_mod
  
  Implicit None
  
  Character(len=*) :: text
  Integer :: pid, noprocs
  Logical :: flag

  Call MPI_Initialized(flag, ierr)
  
  If (flag) Then
     Call MPI_COMM_RANK( COMM_EXAMD, pid, ierr )
     Call MPI_COMM_SIZE( COMM_EXAMD, noprocs, ierr )
  Else
     pid = 0
     noprocs = 1
  End If

  If (pid.Eq.0) Then
     Print *
     Print '(2x,2a)','Stopping.... ', text
  End If

  ! Let's exit as gracefully as we can, clearing up memory as we go

  If (flag) Then
     Call MPI_FINALIZE(ierr)
  End If

  ! Free memory allocated
  If (Allocated(grid)) Deallocate(grid)
  If (Allocated(locblocks)) Deallocate(locblocks)

  If (Allocated(exterior_blks)) Deallocate(exterior_blks)
  If (Allocated(bdry_blks)) Deallocate(bdry_blks)
  If (Allocated(interior_blks)) Deallocate(interior_blks)
  If (Allocated(ext_neigh)) Deallocate(ext_neigh)
  If (noprocs.Gt.1) Then
     If (Allocated(recvbuffer)) Deallocate(recvbuffer)
     If (Allocated(sendbuffer)) Deallocate(sendbuffer)
     If (Allocated(frecvbuffer)) Deallocate(frecvbuffer)
     If (Allocated(fsendbuffer)) Deallocate(fsendbuffer)
     If (Allocated(reqarray)) Deallocate(reqarray)
     If (Allocated(n_frcsendcnt)) Deallocate(n_frcsendcnt)
     If (Allocated(n_frcsendtgt)) Deallocate(n_frcsendtgt)
     If (Allocated(n_pbdryblk_lst)) Deallocate(n_pbdryblk_lst)
     If (Allocated(pbdryblk_lst)) Deallocate(pbdryblk_lst)
     If (Allocated(n_pextblk_lst)) Deallocate(n_pextblk_lst)
     If (Allocated(pextblk_lst)) Deallocate(pextblk_lst)
  End If
  If (Allocated(datasize_s)) Deallocate(datasize_s)
  If (Allocated(datasize_r)) Deallocate(datasize_r)
  If (Allocated(datasize_fs)) Deallocate(datasize_fs)
  If (Allocated(datasize_fr)) Deallocate(datasize_fr)
  If (Allocated(interactions)) Deallocate(interactions)

  Call Exit(0)

End Subroutine exaMD_Finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_domain_create()
!!!!!
!!!!!  * Allocates global array of blocks that are used to describe the domain.
!!!!!
!!!!!  Returns:
!!!!!  (nothing)
!!!!!
!!!!!  Global structure data set
!!!!!  * globalboxi set to 1.0/globalbox
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Subroutine exaMD_domain_create(block_length, nblocks, nblocks_xyz)

  Use exaMD_mod
  
  Implicit None
  
  Real(kind=PRC), Intent(IN) :: block_length ! R_c/2
  Integer, Intent(OUT) :: nblocks, nblocks_xyz(3)

  Integer :: AllocateStatus ! Status returned by Allocate

  ! simple cubic splitting
  nblocks_xyz = Int(globalbox/block_length)
  nblocks = Product(nblocks_xyz)

  Allocate(grid(nblocks_xyz(1),nblocks_xyz(2),nblocks_xyz(3)), stat=AllocateStatus)
  If (AllocateStatus /= 0) Then
     Write(*,*) "Insufficient memory to allocate array grid"
     Call MPI_Abort(COMM_EXAMD, 107, ierr)
  End If

  globalboxi = one / globalbox

End Subroutine exaMD_domain_create

End Module exaMD_utils_mod
