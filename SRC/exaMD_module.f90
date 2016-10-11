!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   exaMD_module.f90
!!!!!
!!!!!   This file provides the main data module for the exaMD proto-app library 
!!!!!   that all other routines reference. The full block structure and lists 
!!!!!   are all contained here.
!!!!!
!!!!!   The routines contained in this module are:
!!!!!   * exaMD_identify_block()
!!!!!   * exaMD_distance2()
!!!!!   * exaMD_check_nearby_blocks()
!!!!!   * exaMD_check_block()
!!!!!   * exaMD_setup_midpoint_interactions()
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module exaMD_mod

  Use exaMD_constants_mod
  Use MPI

  Implicit None

  Integer :: MYFLOAT  ! Precision of floating point arithmetic to work in
  
  Integer :: COMM_EXAMD  ! MPI communicator to use in this library

  Real(KIND=PRC), Parameter :: PROXIMITY=one  ! initial minimum separation of particles 

  Integer, Allocatable :: natms(:) ! number of atoms on task
  Integer, Allocatable :: ltype(:,:) ! chemical species identifier

  ! struct for the information stored about every block
  Type pframe_t
    ! Variables describing the box
    Integer :: host                        ! the rank of the process which hosts this block
    Integer, Allocatable :: neighhosts(:)  ! the hosts of the 26 neighbouring blocks
    Integer, Allocatable :: neighrefs(:,:) ! i, j and k references of the 26 neighbouring blocks

    ! Variables describing the atoms in the box
    Integer :: natms
    Real(kind=PRC), Allocatable :: posn(:,:) ! atomic coords
    Real(kind=PRC), Allocatable :: forces(:,:) ! array of forces
    Integer, Allocatable :: ltype(:) ! chemical species identifier
    Logical :: flag ! Flag to indicate whether interactions with neighbours have been computed 
    Integer :: nsendtoremote ! Flag indicating how many processes this block should be sent to
    Integer, Allocatable :: sendtoremote(:) ! Array indicating which processes this block should be sent to
    Integer :: getfromremote   ! Flag to indicate that this block should be being received
    Integer :: partype ! 0=remote, 1=interior, 2=boundary, 3=off-process but neighbouring
  End Type pframe_t

  Type(pframe_t), Allocatable, Target  :: grid(:,:,:) ! The full domain of blocks
  
  Real(kind=PRC) :: globalbox(3), globalboxi(3) ! global domain size and its inverse

  ! Lists of blocks on this processor
  Integer :: nblocks_loc, n_interior, n_bdry, n_exterior, n_extneigh  ! Counts of blocks in each list
  Integer, Allocatable :: locblocks(:,:)     ! List of block references on this processor (3-d 1st index)
  Integer, Allocatable :: interior_blks(:,:) ! List of blocks entirely inside this processor's partition
  Integer, Allocatable :: bdry_blks(:,:)     ! List of the blocks that form this processor's boundary
  Integer, Allocatable :: exterior_blks(:,:) ! List of the blocks that boundary this processor
  Integer, Allocatable :: ext_neigh(:)       ! List of the processors that boundary this processor
  Integer, Allocatable :: n_frcsendtgt(:), n_frcsendcnt(:)  ! List of the counts of blocks needed by each proc for forces received
  Integer, Allocatable :: n_pbdryblk_lst(:)  ! Counts of boundary blocks with neighbours on each neighbouring process
  Integer, Allocatable :: pbdryblk_lst(:)    ! List of boundary blocks with neighbours on each process 
  Integer, Allocatable :: n_pextblk_lst(:)   ! Counts of exterior blocks with neighbours on each neighbouring process
  Integer, Allocatable :: pextblk_lst(:)     ! List of exterior blocks with neighbours on each process 

  ! MPI arrays - data sizes, buffers, MPI requests and statuses
  Integer, Allocatable :: datasize_r(:), datasize_s(:)
  Integer, Allocatable :: datasize_fr(:), datasize_fs(:)
  Real(kind=PRC), Allocatable :: recvbuffer(:,:,:), sendbuffer(:,:,:)
  Real(kind=PRC), Allocatable :: frecvbuffer(:,:,:), fsendbuffer(:,:,:)
  Integer, Allocatable :: irecvbuffer(:,:), isendbuffer(:,:)
  Integer, Allocatable :: reqarray(:)
  
  ! Error code for MPI routines
  Integer :: ierr

  ! Interactions for midpoint method
  Integer :: n_interactions
  Integer, Allocatable, Dimension(:,:) :: interactions

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_identify_block()
!!!!!
!!!!!  * Identifies the block from a given position
!!!!!
!!!!!  Returns:
!!!!!  * blockxyz - the x,y,z co-ordinates of the block
!!!!!  * blocknum - the global index of the block
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_identify_block(posxyz, blockxyz, blocknum, nblocks, nblocks_xyz, nblocks_xyplane, block_length)

  Implicit None

  Real(KIND=PRC), Intent(IN) :: posxyz(3), block_length
  Integer, Intent(IN) :: nblocks, nblocks_xyz(3), nblocks_xyplane
  Integer, Intent(OUT) :: blockxyz(3), blocknum

  blockxyz = Int(posxyz/block_length + one)

  ! Allow for cases where the end block is bigger,
  ! due to the global length not being a multiple of the block_length.
  Where (blockxyz > nblocks_xyz) blockxyz = nblocks_xyz

  blocknum = blockxyz(1) + (blockxyz(2)-1)*nblocks_xyz(1) &
                         + (blockxyz(3)-1)*nblocks_xyplane

  If (blocknum<1 .Or. blocknum>nblocks) Then
     Print *, 'Block number out of range: ', blocknum
     Print *, 'Block: ', blockxyz
     Print *, 'Position: ', posxyz
     Call MPI_Abort(COMM_EXAMD, 109, ierr)
  End If

End Subroutine exaMD_identify_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_distance2()
!!!!!
!!!!!  * Calculates the square of the Euclidean distance between two points
!!!!!
!!!!!  Returns:
!!!!!  * exaMD_distance2 - the square of the distance
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Real(kind=PRC) Function exaMD_distance2(pos1, pos2)

  Implicit None

  Real(kind=PRC), Dimension(3), Intent(in) :: pos1, pos2

  Real(kind=PRC), Dimension(3) :: disp ! displacement

  disp = pos1 - pos2

  ! Take the periodic boundary conditions into account
  disp = disp - Anint(disp*globalboxi) * globalbox

  exaMD_distance2 = Dot_product(disp, disp)

End Function exaMD_distance2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_check_nearby_blocks()
!!!!!
!!!!!  * Checks all neighbouring blocks to see if a specified point can be 
!!!!!    added to the list of placed particles.
!!!!!
!!!!!  Returns:
!!!!!  * exaMD_check_nearby_blocks - 1 if a particle in one of the blocks is too close
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Integer Function exaMD_check_nearby_blocks(ntot, first_particle, next_particle, posxyz, particle_positions, &
                                           blockxyz, blocknum, nblocks_xyz, nblocks_xyplane)

  Implicit None

  Integer, Intent(in)                     :: blockxyz(3), blocknum, nblocks_xyz(3), nblocks_xyplane
  Integer, Allocatable, Intent(in)        :: ntot(:), first_particle(:), next_particle(:)
  Real(KIND=PRC), Intent(in)              :: posxyz(3)
  Real(KIND=PRC), Allocatable, Intent(in) :: particle_positions(:,:)

  Integer :: ii, jj, kk, i, j, k, blk, num_blocks_x, num_blocks_y, num_blocks_z, partial_sum1, partial_sum2
  
  ! Start by testing own block  
  exaMD_check_nearby_blocks = 0  
  If (exaMD_check_block(ntot, first_particle, next_particle, posxyz, particle_positions, blocknum).Eq.1) Then
     exaMD_check_nearby_blocks = 1
     Return 
  End If
  
  num_blocks_x = nblocks_xyz(1)
  num_blocks_y = nblocks_xyz(2)
  num_blocks_z = nblocks_xyz(3)

  ! Now test neighbouring blocks
  ii = blockxyz(1) + num_blocks_x - 1
  jj = blockxyz(2) + num_blocks_y - 1
  kk = blockxyz(3) + num_blocks_z - 1
  Do k = -1, 1
     partial_sum1 = Mod(kk+k, num_blocks_z) * nblocks_xyplane + 1
     Do j = -1, 1
        partial_sum2 = Mod(jj+j, num_blocks_y) * num_blocks_x + partial_sum1
        Do i = -1, 1
           If (i.Eq.0 .And. j.Eq.0 .And. k.Eq.0) Cycle ! Own block already tested
           blk = Mod(ii+i, num_blocks_x) + partial_sum2
           If (exaMD_check_block(ntot, first_particle, next_particle, posxyz, particle_positions, blk).Eq.1) Then
              exaMD_check_nearby_blocks = 1
              Return
           End If
        End Do  
     End Do  
  End Do  

End Function exaMD_check_nearby_blocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_check_block()
!!!!!
!!!!!  * A potential particle position is tested against
!!!!!    every particle on the specified test block to see if any particles fall
!!!!!    within a distance PROXIMITY.
!!!!!
!!!!!  Returns:
!!!!!  * exaMD_check_block - 1 if a particle is too close
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Integer Function exaMD_check_block(ntot,first_particle,next_particle,posxyz,particle_positions,blk)

  Implicit None

  Integer, Intent(in)                     :: blk
  Integer, Allocatable, Intent(in)        :: ntot(:), first_particle(:), next_particle(:)
  Real(KIND=PRC), Intent(in)              :: posxyz(3)
  Real(KIND=PRC), Allocatable, Intent(in) :: particle_positions(:,:)

  Integer :: i, particle_index

  exaMD_check_block = 0  
  particle_index = first_particle(blk)
  Do i = 1, ntot(blk)
     If (exaMD_distance2(posxyz,particle_positions(:,particle_index)) < PROXIMITY) Then
        exaMD_check_block = 1
        Exit
     End If
     particle_index = next_particle(particle_index)
  End Do

  Return

End Function exaMD_check_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_setup_midpoint_interactions()
!!!!!
!!!!!  * Sets up the interactions list specifying which blocks interactions
!!!!!    are computed for each home block.
!!!!!
!!!!!  * 63 interaction pairs have been identified.  Numbering in comments 
!!!!!    denotes numbers from the user guide
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * interactions array - used for connections in the midpoint method
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_setup_midpoint_interactions()

  Implicit None

  n_interactions = 63
  Allocate(interactions(8,n_interactions))
  interactions(1:2, 1) = (/ 14, 14  /) ! 5 <-> 5
  interactions(1:2, 2) = (/ 14, 15  /) ! 5 <-> 6
  interactions(1:2, 3) = (/ 14, 16  /) ! 5 <-> 7
  interactions(1:2, 4) = (/ 14, 17  /) ! 5 <-> 8
  interactions(1:2, 5) = (/ 14, 18  /) ! 5 <-> 9
  interactions(1:2, 6) = (/ 14, 19  /) ! 5 <-> 10
  interactions(1:2, 7) = (/ 14, 20  /) ! 5 <-> 11
  interactions(1:2, 8) = (/ 14, 21  /) ! 5 <-> 12
  interactions(1:2, 9) = (/ 14, 22  /) ! 5 <-> 13
  interactions(1:2,10) = (/ 14, 23  /) ! 5 <-> 14
  interactions(1:2,11) = (/ 14, 24  /) ! 5 <-> 15
  interactions(1:2,12) = (/ 14, 25  /) ! 5 <-> 16
  interactions(1:2,13) = (/ 14, 26  /) ! 5 <-> 17
  interactions(1:2,14) = (/ 14, 27  /) ! 5 <-> 18
  interactions(1:2,15) = (/ 10, 18  /) ! 1 <-> 9
  interactions(1:2,16) = (/ 10, 27  /) ! 1 <-> 18
  interactions(1:2,17) = (/ 11, 17  /) ! 2 <-> 8
  interactions(1:2,18) = (/ 11, 18  /) ! 2 <-> 9
  interactions(1:2,19) = (/ 11, 25  /) ! 2 <-> 16
  interactions(1:2,20) = (/ 11, 26  /) ! 2 <-> 17
  interactions(1:2,21) = (/ 11, 27  /) ! 2 <-> 18
  interactions(1:2,22) = (/ 12, 16  /) ! 3 <-> 7
  interactions(1:2,23) = (/ 12, 17  /) ! 3 <-> 8
  interactions(1:2,24) = (/ 12, 25  /) ! 3 <-> 16
  interactions(1:2,25) = (/ 13, 15  /) ! 4 <-> 6
  interactions(1:2,26) = (/ 13, 18  /) ! 4 <-> 9
  interactions(1:2,27) = (/ 13, 24  /) ! 4 <-> 15
  interactions(1:2,28) = (/ 13, 27  /) ! 4 <-> 18
  interactions(1:2,29) = (/ 15, 16  /) ! 6 <-> 7
  interactions(1:2,30) = (/ 15, 19  /) ! 6 <-> 10
  interactions(1:2,31) = (/ 15, 22  /) ! 6 <-> 13
  interactions(1:2,32) = (/ 15, 25  /) ! 6 <-> 16
  interactions(1:2,33) = (/ 16, 21  /) ! 7 <-> 12
  interactions(1:2,34) = (/ 16, 24  /) ! 7 <-> 15
  interactions(1:2,35) = (/ 17, 19  /) ! 8 <-> 10
  interactions(1:2,36) = (/ 17, 20  /) ! 8 <-> 11
  interactions(1:2,37) = (/ 17, 21  /) ! 8 <-> 12
  interactions(1:2,38) = (/ 18, 19  /) ! 9 <-> 10
  interactions(1:2,39) = (/ 19, 5  /) ! 10 <-> 23
  interactions(1:2,40) = (/ 19, 9  /) ! 10 <-> 27
  interactions(1:2,41) = (/ 20, 5  /) ! 11 <-> 23
  interactions(1:2,42) = (/ 20, 7  /) ! 11 <-> 25
  interactions(1:2,43) = (/ 20, 8  /) ! 11 <-> 26
  interactions(1:2,44) = (/ 20, 9  /) ! 11 <-> 27
  interactions(1:2,45) = (/ 21, 5  /) ! 12 <-> 23
  interactions(1:2,46) = (/ 21, 7  /) ! 12 <-> 25
  interactions(1:2,47) = (/ 22, 5  /) ! 13 <-> 23
  interactions(1:2,48) = (/ 22, 6  /) ! 13 <-> 24
  interactions(1:2,49) = (/ 22, 9  /) ! 13 <-> 27
  interactions(1:2,50) = (/ 23, 5  /) ! 14 <-> 23
  interactions(1:2,51) = (/ 24, 5  /) ! 15 <-> 23
  interactions(1:2,52) = (/ 24, 4  /) ! 15 <-> 22
  interactions(1:2,53) = (/ 24, 7  /) ! 15 <-> 25
  interactions(1:2,54) = (/ 25, 5  /) ! 16 <-> 23
  interactions(1:2,55) = (/ 25, 3  /) ! 16 <-> 21
  interactions(1:2,56) = (/ 25, 6  /) ! 16 <-> 24
  interactions(1:2,57) = (/ 26, 5  /) ! 17 <-> 23
  interactions(1:2,58) = (/ 26, 2  /) ! 17 <-> 20
  interactions(1:2,59) = (/ 26, 3  /) ! 17 <-> 21
  interactions(1:2,60) = (/ 27, 1  /) ! 18 <-> 19
  interactions(1:2,61) = (/ 27, 2  /) ! 18 <-> 20
  interactions(1:2,62) = (/ 27, 4  /) ! 18 <-> 22
  interactions(1:2,63) = (/ 27, 5  /) ! 18 <-> 23

End Subroutine exaMD_setup_midpoint_interactions

End Module exaMD_mod
