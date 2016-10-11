!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   exaMD_input_formats.f90
!!!!!
!!!!!   This file provides the routines for the exaMD proto-app library 
!!!!!   that handle particle input and output.
!!!!!
!!!!!   The routines contained in this module are:
!!!!!   * exaMD_particles_generate()
!!!!!   * exaMD_read_params()
!!!!!   * exaMD_read_particles()
!!!!!   * exaMD_write_particles()
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module exaMD_input_formats
  
  Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_particles_generate()
!!!!!
!!!!!  * Randomly generates particle locations within the specified domain.
!!!!!
!!!!!  * Assignment done into existing block structure
!!!!!
!!!!!  * Particles tested to ensure they are not too close 
!!!!!    to other particles already placed.
!!!!!
!!!!!  Returns:
!!!!!  * particle_positions - x,y,z particle co-ordinates
!!!!!  * ntot - the number of particles in each block
!!!!!  * first_particle - the indices of the first particle in each block
!!!!!  * next_particle - the indices of the next particle in the same block
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_particles_generate(particle_positions, natms_tot, ntot, first_particle, &
                                    next_particle, nblocks, nblocks_xyz, block_length)

  Use exaMD_mod

  Implicit None
  
  Real(KIND=PRC), Intent(INOUT), Allocatable :: particle_positions(:,:)
  Integer, Intent(INOUT), Allocatable :: ntot(:), first_particle(:), next_particle(:)
  Integer, Intent(IN) :: natms_tot, nblocks, nblocks_xyz(3)
  Real(KIND=PRC), Intent(IN) :: block_length

  Real(KIND=PRC) :: rnd(3), posxyz(3)
  Integer, Allocatable :: last_particle(:)
  Integer :: nblocks_xyplane, blockxyz(3), blocknum, i

  Allocate(last_particle(nblocks))

  ntot(:) = 0
  first_particle(:) = 0
  next_particle(:) = 0
  last_particle(:) = 0

  nblocks_xyplane = nblocks_xyz(1)*nblocks_xyz(2)

  Do i=1, natms_tot

     ! Generate a random particle position within the domain.
1    Call Random_number(rnd)
     posxyz = rnd * globalbox

     ! Identify the block this particle is in.
     Call exaMD_identify_block(posxyz, blockxyz, blocknum, nblocks, nblocks_xyz, nblocks_xyplane, block_length)

     ! Avoid having atoms too close to each other.
     If (exaMD_check_nearby_blocks(ntot, first_particle, next_particle, posxyz, particle_positions, &
                                   blockxyz, blocknum, nblocks_xyz, nblocks_xyplane) .Eq. 1) Goto 1

     ! Record the details of this particle.
     ntot(blocknum) = ntot(blocknum) + 1

     If (ntot(blocknum) .Eq. 1) Then
        first_particle(blocknum) = i
     Else
        next_particle(last_particle(blocknum)) = i
     End If
     last_particle(blocknum) = i
     particle_positions(:,i) = posxyz

     If(Mod(i,10000)==0) Print *,'Number of particles generated:', i

  End Do

  Deallocate(last_particle)

End Subroutine exaMD_particles_generate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_read_params()
!!!!!
!!!!!  * Read in the values of natms_tot and globalbox.
!!!!!
!!!!!  Returns:
!!!!!  * natms_tot - the total number of particles in the simulation
!!!!!
!!!!!  Sets the global variable:
!!!!!  * globalbox - the dimensions of the simulation box
!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_read_params(natms_tot)

  Use exaMD_mod

  Implicit None

  Integer, Intent(INOUT) :: natms_tot
  
  Open(50, file='./exaMD_particles.dat', status='old')
  Read(50,*) natms_tot
  Read(50,*) globalbox

End Subroutine exaMD_read_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_read_particles()
!!!!!
!!!!!  * Loads the set of particle positions from the file exaMD_particles.dat
!!!!!    and places them within the block structure. 
!!!!!
!!!!!  Returns:
!!!!!  * particle_positions - x,y,z particle co-ordinates
!!!!!  * ntot - the number of particles in each block
!!!!!  * first_particle - the indices of the first particle in each block
!!!!!  * next_particle - the indices of the next particle in the same block
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_read_particles(particle_positions, natms_tot, ntot, first_particle, &
                                next_particle, nblocks, nblocks_xyz, block_length)

  Use exaMD_mod

  Implicit None

  Real(KIND=PRC), Intent(INOUT), Allocatable :: particle_positions(:,:)
  Integer, Intent(INOUT), Allocatable :: ntot(:), first_particle(:), next_particle(:)
  Integer, Intent(IN) :: natms_tot, nblocks, nblocks_xyz(3)
  Real(KIND=PRC), Intent(IN) :: block_length

  Integer :: i, nblocks_xyplane, blockxyz(3), blocknum
  Integer, Allocatable :: last_particle(:)
  Real(KIND=PRC) :: posxyz(3), temp(3)
  
  Allocate(last_particle(nblocks))

  ntot(:) = 0
  first_particle(:) = 0
  next_particle(:) = 0
  last_particle(:) = 0

  nblocks_xyplane = nblocks_xyz(1)*nblocks_xyz(2)

  Do i = 1, natms_tot

     ! Read the position of a particle
     Read(50,*) posxyz

     ! Identify the block this particle is in.
     Call exaMD_identify_block(posxyz, blockxyz, blocknum, nblocks, nblocks_xyz, nblocks_xyplane, block_length)

     ! N.B. The particles were checked for closeness when they were initially generated.

     ! Record the details of this particle.
     ntot(blocknum) = ntot(blocknum) + 1

     If (ntot(blocknum) .Eq. 1) Then
        first_particle(blocknum) = i
     Else
        next_particle(last_particle(blocknum)) = i
     End If
     last_particle(blocknum) = i
     particle_positions(:,i) = posxyz

     If(Mod(i,10000)==0) Print *,'Number of particles loaded:', i

  End Do
  
  Close(50)
  
  Deallocate(last_particle)

End Subroutine exaMD_read_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_write_particles()
!!!!!
!!!!!  * Writes out the list of particle positions to the file 
!!!!!    exaMD_particles.dat for subsequent re-reading
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_write_particles(particle_positions, natms_tot)

  Use exaMD_mod

  Implicit None

  Real(KIND=PRC), Intent(IN), Allocatable :: particle_positions(:,:)
  Integer, Intent(IN) :: natms_tot

  Integer :: i
  
  Open(50, file='exaMD_particles.dat', status='replace')
  Write(50,*) natms_tot
  Write(50,*) globalbox

  Do i = 1, natms_tot
     Write(50,*) particle_positions(1:3,i)
  End Do

  Close(50)

End Subroutine exaMD_write_particles

End Module exaMD_input_formats
