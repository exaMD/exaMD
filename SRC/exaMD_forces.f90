!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   exaMD_forces.f90
!!!!!
!!!!!   This file provides the functions for the exaMD proto-app library 
!!!!!   that handle particle interactions between blocks.
!!!!!
!!!!!   The routines contained in this module are:
!!!!!   * exaMD_forces_midpoint_overlap()
!!!!!   * exaMD_forces_midpoint_overlap_bdry()
!!!!!   * exaMD_forces_midpoint_all()
!!!!!   * exaMD_forces_midpoint_all_bdry()
!!!!!   * exaMD_particle_interaction()
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module exaMD_forces_mod

    Use exaMD_mod

    Implicit None

    Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_forces_midpoint_overlap()
!!!!!
!!!!!  * Top level algorithm for computing neighbour interactions.
!!!!!    This routine loops over the boundary blocks calling a separate routine 
!!!!!    to compute the interactions on those. This routine has boundary blocks 
!!!!!    received and processed as the messages arrive.
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_forces_midpoint_overlap(user_func, cutoff, pid, noprocs)

  Use exaMD_comms
  Use exaMD_comms_overlap
  
  Implicit None

  Real(kind=PRC), Intent(in) :: cutoff
  Integer, Intent(in) :: pid, noprocs
  External :: user_func
  
  Integer :: i, j, k, b, i1, i2, j1, j2, k1, k2, pair, natm1, natm2, a1, a2, chunk
  Logical :: flag  ! Test for all mesages having been received
  Logical :: sentflag  ! Test for all mesages having been sent
  Type(pframe_t), Pointer :: block
  Integer, Pointer, Dimension(:,:) :: refs
  Integer, Allocatable :: arrived(:), justin(:)
  Integer :: n_justarrived, n_arrived

  Allocate(arrived(n_extneigh), justin(n_extneigh))
  arrived(:) = 0
  n_arrived = 0

  ! Set up variables for testing computation/communication overlap
  chunk = n_interior/10
  If (n_interior.Lt.10) chunk = n_interior

  sentflag = .False.
  flag = .False.

  ! Loop over interior blocks
  Do b = 1, n_interior

     If (b .Eq. chunk) Then
        If (.Not.sentflag) Then
           Call MPI_Waitall(n_extneigh, reqarray(n_extneigh:2*n_extneigh-1), MPI_STATUSES_IGNORE, ierr)
           sentflag = .True.
        End If

        ! Test to see if messages have all been received
        If (.Not.flag) Then
           n_justarrived = exaMD_testsome_comms_posn(arrived, justin, noprocs)
           If (n_justarrived.Gt.0) Then
              n_arrived = n_arrived + n_justarrived
              Call exaMD_forces_midpoint_overlap_bdry(user_func, cutoff, arrived, n_justarrived, justin, pid, noprocs)
              If (n_arrived.Ge.n_extneigh) Then
                 flag = .True.
              End If
           End If
           chunk = chunk + n_interior/10  ! Do test after every 10%
        End If
     End If

     i = interior_blks(1,b)
     j = interior_blks(2,b)
     k = interior_blks(3,b)
     block => grid(i,j,k)
     refs=> block%neighrefs
     ! Do local interactions on same block first
     natm1 = block%natms
     If (natm1.Gt.1) Then
        Do a1 = 1, natm1-1
           Do a2 = a1+1, natm1
              Call exaMD_particle_interaction(user_func,i,j,k,i,j,k,a1,a2)
           End Do  
        End Do  
     Endif
     
     ! Now do all other interactions
     Do pair = 2, n_interactions
       ! Work out global numbering using pointer to neighrefs in correct grid
        i1 = refs(1,interactions(1,pair))
        i2 = refs(1,interactions(2,pair))
        j1 = refs(2,interactions(1,pair))
        j2 = refs(2,interactions(2,pair))
        k1 = refs(3,interactions(1,pair))
        k2 = refs(3,interactions(2,pair))
        
        natm1 = grid(i1,j1,k1)%natms
        natm2 = grid(i2,j2,k2)%natms

        If (natm1.Ge.1.And.natm2.Ge.1) Then
           Do a1 = 1, natm1
              Do a2 = 1, natm2
                 Call exaMD_particle_interaction(user_func, i1,j1,k1,i2,j2,k2,a1,a2)
              End Do  
           End Do
        Endif  
     End Do  
  End Do

  ! If not all communications have been received, wait for them and then do the boundary points
  If (.Not.flag) Then
     n_justarrived = n_extneigh - Sum(arrived(:))
     Call exaMD_waitrest_comms_posn(arrived, justin, noprocs)
     Call exaMD_forces_midpoint_overlap_bdry(user_func, cutoff, arrived, n_justarrived, justin, pid, noprocs)
  Endif

  Deallocate(arrived, justin)

End Subroutine exaMD_forces_midpoint_overlap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_forces_midpoint_overlap_bdry()
!!!!!
!!!!!  * Called from exaMD_forces_midpoint_overlap this routine calculates the forces 
!!!!!    with home blocks on the boundary between processors.  This operates 
!!!!!    only on messages received, hence will typically be called several 
!!!!!    times.
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * block%flag - that this block has been calculated 
!!!!!  * n_frcsendcnt(i) - count of blocks completed to send to process i
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_forces_midpoint_overlap_bdry(user_func, cutoff, arrived, n_justarrived, justin, pid, noprocs)

  Use exaMD_comms
  Use exaMD_comms_overlap

  Implicit None
  
  Real(kind=PRC),Intent(in) :: cutoff
  Integer,Intent(in) :: pid, noprocs, n_justarrived
  Integer,Allocatable :: arrived(:), justin(:)
  External :: user_func
  
  Integer :: i, j, k, jj, b, i1, i2, j1, j2, k1, k2, pair, natm1, natm2, a1, a2, p, pe, pp, locsendcnt
  Integer :: jcnt, bcnt
  Logical :: skipflag
  Type(pframe_t),Pointer :: block
  Integer,Pointer,Dimension(:,:) :: refs
  Integer,Allocatable :: locsendlist(:)

  If (n_bdry.Eq.0) Return
  
  ! Loop over the processes from which messages have arrived
  Do jcnt = 1, n_justarrived
     
     ! Loop over the boundary blocks related to the processes from which messages have arrived
     Do bcnt = n_pbdryblk_lst(justin(jcnt)), n_pbdryblk_lst(justin(jcnt)+1)-1
        b = pbdryblk_lst(bcnt)
        i = bdry_blks(1,b)
        j = bdry_blks(2,b)
        k = bdry_blks(3,b)
        block => grid(i,j,k)

        ! Have we already calculated this block?
        If (block%flag) Cycle

        ! Loop over connections to see if we've already received all the neighbours to this block
        skipflag = .False.
        Do p = 1, block%nsendtoremote
           ! Do reverse lookup on where in ext_neigh the required processor is
           pe = -1
           Do jj = 1, n_extneigh
              If (block%sendtoremote(p).Gt.-1.And.ext_neigh(jj).Eq.block%sendtoremote(p)) Then
                 pe = jj
                 Exit
              End If
           End Do
           If (pe.Eq.-1) Print *,'ERROR at exaMD_forces.f90 reverse lookup'
           If (1.Ne.arrived(pe)) Then
              skipflag = .True.
              Exit
           Endif
        End Do
        If (skipflag) Cycle  ! If not all ready, move to next block

        ! Do local interactions on same block first
        refs=> block%neighrefs
        natm1 = block%natms
        If (natm1.Gt.1) Then
           Do a1 = 1, natm1-1
              Do a2 = a1+1, natm1
                 Call exaMD_particle_interaction(user_func, i,j,k,i,j,k,a1,a2)
              End Do  
           End Do  
        Endif

        ! Now do all other interactions
        Do pair = 2, n_interactions
           ! Work out global numbering using pointer to neighrefs in correct grid
           i1 = refs(1,interactions(1,pair))
           i2 = refs(1,interactions(2,pair))
           j1 = refs(2,interactions(1,pair))
           j2 = refs(2,interactions(2,pair))
           k1 = refs(3,interactions(1,pair))
           k2 = refs(3,interactions(2,pair))
        
           natm1 = grid(i1,j1,k1)%natms
           natm2 = grid(i2,j2,k2)%natms

           If (natm1.Ge.1.And.natm2.Ge.1) Then
              Do a1 = 1, natm1
                 Do a2 = 1, natm2
                    Call exaMD_particle_interaction(user_func, i1,j1,k1,i2,j2,k2,a1,a2)
                 End Do  
              End Do
           Endif  
        End Do  

        ! Mark this block as done
        block%flag = .True.
     
        Do p = 1, block%nsendtoremote
           pe = block%sendtoremote(p)
           If (pe.Gt.-1.And.pe.Ne.pid) Then
              Do pp = 1, n_extneigh
                 If (ext_neigh(pp).Eq.pe) Then
                    n_frcsendcnt(pp) = n_frcsendcnt(pp) + 1
                    Exit
                 Endif
              End Do
           Endif
        End Do
     
     End Do
  End Do
  
  Allocate(locsendlist(n_extneigh))
  locsendcnt = 0

  ! Now see if we are ready to pass all the forces back to any of the neighbours   
  Do p = 1, n_extneigh
     If (n_frcsendcnt(p).Eq.n_frcsendtgt(p)) Then
        n_frcsendcnt(p) = -1
        locsendcnt = locsendcnt + 1
        locsendlist(locsendcnt) = p
     Endif
  End Do

  ! set up array ready for next iteration
  If (locsendcnt.Gt.0) Then
     Call exaMD_start_comms_forcessome(locsendcnt, locsendlist, noprocs)
  Endif
  Deallocate(locsendlist)

End Subroutine exaMD_forces_midpoint_overlap_bdry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_forces_midpoint_all()
!!!!!
!!!!!  * Top level algorithm for computing neighbour interactions.
!!!!!    This routine loops over the boundary blocks calling a separate routine 
!!!!!    to compute the interactions on those.  This routine has boundary blocks 
!!!!!    received and processed only once all messages have arrived.
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_forces_midpoint_all(user_func, cutoff, noprocs)

  Use exaMD_comms

  Implicit None
  
  Real(kind=PRC),Intent(in) :: cutoff
  Integer, Intent(in) :: noprocs

  Integer, Allocatable :: arrived(:), justin(:)
  External :: user_func
  
  Integer :: i, j, k, b, i1, i2, j1, j2, k1, k2, pair, natm1, natm2, a1, a2, chunk
  Logical :: flag
  Type(pframe_t), Pointer :: block
  Integer, Pointer, Dimension(:,:) :: refs
  
  ! Set up variables for testing computation/communication overlap
  chunk = n_interior/10
  If (n_interior.Lt.10) chunk = n_interior
  flag = .False.

  ! Loop over interior blocks
  Do b = 1, n_interior

     ! Test to see if messages have all been received
     If (b.Eq.chunk.And.(.Not.flag)) Then
        flag = exaMD_test_comms_posn(noprocs)
        If (flag) Then
           Call exaMD_forces_midpoint_all_bdry(user_func, cutoff, noprocs)
        End If
        chunk = chunk + n_interior/10
     Endif

     i = interior_blks(1,b)
     j = interior_blks(2,b)
     k = interior_blks(3,b)
     block => grid(i,j,k)
     refs=> block%neighrefs

     ! Do local interactions on same block first
     natm1 = block%natms
     If (natm1.Gt.1) Then
        Do a1 = 1, natm1-1
           Do a2 = a1+1, natm1
              Call exaMD_particle_interaction(user_func,i,j,k,i,j,k,a1,a2)
           End Do  
        End Do  
     Endif
     
     ! Now do all other interactions
     Do pair = 2, n_interactions
       ! Work out global numbering using pointer to neighrefs in correct grid
        i1 = refs(1,interactions(1,pair))
        i2 = refs(1,interactions(2,pair))
        j1 = refs(2,interactions(1,pair))
        j2 = refs(2,interactions(2,pair))
        k1 = refs(3,interactions(1,pair))
        k2 = refs(3,interactions(2,pair))
        
        natm1 = grid(i1,j1,k1)%natms
        natm2 = grid(i2,j2,k2)%natms

        If (natm1.Ge.1.And.natm2.Ge.1) Then
           Do a1 = 1, natm1
              Do a2 = 1, natm2
                 Call exaMD_particle_interaction(user_func,i1,j1,k1,i2,j2,k2,a1,a2)
              End Do  
           End Do
        Endif  
     End Do  
  End Do

  ! If not all communications have been received then let's wait for them and finally do the boundary points
  If (.Not.flag) Then

     Call exaMD_wait_comms_posn(noprocs)

     Call exaMD_forces_midpoint_all_bdry(user_func, cutoff, noprocs)

  Endif

End Subroutine exaMD_forces_midpoint_all

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_forces_midpoint_all_bdry()
!!!!!
!!!!!  * Called from exaMD_forces_midpoint_all this routine calculates the forces 
!!!!!    with home blocks on the boundary between processors.  This operates 
!!!!!    only once all messages have been received, hence will only be called  
!!!!!    once per iteration.
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_forces_midpoint_all_bdry(user_func, cutoff, noprocs)

  Use exaMD_comms
  
  Implicit None
  
  Real(kind=PRC),Intent(in) :: cutoff
  Integer,Intent(in) :: noprocs
  External :: user_func
  
  Integer :: i, j, k, b, i1, i2, j1, j2, k1, k2, pair, natm1, natm2, a1, a2
  Logical :: flag
  Type(pframe_t),Pointer :: block
  Integer,Pointer,Dimension(:,:) :: refs
  
  ! Receive from other processors 
  Call exaMD_recv_comms_posn(noprocs)

  ! Loop over boundary blocks
  Do b = 1, n_bdry
     i = bdry_blks(1,b)
     j = bdry_blks(2,b)
     k = bdry_blks(3,b)
     block => grid(i,j,k)
     refs=> block%neighrefs
     ! Do local interactions on same block first
     natm1 = block%natms
     If (natm1.Gt.1) Then
        Do a1 = 1, natm1-1
           Do a2 = a1+1, natm1
              Call exaMD_particle_interaction(user_func, i,j,k,i,j,k,a1,a2)
           End Do  
        End Do  
     Endif
     
     ! Now do all other interactions
     Do pair = 2, n_interactions
       ! Work out global numbering using pointer to neighrefs in correct grid
        i1 = refs(1,interactions(1,pair))
        i2 = refs(1,interactions(2,pair))
        j1 = refs(2,interactions(1,pair))
        j2 = refs(2,interactions(2,pair))
        k1 = refs(3,interactions(1,pair))
        k2 = refs(3,interactions(2,pair))
        
        natm1 = grid(i1,j1,k1)%natms
        natm2 = grid(i2,j2,k2)%natms

        If (natm1.Ge.1.And.natm2.Ge.1) Then
           Do a1 = 1, natm1
              Do a2 = 1, natm2
                 Call exaMD_particle_interaction(user_func, i1,j1,k1,i2,j2,k2,a1,a2)
              End Do  
           End Do
        Endif  
     End Do  
  End Do
  
  ! Now pass forces back to neighbours   
  Call exaMD_start_comms_forces(noprocs)

End Subroutine exaMD_forces_midpoint_all_bdry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_particle_interaction()
!!!!!
!!!!!  * Computes the distance between two particles and calls the
!!!!!    user-supplied function to calculate the force between them.
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * force contributions on each of the specified particles
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!DEC$ ATTRIBUTES FORCEINLINE :: exaMD_particle_interaction
Subroutine exaMD_particle_interaction(user_func,i1,j1,k1,i2,j2,k2,a1,a2)

  Implicit None
  
  Integer,Intent(in) :: i1,j1,k1,i2,j2,k2,a1,a2
  External :: user_func

  Integer :: ityp, jtyp, kk, keyfce
  Integer :: offset, retval
  Real(kind=PRC) :: dx, dy, dz, bt(3), d2
  Real(kind=PRC) :: f12(3), e12
  Type(pframe_t),Pointer :: block1
  Type(pframe_t),Pointer :: block2

  block1 => grid(i1,j1,k1)
  block2 => grid(i2,j2,k2)
  ! Calculate distance between points
  dx = block2%posn(1,a2) - block1%posn(1,a1)
  dy = block2%posn(2,a2) - block1%posn(2,a1)
  dz = block2%posn(3,a2) - block1%posn(3,a1)

  ! Fix periodic points
  bt(1) = Nint(dx*globalboxi(1))
  bt(2) = Nint(dy*globalboxi(2))
  bt(3) = Nint(dz*globalboxi(3))

  dx = dx - bt(1)*globalbox(1)
  dy = dy - bt(2)*globalbox(2)
  dz = dz - bt(3)*globalbox(3)

  ! Calculate distance^2 between particles
  d2 = dx**2 + dy**2 + dz**2

  If (Abs(d2).Lt.0.001) Then
     Print'(A,4(I3,X,I3,4X))','ZERO at ',a1,a2,i1,i2,j1,j2,k1,k2
  Endif

  ! Call user supplied function that handles all interactions based on distance alone, supplied via dx, dy and dz
  ! Outputs are the f12 array and a real number e12 along with a return value to indicate problems or methods to use

  Call user_func(dx, dy, dz, d2, f12, e12, retval)

  ! Handle return codes from user supplied subroutine:
  ! 0 - update forces
  ! >0 - do nothing
  ! <0 - error
  If (retval.Eq.0) Then

! 2-d as arrays
     block1%forces(1:3,a1) = block1%forces(1:3,a1) + f12(1:3)
     block2%forces(1:3,a2) = block2%forces(1:3,a2) - f12(1:3)

  Else If (retval.Lt.0) Then
     Print*,'Error from user_func: ', retval

!  else ! retval was positive, hence no update

  Endif

End Subroutine exaMD_particle_interaction

End Module exaMD_forces_mod
