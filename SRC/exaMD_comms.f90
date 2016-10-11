!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   exaMD_comms.f90
!!!!!
!!!!!   This file provides the functions for the exaMD proto-app library 
!!!!!   that handle communication of full sets of the boundary and 
!!!!!   exterior block lists for positions and calculated forces.
!!!!!
!!!!!   The routines contained in this module are:
!!!!!   * exaMD_start_comms_posn
!!!!!   * exaMD_test_comms_posn
!!!!!   * exaMD_wait_comms_posn
!!!!!   * exaMD_recv_comms_posn
!!!!!   * exaMD_start_comms_forces
!!!!!   * exaMD_wait_comms_forces
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module exaMD_comms
  
  Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_start_comms_posn()
!!!!!
!!!!!  * Starts communication of particle positions on boundary blocks to all 
!!!!!    neighbouring processes that need this information.  Buffers declared 
!!!!!    and filled for the sending part, and appropriately sized for receiving 
!!!!!    the forces.  MPI_Isend and Irecvs for these messages posted, along with 
!!!!!    Irecvs for the final forces back too.
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * sendbuffer - array for the position data to be sent to each process
!!!!!  * recvbuffer - array for the position data to be received on each process
!!!!!  * fsendbuffer - array for the force data to be received to each process
!!!!!  * frecvbuffer - array for the force data to be sent on each process
!!!!!  * reqarray - array of MPI requests
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_start_comms_posn()

  Use exaMD_mod

  Implicit None

  Integer :: i, ii, jj, kk, b, p, pp, reqcnt
  Type(pframe_t), Pointer :: block

  ! Streamlined communication of blocks to other processes can be broken down into:
  ! * Calculate datsizes for send and receives
  ! * Post MPI_Irecvs for all exterior_blks()
  ! * Post MPI_Isends for all bdry_blks()
  ! * At some point before calculation check for appropriate receive completion
  ! * Eventually check for send completion
  ! Note : if to be run on very large processor counts need to avoid noprocs length arrays

  ! Fill sendbuffer
  datasize_s(:)=0
  Do b = 1, n_bdry
     ii = bdry_blks(1,b)
     jj = bdry_blks(2,b)
     kk = bdry_blks(3,b)
     block => grid(ii,jj,kk)
     ! Now there may be multiple processes wanting this same block
     Do i = 1, block%nsendtoremote
        p = block%sendtoremote(i)
        If (p.Gt.-1) Then
           Do pp = 1, n_extneigh
              If (p.Eq.ext_neigh(pp)) Exit
           End Do
           datasize_s(p) = datasize_s(p) + block%natms  
        Endif
     End Do
  End Do
  ! Allocate memory for sending and receiving arrays
  Allocate(sendbuffer(1:3,0:Maxval(datasize_s(:)),1:n_extneigh))
  Allocate(recvbuffer(1:3,0:Maxval(datasize_r(:)),1:n_extneigh))
  Allocate(reqarray(0:4*n_extneigh-1))
  Allocate(fsendbuffer(1:3,0:Maxval(datasize_s(:)),1:n_extneigh))
  Allocate(frecvbuffer(1:3,0:Maxval(datasize_r(:)),1:n_extneigh))
  frecvbuffer(:,:,:) = zero
  fsendbuffer(:,:,:) = zero
  
  datasize_s(:)=0
  Do b = 1, n_bdry
     ii = bdry_blks(1,b)
     jj = bdry_blks(2,b)
     kk = bdry_blks(3,b)
     block => grid(ii,jj,kk)
     ! Now there may be multiple processes wanting this same block
     Do i = 1, block%nsendtoremote
        p = block%sendtoremote(i)
        If (p.Gt.-1) Then
           Do pp = 1, n_extneigh
              If (p.Eq.ext_neigh(pp)) Exit
           End Do
           If (pp.Eq.n_extneigh+1) Then
              Print*,'ERROR - pp not found (2)'
           Endif
           sendbuffer(1:3,datasize_s(p):datasize_s(p)+block%natms-1,pp) = block%posn(1:3,1:block%natms)
           datasize_s(p) = datasize_s(p) + block%natms  
        Endif
     End Do
  End Do
  

  ! Initiate MPI_Isends and Irecvs
  reqcnt = 0
  Do b = 1, n_extneigh
     p = ext_neigh(b)
     If (p.Eq.-1) Print*,'Help! -1 in ext_neigh'
     Call MPI_Irecv(recvbuffer(1,0,b), 3*datasize_r(p), MYFLOAT, p, 200, COMM_EXAMD, reqarray(reqcnt), ierr)
     Call MPI_Isend(sendbuffer(1,0,b), 3*datasize_s(p), MYFLOAT, p, 200, COMM_EXAMD, reqarray(reqcnt+n_extneigh), ierr)

     ! Set up receive for forces back too
     Call MPI_Irecv(fsendbuffer(1,0,b), 3*datasize_fs(p), MYFLOAT, p, 300, COMM_EXAMD, &
                    reqarray(2*n_extneigh+reqcnt), ierr)
     reqcnt = reqcnt + 1
  End Do  

End Subroutine exaMD_start_comms_posn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_test_comms_posn()
!!!!!
!!!!!  * Tests for all incoming position messages having been received
!!!!!
!!!!!  Returns:
!!!!!  * exaMD_test_comms_posn - true/false, depending on whether they've been received or not
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Logical Function exaMD_test_comms_posn(noprocs)

  Use exaMD_mod

  Implicit None

  Integer, Intent(IN) :: noprocs

  If (noprocs.Eq.1) Then
     exaMD_test_comms_posn = .True.
     Return
  Endif
  
  Call MPI_Testall(n_extneigh, reqarray(0:n_extneigh-1), exaMD_test_comms_posn, MPI_STATUSES_IGNORE, ierr)

End Function exaMD_test_comms_posn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_wait_comms_posn()
!!!!!
!!!!!  * Waits for all incoming position messages to have been received
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_wait_comms_posn(noprocs)

  Use exaMD_mod

  Implicit None

  Integer, Intent(IN) :: noprocs

  If (noprocs.Gt.1) &
     Call MPI_Waitall(n_extneigh, reqarray(0:n_extneigh-1), MPI_STATUSES_IGNORE, ierr)

End Subroutine exaMD_wait_comms_posn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_recv_comms_posn()
!!!!!
!!!!!  * Unpacks the messages received containing position information
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * block%posn data from the receiving processes
!!!!!  * block%forces - zeroed
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_recv_comms_posn(noprocs)

  Use exaMD_mod

  Implicit None

  Integer, Intent(IN) :: noprocs

  Integer :: ii, jj, kk, b, p, pp
  Type(pframe_t),Pointer :: block

  If (noprocs.Eq.1) Return

  ! Now unpack received messages
  datasize_r(:)=0
  Do b = 1, n_exterior
     ii = exterior_blks(1,b)
     jj = exterior_blks(2,b)
     kk = exterior_blks(3,b)
     block => grid(ii,jj,kk)
     p = block%host
     If (p.Gt.-1) Then
        Do pp = 1, n_extneigh
           If (p.Eq.ext_neigh(pp)) Exit
        End Do
        If (block%natms.Gt.0) Then
           block%posn(1:3,1:block%natms) = recvbuffer(1:3, datasize_r(p):datasize_r(p)+block%natms-1, pp)
           ! Zero forces array for calculation
           block%forces(:,:) = zero
        Endif
        datasize_r(p) = datasize_r(p) + block%natms  
     Endif

  End Do

  ! This is only testing for the sends as receives already done
  Call MPI_Waitall(n_extneigh, reqarray(n_extneigh:2*n_extneigh-1), MPI_STATUSES_IGNORE, ierr)

End Subroutine exaMD_recv_comms_posn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_start_comms_forces()
!!!!!
!!!!!  * Sends out contributions to forces on all exterior blocks  
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * frecvbuffer - buffer used for sending forces
!!!!!  * reqarray - array of MPI requests
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_start_comms_forces(noprocs)

  Use exaMD_mod

  Implicit None

  Integer, Intent(IN) :: noprocs

  Integer :: ii, jj, kk, b, p, pp, reqcnt
  Type(pframe_t),Pointer :: block

  ! This routine is based on a copy of the posn one above, however the changes are:
  ! * sendbuffer and recvbuffer are swapped as this time we are receiving from many rather than sending to many
  !   + this means moving the double loop to the receiving phase
  !   + when the results are received they are added to the forces array rather than replacing them

  If (noprocs.Eq.1) Return

  ! Fill recvbuffer with forces to send
  Do b = 1, n_exterior
     ii = exterior_blks(1,b)
     jj = exterior_blks(2,b)
     kk = exterior_blks(3,b)
     block => grid(ii,jj,kk)

     ! Only one process owns this block, hence wants the results
     p = block%host
     If (p.Gt.-1) Then
        Do pp = 1, n_extneigh
           If (p.Eq.ext_neigh(pp)) Exit
        End Do
        frecvbuffer(1:3,datasize_fr(p):datasize_fr(p)+block%natms-1,pp) = block%forces(1:3,1:block%natms)
        datasize_fr(p) = datasize_fr(p) + block%natms  
     Endif
  End Do
  
  ! Initiate MPI_Isends and Irecvs
  reqcnt = 0
  Do b = 1, n_extneigh
     p = ext_neigh(b)
     If (p.Eq.-1) Print*,'why p=-1 in ext_neigh?'
     Call MPI_Isend(frecvbuffer(1,0,b), 3*datasize_fr(p), MYFLOAT, p, 300, COMM_EXAMD, &
                    reqarray(3*n_extneigh+reqcnt), ierr)
     reqcnt = reqcnt + 1
  End Do  

End Subroutine exaMD_start_comms_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_wait_comms_forces()
!!!!!
!!!!!  * This function waits for all forces to have been received from all 
!!!!!    neighbouring processes since there is no point in continuing until
!!!!!    we have full information on all boundary blocks.
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * block%forces data for all boundary blocks
!!!!!  * datasize_fs - data reference through received array from each process
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_wait_comms_forces(noprocs)

  Use exaMD_mod

  Implicit None

  Integer, Intent(IN) :: noprocs

  Integer :: i, ii, jj, kk, b, p, pp
  Type(pframe_t),Pointer :: block

  ! This is only testing for the receives - sends tested later
  Call MPI_Waitall(n_extneigh, reqarray(2*n_extneigh:3*n_extneigh-1), MPI_STATUSES_IGNORE, ierr)

  ! Now unpack received messages
  datasize_fs(:)=0
  Do b = 1, n_bdry
     ii = bdry_blks(1,b)
     jj = bdry_blks(2,b)
     kk = bdry_blks(3,b)
     block => grid(ii,jj,kk)
     
     ! Could be receiving contributions from multiple
     Do i = 1, block%nsendtoremote
        p = block%sendtoremote(i)
        If (p.Gt.-1) Then
           Do pp = 1, n_extneigh
              If (p.Eq.ext_neigh(pp)) Exit
           End Do
           If (pp.Eq.n_extneigh+1) Then
              Print*,'ERROR - pp not found (5)'
              Call MPI_Abort(COMM_EXAMD, 102, ierr)
           Endif

           block%forces(1:3,1:block%natms) = block%forces(1:3,1:block%natms) + &
                                             fsendbuffer(1:3, datasize_fs(p):datasize_fs(p)+block%natms-1, pp)
           datasize_fs(p) = datasize_fs(p) + block%natms  
        Endif
     End Do
     
  End Do

  ! This is only testing for the sends as receives already done
  Call MPI_Waitall(n_extneigh, reqarray(n_extneigh:2*n_extneigh-1), MPI_STATUSES_IGNORE, ierr)
  Call MPI_Waitall(n_extneigh, reqarray(3*n_extneigh:4*n_extneigh-1), MPI_STATUSES_IGNORE, ierr)

End Subroutine exaMD_wait_comms_forces

End Module exaMD_comms
