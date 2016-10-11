!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   exaMD_comms_overlap.f90
!!!!!
!!!!!   This file provides the functions for the exaMD proto-app library 
!!!!!   that handle overlapping communication of subsets of the boundary and 
!!!!!   exterior block lists for positional and calculated forces.
!!!!!
!!!!!   The routines contained in this module are:
!!!!!   * exaMD_testsome_comms_posn()
!!!!!   * exaMD_recvsome_comms_posn()
!!!!!   * exaMD_waitrest_comms_posn()
!!!!!   * exaMD_start_comms_forcessome()
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module exaMD_comms_overlap
  
  Contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_testsome_comms_posn()
!!!!!
!!!!!  * Tests to see how many messages from neighbouring processors with force 
!!!!!    calculations have been received.  Calls exaMD_recvsome_comms_posn()
!!!!!    to actually unpack the received data.
!!!!!
!!!!!  Returns:
!!!!!  * exaMD_testsome_comms_posn - number of messages newly received
!!!!!  * justin - an array of the indices into the ext_neigh list of the newly 
!!!!!             arrived messages
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Integer Function exaMD_testsome_comms_posn(arrived, justin, noprocs)

  Use exaMD_mod

  Implicit None

  Integer, Intent(IN) :: noprocs
  Integer, Allocatable, Intent(INOUT) :: arrived(:), justin(:)

  Integer :: p

  ! Return if running in serial
  If (noprocs.Eq.1) Then
     exaMD_testsome_comms_posn = 1
     Return
  End If
  
  Call MPI_Testsome(n_extneigh, reqarray(0:n_extneigh-1), exaMD_testsome_comms_posn, justin, MPI_STATUSES_IGNORE, ierr)

  If (exaMD_testsome_comms_posn.Gt.0) Then
     ! Mark read messages as such in arrived array - needed for subsequent calls
     Do p = 1, exaMD_testsome_comms_posn
        arrived(justin(p)) = 1
     End Do
  End If

  ! Copy the data from the completed messages
  Call exaMD_recvsome_comms_posn(justin, exaMD_testsome_comms_posn, noprocs)

End Function exaMD_testsome_comms_posn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_recvsome_comms_posn()
!!!!!
!!!!!  * Unpacks the messages received from processes listed in the justin list
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * block%posn data from the receiving processes
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_recvsome_comms_posn(justin, n_justin, noprocs)

  Use exaMD_mod

  Implicit None

  Integer, Intent(IN) :: n_justin, noprocs
  Integer, Allocatable, Intent(IN) :: justin(:)

  Integer :: i, ii, jj, kk, b, p, pp, jcnt, bcnt
  Logical :: testflag
  Type(pframe_t), Pointer :: block

  If (noprocs.Eq.1) Return

  ! Now unpack received messages
  datasize_r(:)=0
  ! do b = 1, n_exterior
  ! Loop over the processes from which messages have arrived
  Do jcnt = 1, n_justin
     
     ! Loop over the boundary blocks related to the processes from which messages have arrived
     Do bcnt = n_pextblk_lst(justin(jcnt)), n_pextblk_lst(justin(jcnt)+1)-1
        b = pextblk_lst(bcnt)
        ii = exterior_blks(1,b)
        jj = exterior_blks(2,b)
        kk = exterior_blks(3,b)
        block => grid(ii,jj,kk)
        p = block%host

        If (p.Gt.-1) Then
           Do pp = 1, n_extneigh
              If (p.Eq.ext_neigh(pp)) Exit
           End Do
           ! check if this process's message has just arrived
           testflag = .False.
           Do i = 1, n_justin
              If (justin(i).Eq.pp) Then
                 testflag = .True.
                 Exit
              Endif
           End Do
           If (.Not.testflag) Cycle

           If (block%natms.Gt.0) Then
              block%posn(1:3,1:block%natms) = recvbuffer(1:3, datasize_r(p):datasize_r(p)+block%natms-1, pp)
              ! Zero forces array for calculation
              block%forces(:,:) = zero
           Endif
           datasize_r(p) = datasize_r(p) + block%natms  
        Endif

     End Do
  End Do

End Subroutine exaMD_recvsome_comms_posn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_waitrest_comms_posn()
!!!!!
!!!!!  * At end of interior block loop there is no point in waiting for any 
!!!!!    individual messages, so may as well wait for them all.
!!!!!
!!!!!  Returns:
!!!!!  * justin - array of processes handled in this clean-up section
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_waitrest_comms_posn(arrived, justin, noprocs)

  Use exaMD_mod

  Implicit None

  Integer, Intent(IN) :: noprocs
  Integer, Allocatable :: justin(:), arrived(:)

  Integer :: p, n_justin

  ! Mark the remaining messages as received and put entries into justin array ready for the actual unpacking
  n_justin = 0
  Do p = 1, n_extneigh
     If (arrived(p).Eq.0) Then
        Call MPI_Wait(reqarray(p-1), MPI_STATUS_IGNORE, ierr)
        n_justin = n_justin + 1
        justin(n_justin) = p
        arrived(p) = 1
     Endif
  End Do

  ! Actually do the receives from these processors 
  Call exaMD_recvsome_comms_posn(justin, n_justin, noprocs)

End Subroutine exaMD_waitrest_comms_posn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_start_comms_forcessome()
!!!!!
!!!!!  * Sends out blocks to particular processors who have already had all the 
!!!!!    relevant contributions calculated
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * frecvbuffer - buffer used for sending forces
!!!!!  * reqarray - array of MPI requests
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_start_comms_forcessome(sendcnt, sendlist, noprocs)

  Use exaMD_mod

  Implicit None

  Integer, Intent(IN) :: noprocs, sendcnt
  Integer, Allocatable, Intent(IN) :: sendlist(:)

  Integer :: i, ii, jj, kk, b, p, pp, pe, jcnt, bcnt
  Type(pframe_t),Pointer :: block
  Logical :: locflag

  ! This routine is based on a copy of the posn one, however the changes are:
  ! * sendbuffer and recvbuffer are swapped as this time we are receiving from many rather than sending to many
  !   + this means moving the double loop to the receiving phase
  !   + when the results are received they are added to the forces array rather than replacing them

  If (noprocs.Eq.1) Return

  ! Fill recvbuffer with forces to send

  ! Loop over the processes from which messages have arrived
  Do jcnt = 1, sendcnt
     
     ! Loop over the boundary blocks related to the processes from which messages have arrived
     Do bcnt = n_pextblk_lst(sendlist(jcnt)), n_pextblk_lst(sendlist(jcnt)+1)-1
        b = pextblk_lst(bcnt)
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
        
           ! See if this is a process to send to
           locflag = .False.
           Do pe = 1, sendcnt
              If (sendlist(pe).Eq.pp) Then
                 locflag = .True.
                 Exit
              Endif
           End Do
        
           ! If so add it into the send buffer
           If (locflag) Then
              frecvbuffer(1:3,datasize_fr(p):datasize_fr(p)+block%natms-1,pp) = block%forces(1:3,1:block%natms)
              datasize_fr(p) = datasize_fr(p) + block%natms  
           Endif
        Endif
     End Do
  End Do

  ! Initiate MPI_Isends
  Do b = 1, sendcnt
     p = ext_neigh(sendlist(b))

     If (p.Eq.-1) Print*,'why p=-1 in ext_neigh? (start_comms_forcessome)'
     Call MPI_Isend(frecvbuffer(1,0,sendlist(b)), 3*datasize_fr(p), MYFLOAT, p, 300, COMM_EXAMD, &
                    reqarray(3*n_extneigh+sendlist(b)-1), ierr)
  End Do  

End Subroutine exaMD_start_comms_forcessome

End Module exaMD_comms_overlap
