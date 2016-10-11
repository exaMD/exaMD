!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   exaMD_partitioning.f90
!!!!!
!!!!!   This file provides the functions for the exaMD proto-app library 
!!!!!   that handle the partitioning of the domain between processes and
!!!!!   the distribution of data between them. It also sets up the lists 
!!!!!   used in later parts of the code.
!!!!!
!!!!!   The routines contained in this module are:
!!!!!   * exaMD_partition_domain()
!!!!!   * exaMD_print_partition()
!!!!!   * exaMD_find_partition_extent()
!!!!!   * exaMD_partition_receive()
!!!!!   * exaMD_setup_partitioned_grid()
!!!!!   * exaMD_partition_work()
!!!!!   * exaMD_partition_neighbours()
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Partitioning of blocks for arbitrary processor counts
! Choosing to do this via METIS
!   * Using the multilevel k-way partitioning method
!   * the graph really has diagonal communications in there but choosing not to use these as:
!     + diagonals make alignment on diagonals as important as on axial direction (bad for memory allocation)
!     + Seems to keep partitions more contiguous for larger core counts
!     + Can be disabled using -DDODIAGS compiler flag

Module exaMD_partitioning

! Use this option to avoid empty blocks with empty neighbours from being skipped in partitioning
#define NONULLBLOCK 

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_partition_domain()
!!!!!
!!!!!  * This routine is only run on a single processor and sets up the arrays
!!!!!    used in the partitioning of the blocked domain into an arbitrary number 
!!!!!    of partitions. It uses METIS to perform this partitioning, using the 
!!!!!    multilevel k-way partitioning method although other packages could 
!!!!!    easily be substituted in instead with limited changes to the code.
!!!!!    
!!!!!  * The partitioning uses the number of atoms in blocks as weights in the 
!!!!!    partitioning process. Edge weights are formed of the sum of the number 
!!!!!    of particles at each end.  This is an estimate, ignoring the other 
!!!!!    connections the block would need to receive.
!!!!!    
!!!!!  * The true graph would have diagonal connectivities included. However, 
!!!!!    we do not use them since:
!!!!!     + diagonals make alignment on diagonals as important as in axial 
!!!!!       direction (bad for memory allocation)
!!!!!     + Seems to keep partitions more contiguous for larger core counts
!!!!!     + Can be disabled using -DDODIAGS compiler flag
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * most of the block data is set in this routine for the root process
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_partition_domain(nblocks, nblocks_xyz, ntot, pid, noprocs, vtxdist)

  Use exaMD_mod

  Implicit None

  Integer, Intent(IN) :: nblocks, nblocks_xyz(3), pid, noprocs
  Integer, Intent(IN), Allocatable :: ntot(:)
  Integer, Intent(INOUT), Allocatable :: vtxdist(:,:,:)

  Integer, Allocatable :: xadj(:), adjncy(:), vwgt(:), vsize(:), adjwgt(:)
  Integer, Allocatable :: options(:), part(:)
  Integer :: nvtxs, ncon, nparts, objval
  Integer :: num_blocks_x, num_blocks_y, num_blocks_z
  Integer :: i, j, k, im, ip, jm, jp, km, kp, bcnt, cnt, split(3)
  Real(kind=PRC) :: tgtwgt, avgblksInv, avgblksSqInv
  Real(kind=PRC), Allocatable :: tpwgts(:), ubvec(:)

  num_blocks_x = nblocks_xyz(1)
  num_blocks_y = nblocks_xyz(2)
  num_blocks_z = nblocks_xyz(3)

! Work out weight distribution
  If (Allocated(ntot)) Then
     cnt = 0
     Do i = 1, nblocks
        cnt = cnt+ntot(i)
     End Do
     avgblksInv = Real(nblocks,kind=PRC)/(Real(cnt,kind=PRC))
     avgblksSqInv = avgblksInv*avgblksInv
  Else
     Print*,'ERROR : ntot not allocated'
     Call MPI_Abort(COMM_EXAMD, 101, ierr)
  Endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! METIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set basic parameters for partitioning
  nvtxs = nblocks
  ncon = 1
  nparts = noprocs

  ! Number blocks
#ifdef DODIAGS
  cnt=0
  Do k = 1, num_blocks_z
     Do j = 1, num_blocks_y
        Do i = 1, num_blocks_x
           vtxdist(i,j,k) = cnt
           cnt = cnt+1
        End Do
     End Do
  End Do
#endif

  ! Allocate memory for connections in CSR format
  Allocate(xadj(0:nblocks))         ! List of start point in adjncy for connections
#ifdef DODIAGS
  Allocate(adjncy(0:26*nblocks))    ! List of connections for each node
  Allocate(adjwgt(1:26*nblocks))    ! List of edge weights
#else
  Allocate(adjncy(0:6*nblocks))     ! List of connections for each node
  Allocate(adjwgt(1:6*nblocks))     ! List of edge weights
#endif

  !Allocate memory for weights
  Allocate(vwgt(1:nblocks))      ! List of weights of vertices (related to computational cost)

  ! Assign connections
  cnt=0
  bcnt=0
#ifdef DODIAGS
  Do i = 0, nblocks
     xadj(i) = i*26
  End Do
#else
!  do i = 0, nblocks
!     xadj(i) = i*6
!     xadj(i) = 0
!  end do
  xadj(0) = 0
#endif

  Do k = 1, num_blocks_z
     km = k-1
     If (km.Eq.0) km=num_blocks_z
     kp = k+1
     If (kp.Eq.num_blocks_z+1) kp=1
     Do j = 1, num_blocks_y
        jm = j-1
        If (jm.Eq.0) jm=num_blocks_y
        jp = j+1
        If (jp.Eq.num_blocks_y+1) jp=1
        Do i = 1, num_blocks_x
           im = i-1
           If (im.Eq.0) im=num_blocks_x
           ip = i+1
           If (ip.Eq.num_blocks_x+1) ip=1
           ! Nasty big if statement to check all neighbouring blocks to see if they have any atoms
#if defined(NONULLBLOCK)
           If (.True.) Then
#else
           If (&
               ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jm-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jm-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jm-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jp-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jp-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jp-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &

               ntot(i+(j-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(j-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(j-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jm-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jm-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jm-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jp-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jp-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jp-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &

               ntot(i+(j-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(j-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(j-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jm-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jm-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jm-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jp-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jp-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jp-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 &
           ) Then
#endif
              vtxdist(i,j,k) = bcnt
              bcnt = bcnt+1
              xadj(bcnt) = xadj(bcnt-1)
              vwgt(bcnt) = Nint(100*ntot(i+(j-1)*num_blocks_x+(k-1)*  &
                  (num_blocks_x*num_blocks_y))*avgblksInv)
           Else
              vtxdist(i,j,k) = -1
           Endif
        End Do
     End Do
  End Do

  cnt=0
  bcnt=0
  Do k = 1, num_blocks_z
     km = k-1
     If (km.Eq.0) km=num_blocks_z
     kp = k+1
     If (kp.Eq.num_blocks_z+1) kp=1
     Do j = 1, num_blocks_y
        jm = j-1
        If (jm.Eq.0) jm=num_blocks_y
        jp = j+1
        If (jp.Eq.num_blocks_y+1) jp=1
        Do i = 1, num_blocks_x
           im = i-1
           If (im.Eq.0) im=num_blocks_x
           ip = i+1
           If (ip.Eq.num_blocks_x+1) ip=1

#ifdef DODIAGS
           adjncy(cnt+00) = vtxdist(im,jm,km)
           adjncy(cnt+01) = vtxdist(i, jm,km)
           adjncy(cnt+02) = vtxdist(ip,jm,km)
           adjncy(cnt+03) = vtxdist(im,j, km)
           adjncy(cnt+04) = vtxdist(i, j, km)
           adjncy(cnt+05) = vtxdist(ip,j, km)
           adjncy(cnt+06) = vtxdist(im,jp,km)
           adjncy(cnt+07) = vtxdist(i, jp,km)
           adjncy(cnt+08) = vtxdist(ip,jp,km)
           adjncy(cnt+09) = vtxdist(im,jm,k)
           adjncy(cnt+10) = vtxdist(i, jm,k)
           adjncy(cnt+11) = vtxdist(ip,jm,k)
           adjncy(cnt+12) = vtxdist(im,j, k)
!           adjncy(cnt+13) = vtxdist(i, j, k)
           adjncy(cnt+13) = vtxdist(ip,j, k)
           adjncy(cnt+14) = vtxdist(im,jp,k)
           adjncy(cnt+15) = vtxdist(i, jp,k)
           adjncy(cnt+16) = vtxdist(ip,jp,k)
           adjncy(cnt+17) = vtxdist(im,jm,kp)
           adjncy(cnt+18) = vtxdist(i, jm,kp)
           adjncy(cnt+19) = vtxdist(ip,jm,kp)
           adjncy(cnt+20) = vtxdist(im,j,kp)
           adjncy(cnt+21) = vtxdist(i, j,kp)
           adjncy(cnt+22) = vtxdist(ip,j,kp)
           adjncy(cnt+23) = vtxdist(im,jp,kp)
           adjncy(cnt+24) = vtxdist(i, jp,kp)
           adjncy(cnt+25) = vtxdist(ip,jp,kp)
           cnt=cnt+26
#else
#if defined(NONULLBLOCK)
           If (.True.) Then
#else
           If (&
               ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jm-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jm-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jm-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jp-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jp-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jp-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &

               ntot(i+(j-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(j-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(j-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jm-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jm-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jm-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jp-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jp-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jp-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &

               ntot(i+(j-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(j-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(j-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jm-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jm-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jm-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jp-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(jp-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(jp-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0 &
           ) Then
#endif
              bcnt = bcnt+1
              If (bcnt.Gt.0) xadj(bcnt) = xadj(bcnt-1)
           Endif
#if defined(NONULLBLOCK)
           If (.True.) Then
#else
           If (ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(j-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y)).Gt.0 ) Then
#endif
              adjncy(cnt) = vtxdist(i,j,km)
              xadj(bcnt) = xadj(bcnt)+1
              adjwgt(cnt+1) = ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y))+ &
                                 ntot(i+(j-1)*num_blocks_x+(km-1)*(num_blocks_x*num_blocks_y))
              cnt=cnt+1
           Endif
#if defined(NONULLBLOCK)
           If (.True.) Then
#else
           If (ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(im+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0) Then
#endif
              adjncy(cnt) = vtxdist(im,j,k)
              xadj(bcnt) = xadj(bcnt)+1
              adjwgt(cnt+1) = ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y))+ &
                                 ntot(im+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y))
              cnt=cnt+1
           Endif
#if defined(NONULLBLOCK)
           If (.True.) Then
#else
           If (ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(ip+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0) Then
#endif
              adjncy(cnt) = vtxdist(ip,j,k)
              xadj(bcnt) = xadj(bcnt)+1
              adjwgt(cnt+1) = ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y))+ &
                                 ntot(ip+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y))
              cnt=cnt+1
           Endif
#if defined(NONULLBLOCK)
           If (.True.) Then
#else
           If (ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jm-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0) Then
#endif
              adjncy(cnt) = vtxdist(i,jm,k)
              xadj(bcnt) = xadj(bcnt)+1
              adjwgt(cnt+1) = ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y))+ &
                                 ntot(i+(jm-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y))
              cnt=cnt+1
           Endif
#if defined(NONULLBLOCK)
           If (.True.) Then
#else
           If (ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(jp-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0) Then
#endif
              adjncy(cnt) = vtxdist(i,jp,k)
              xadj(bcnt) = xadj(bcnt)+1
              adjwgt(cnt+1) = ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y))+ &
                                 ntot(i+(jp-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y))
              cnt=cnt+1
           Endif
#if defined(NONULLBLOCK)
           If (.True.) Then
#else
           If (ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y)).Gt.0 .Or. &
               ntot(i+(j-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y)).Gt.0) Then
#endif
              adjncy(cnt) = vtxdist(i,j,kp)
              xadj(bcnt) = xadj(bcnt)+1
              adjwgt(cnt+1) = ntot(i+(j-1)*num_blocks_x+(k-1)*(num_blocks_x*num_blocks_y))+ &
                                 ntot(i+(j-1)*num_blocks_x+(kp-1)*(num_blocks_x*num_blocks_y))
              cnt=cnt+1
           Endif
!           adjncy(cnt+00) = vtxdist(i,j,km)
!           adjncy(cnt+01) = vtxdist(im,j,k)
!           adjncy(cnt+02) = vtxdist(ip,j,k)
!           adjncy(cnt+03) = vtxdist(i,jm,k)
!           adjncy(cnt+04) = vtxdist(i,jp,k)
!           adjncy(cnt+05) = vtxdist(i,j, kp)
!           cnt=cnt+6
#endif
        End Do
     End Do
  End Do

  ! Now reassign number of vertices based on number of blocks we've put in the graph
  nvtxs = bcnt

  !Allocate memory for weights
  Allocate(vsize(1:nvtxs))      ! List of cost of transmitting (related to communication cost)
  Allocate(tpwgts(1:nparts*ncon))      ! List of target weights - must add to 1.0
  Allocate(ubvec(1:ncon))      ! List of target tolerances for imbalance

  ! Set default METIS options
  Allocate(options(0:39))          ! ParMETIS options list
  Call METIS_SetDefaultOptions(options)

  ! Assign weights
  tgtwgt=1.0/noprocs
  vsize(:)=1
  split(:)=0
  Do i=1, nblocks
#ifdef DODIAGS
     vwgt(i)=Nint(100.0_PRC*ntot(i)*avgblksInv)
#endif

! Do a quick count of how many blocks have some many particles in 
     If (ntot(i).Eq.0)  Then
        split(1) = split(1)+1
     Else If (ntot(i).Le.5) Then
        split(2) = split(2)+1
     Else
        split(3) = split(3) + 1
     Endif
  End Do
!  print*, 'Split [0, 1->5, 6+]', split
  tpwgts(:)=tgtwgt
  ubvec(:)=1.001

  ! Allocate memory for output
  Allocate(part(1:nblocks))      ! List of 

  ! Now call METIS
  If (noprocs.Gt.1) Then
     Print *, '-----------------------------------------------------------------------'
     Print *, 'Calling METIS'
     options(2) = 1 ! objtype - edgecut minimization or comm volume minimization
     options(3) = 1 ! ctype
     options(4) = 0 ! iptype
     options(5) = 2 ! rtype
     options(8) = 1  ! ncuts - i.e. the number of partitions the best is chosen from
     options(7) = 30 ! niter
     options(12) = 1 ! contig - keeps partitions whole
     Call METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt,nparts,tpwgts,ubvec,options,objval,part)
     Print *, '-----------------------------------------------------------------------'

     ! Finally put part back onto full block structure including any null blocks (already set as -1 in vtxdist)
     bcnt = 1
     Do k = 1, num_blocks_z
        Do j = 1, num_blocks_y
           Do i = 1, num_blocks_x
              If (vtxdist(i,j,k).Gt.-1) Then
                 vtxdist(i,j,k) = part(bcnt)
                 bcnt = bcnt + 1
              Endif
           End Do
        End Do
     End Do

!     call exaMD_print_partition(num_blocks_x, num_blocks_y, num_blocks_z, vtxdist, noprocs)
!     call exaMD_find_partition_extent(num_blocks_x,num_blocks_y,num_blocks_z,part,noprocs)

  Else
     Print*, 'Serial case'
     vtxdist(:,:,:)=0
  Endif

  Deallocate(xadj,adjncy,vwgt,vsize,adjwgt,tpwgts,ubvec,options,part)

End Subroutine exaMD_partition_domain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_print_partition()
!!!!!
!!!!!  * Writes a BOV file that may be opened in VisIt, etc, to display
!!!!!    partitioning information.  Location currently hardwired.
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_print_partition(num_blocks_x, num_blocks_y, num_blocks_z, part, noprocs)

  Implicit None

  Integer,Intent(in) :: num_blocks_x,num_blocks_y,num_blocks_z, noprocs
  Integer,Intent(in) :: part(:,:,:)

  Integer :: i, j, k, output
  Print*,'part'
  output=63

  Open (unit=output, file='MUPHY-partition.bov', status='replace',form='unformatted')
  Write(output) part(1:num_blocks_x,1:num_blocks_y,1:num_blocks_z)
  Close (output)
        
End Subroutine exaMD_print_partition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_find_partition_extent()
!!!!!
!!!!!  * Loops over the mesh to find the extent covered on each process.
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_find_partition_extent(num_blocks_x, num_blocks_y, num_blocks_z, part, noprocs)

  Implicit None

  Integer,Intent(in) :: num_blocks_x, num_blocks_y, num_blocks_z, noprocs
  Integer,Intent(in),Allocatable :: part(:)

  Integer :: i, j, k, kk, jj, ii, b, p, start
  Integer,Allocatable :: lastzero(:,:), lastnonzero(:,:), found(:)
  
  Allocate(lastzero(0:noprocs-1,3),lastnonzero(0:noprocs-1,3),found(0:noprocs-1))
  lastzero(:,:)=0
  lastnonzero(:,:)=0

  ! z planes
  Do k=1, num_blocks_z
     found(:)=0
     Do j=1, num_blocks_y
        Do i=1, num_blocks_x
           b = i + (j-1)*num_blocks_x + (k-1)*(num_blocks_x*num_blocks_y)
           found(part(b))=1
        End Do
     End Do
     Do p = 0, noprocs-1 
        If (found(p).Eq.0) lastzero(p,3) = k
     End Do
  End Do

  Do p = 0, noprocs-1
     start = lastzero(p,3)
     If (start.Eq.0) Then
        lastnonzero(p,3) = num_blocks_z
        Cycle
     Endif
     Do kk=start+1, start+num_blocks_z
        k = Mod(kk-1,num_blocks_z)+1
        found(:)=0
        Do j=1, num_blocks_y
           Do i=1, num_blocks_x
              b = i + (j-1)*num_blocks_x + (k-1)*(num_blocks_x*num_blocks_y)
              found(part(b))=1
           End Do
        End Do
        If (found(p).Eq.1) lastnonzero(p,3) = k
     End Do
  End Do

  ! y planes
  Do j=1, num_blocks_y
     found(:)=0
     Do k=1, num_blocks_z
        Do i=1, num_blocks_x
           b = i + (j-1)*num_blocks_x + (k-1)*(num_blocks_x*num_blocks_y)
           found(part(b))=1
        End Do
     End Do
     Do p = 0, noprocs-1 
        If (found(p).Eq.0) lastzero(p,2) = j
     End Do
  End Do

  Do p = 0, noprocs-1
     start = lastzero(p,2)
     If (start.Eq.0) Then
        lastnonzero(p,2) = num_blocks_y
        Cycle
     Endif
     Do jj=start+1, start+num_blocks_y
        j = Mod(jj-1,num_blocks_y)+1
        found(:)=0
        Do k=1, num_blocks_z
           Do i=1, num_blocks_x
              b = i + (j-1)*num_blocks_x + (k-1)*(num_blocks_x*num_blocks_y)
              found(part(b))=1
           End Do
        End Do
        If (found(p).Eq.1) lastnonzero(p,2) = j
     End Do
  End Do

  ! x planes
  Do i=1, num_blocks_x
     found(:)=0
     Do k=1, num_blocks_z
        Do j=1, num_blocks_y
           b = i + (j-1)*num_blocks_x + (k-1)*(num_blocks_x*num_blocks_y)
           found(part(b))=1
        End Do
     End Do
     Do p = 0, noprocs-1 
        If (found(p).Eq.0) lastzero(p,1) = i
     End Do
  End Do

  Do p = 0, noprocs-1
     start = lastzero(p,1)
     If (start.Eq.0) Then
        lastnonzero(p,1) = num_blocks_x
        Cycle
     Endif
     Do ii=start+1, start+num_blocks_x
        i = Mod(ii-1,num_blocks_x)+1
        found(:)=0
        Do k=1, num_blocks_z
           Do j=1, num_blocks_y
              b = i + (j-1)*num_blocks_x + (k-1)*(num_blocks_x*num_blocks_y)
              found(part(b))=1
           End Do
        End Do
        If (found(p).Eq.1) lastnonzero(p,1) = i
     End Do
  End Do


  Write(6,'(A,I3,A,I3,A,I3,A)') 'Ranges: out of (',num_blocks_x,',',num_blocks_y,',',num_blocks_z,')'
  Do p = 0, noprocs-1
     Write(6,'(I4,2X)',advance='no') p
  End Do
  Write(6,*) ' ' 
  Do p = 0, noprocs-1
     Write(6,'(I2,A,I2,X)',advance='no') Mod(lastzero(p,1)+1-1,num_blocks_x)+1,'-',lastnonzero(p,1)
  End Do
  Write(6,*) ' ' 
  Do p = 0, noprocs-1
     Write(6,'(I2,A,I2,X)',advance='no') Mod(lastzero(p,2)+1-1,num_blocks_y)+1,'-',lastnonzero(p,2)
  End Do
  Write(6,*) ' ' 
  Do p = 0, noprocs-1
     Write(6,'(I2,A,I2,X)',advance='no') Mod(lastzero(p,3)+1-1,num_blocks_z)+1,'-',lastnonzero(p,3)
  End Do
  Write(6,*) ' ' 
  Do p = 0, noprocs-1
     Write(6,'(F5.2,X)',advance='no')&
     Real(Mod(lastnonzero(p,1)+num_blocks_x-(Mod(lastzero(p,1)+1-1,num_blocks_x)+1),num_blocks_x)+1,KIND=PRC)  *&
     (Mod(lastnonzero(p,2)+num_blocks_y-(Mod(lastzero(p,2)+1-1,num_blocks_y)+1),num_blocks_y)+1)*&
     (Mod(lastnonzero(p,3)+num_blocks_z-(Mod(lastzero(p,3)+1-1,num_blocks_z)+1),num_blocks_z)+1)/&
     Real(0.01*num_blocks_x*num_blocks_y*num_blocks_z,KIND=PRC)
  End Do
  Write(6,*) ' ' 


  Deallocate(lastzero,lastnonzero,found)

End Subroutine exaMD_find_partition_extent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!! 
!!!!!  exaMD_setup_partitioned_grid() 
!!!!! 
!!!!!  * Run by all processes this routine sets up all information on the 
!!!!!    partitioned grid related to parallel partitioning.  Nothing to do with
!!!!!    block data is set here.
!!!!! 
!!!!!  * Lists are formed of interior, boundary and exterior blocks.
!!!!! 
!!!!!  * Lists of neighbouring processes and which blocks neighbour 
!!!!!    processes also formed. 
!!!!! 
!!!!!  Returns: 
!!!!!  * (none) 
!!!!!  
!!!!!  Global structure data set 
!!!!!  * Most block structure and lists set up on each process 
!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!

Subroutine exaMD_setup_partitioned_grid(part, nblocks, nblocks_xyz, pid)

  Use exaMD_mod

  Implicit None

  Integer, Intent(IN) :: nblocks, nblocks_xyz(3), pid
  Integer, Allocatable, Intent(IN) :: part(:,:,:)
  
  Integer :: i, j, k, ii, jj, kk, iii, jjj, kkk, b, cnt, owned, remoteproc, found, r
  Integer :: offprocneigh, nullproccount, num_blocks_x, num_blocks_y, num_blocks_z
  Integer, Allocatable :: partnerprocs(:)
  Type(pframe_t), Pointer :: block

  num_blocks_x = nblocks_xyz(1)
  num_blocks_y = nblocks_xyz(2)
  num_blocks_z = nblocks_xyz(3)

  ! Allocate memory for array detailing which processors we need to communicate with
  
  ! Put part data into the grid structure and count how many blocks are ours
!  cnt = 1
  Do k = 1, num_blocks_z
     Do j = 1, num_blocks_y
        Do i = 1, num_blocks_x
           ! Use a POINTER to speed up accesses
           block => grid(i,j,k)
           block%host = part(i,j,k)
           block%nsendtoremote = 0
        End Do
     End Do
  End Do

  cnt = 1
  nblocks_loc = 0
  n_interior = 0
  n_bdry = 0
  n_exterior = 0
  Do k = 1, num_blocks_z
     Do j = 1, num_blocks_y
        Do i = 1, num_blocks_x
           ! Use a POINTER to speed up accesses
           block => grid(i,j,k)
           If (part(i,j,k).Eq.pid) nblocks_loc = nblocks_loc+1

           offprocneigh = 0
           nullproccount = 0
           Do kk=-1,1
              Do jj=-1,1
                 Do ii=-1,1
                    iii=Mod(i+ii+num_blocks_x-1,num_blocks_x)+1
                    jjj=Mod(j+jj+num_blocks_y-1,num_blocks_y)+1
                    kkk=Mod(k+kk+num_blocks_z-1,num_blocks_z)+1
                    remoteproc = part(iii,jjj,kkk)

                    ! Mark that the host of this neighbour is on a different process
                    If (remoteproc.Eq.-1) Then
                       nullproccount = nullproccount + 1
                    Else If (remoteproc.Ne.pid) Then
                       offprocneigh = offprocneigh + 1
                    Endif                    
                  End Do
              End Do
           End Do

           If (offprocneigh+nullproccount.Eq.27) Then
              block%partype = 0
           Else If (offprocneigh.Eq.0.And.block%host.Eq.pid) Then
              block%partype = 1
              n_interior = n_interior + 1
           Else If (block%host.Eq.pid) Then
              block%partype = 2
              n_bdry = n_bdry + 1
           Else 
              block%partype = 3
              n_exterior = n_exterior + 1
           Endif
           cnt = cnt + 1
        End Do
     End Do
  End Do

  Allocate(locblocks(3,nblocks_loc))
  Allocate(interior_blks(3,n_interior))
  Allocate(bdry_blks(3,n_bdry))
  Allocate(exterior_blks(3,n_exterior))

  ! Now go through and add the information to the nodes we own
  nblocks_loc = 0
  n_interior = 0
  n_bdry = 0
  n_exterior = 0

  b = 1
  Do k = 1, num_blocks_z
     Do j = 1, num_blocks_y
        Do i = 1, num_blocks_x
           ! Use a POINTER to speed up accesses
           block => grid(i,j,k)
           If (block%partype.Gt.0) Then        ! i.e. needs something to do with this block
              ! Add this block to the list of blocks on this processor
              If (block%partype.Eq.1) Then     ! Interior block
                 nblocks_loc = nblocks_loc+1
                 locblocks(1,nblocks_loc) = i
                 locblocks(2,nblocks_loc) = j
                 locblocks(3,nblocks_loc) = k

                 n_interior = n_interior+1
                 interior_blks(1,n_interior) = i
                 interior_blks(2,n_interior) = j
                 interior_blks(3,n_interior) = k

                 ! Now fill list of adjacent neighbours
                 Allocate(block%neighrefs(1:3,1:27))
                 cnt=1
                 offprocneigh = 0
                 Do kk=-1,1
                    Do jj=-1,1
                       Do ii=-1,1
                          iii=Mod(i+ii+num_blocks_x-1,num_blocks_x)+1
                          jjj=Mod(j+jj+num_blocks_y-1,num_blocks_y)+1
                          kkk=Mod(k+kk+num_blocks_z-1,num_blocks_z)+1
                          block%neighrefs(1,cnt) = iii
                          block%neighrefs(2,cnt) = jjj
                          block%neighrefs(3,cnt) = kkk
                          cnt = cnt+1
                       End Do  
                    End Do  
                 End Do  
              Else If (block%partype.Eq.3) Then ! Exterior block
                 n_exterior = n_exterior+1
                 exterior_blks(1,n_exterior) = i
                 exterior_blks(2,n_exterior) = j
                 exterior_blks(3,n_exterior) = k

              Else If (block%partype.Eq.2) Then ! Boundary block
                 nblocks_loc = nblocks_loc+1
                 locblocks(1,nblocks_loc) = i
                 locblocks(2,nblocks_loc) = j
                 locblocks(3,nblocks_loc) = k

                 n_bdry = n_bdry+1
                 block%flag = .False.
                 bdry_blks(1,n_bdry) = i
                 bdry_blks(2,n_bdry) = j
                 bdry_blks(3,n_bdry) = k
              
                 Allocate(block%neighhosts(27))
                 Allocate(block%neighrefs(1:3,1:27))

                 ! Now fill list of adjacent neighbours
                 cnt=1
                 offprocneigh = 0
                 Do kk=-1,1
                    Do jj=-1,1
                       Do ii=-1,1
                          iii=Mod(i+ii+num_blocks_x-1,num_blocks_x)+1
                          jjj=Mod(j+jj+num_blocks_y-1,num_blocks_y)+1
                          kkk=Mod(k+kk+num_blocks_z-1,num_blocks_z)+1
                          block%neighrefs(1,cnt) = iii
                          block%neighrefs(2,cnt) = jjj
                          block%neighrefs(3,cnt) = kkk
                          block%neighhosts(cnt) = grid(iii,jjj,kkk)%host

                          ! Mark that the host of this neighbour is on a different processor
                          remoteproc = grid(iii,jjj,kkk)%host
                          If (remoteproc.Ne.pid.And.remoteproc.Ge.0) Then
                             grid(iii,jjj,kkk)%getfromremote = 1
                             If (block%nsendtoremote.Gt.0) Then
                                found = 0
                                Do r=1, block%nsendtoremote
                                   If (block%sendtoremote(r).Eq.remoteproc) Then
                                      found = 1
                                      Exit
                                   End If
                                End Do
                                If (found.Eq.0) Then
                                   block%nsendtoremote = block%nsendtoremote+1
                                   block%sendtoremote(block%nsendtoremote) = remoteproc
                                Endif
                             Else
                                Allocate(block%sendtoremote(26))
                                block%sendtoremote(2:26) = -1
                                block%nsendtoremote = 1
                                block%sendtoremote(block%nsendtoremote) = remoteproc
                             Endif
                          Endif
                          cnt = cnt+1
                       End Do  
                    End Do  
                 End Do  
              Endif
           End If
           b = b + 1
        End Do  
     End Do  
  End Do  

End Subroutine exaMD_setup_partitioned_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_partition_work()
!!!!!
!!!!!  * Calculates data sizes for transmission to neighbouring processes
!!!!!    for both positions and forces.
!!!!!  
!!!!!  * Blocks are also fully set up for use, initialising posn and force 
!!!!!    arrays.
!!!!! 
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * More block data
!!!!!  * send and receive buffers
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_partition_work(particle_positions, ntot, first_particle, next_particle, nblocks, nblocks_xyz, pid, noprocs)

  Use exaMD_mod
  Implicit None

  Integer, Intent(IN) :: nblocks, nblocks_xyz(3), pid, noprocs
  Integer, Allocatable, Intent(IN) :: ntot(:), first_particle(:), next_particle(:)
  Real(KIND=PRC), Allocatable, Intent(IN) :: particle_positions(:,:)

  Integer :: i, j, k, ii, jj, kk, b, cnt, tgt, tag, p, pp, loop, particle_index
  Integer :: num_blocks_x, num_blocks_y, num_blocks_z
  Integer, Allocatable :: datasize(:)
  Integer, Allocatable :: req(:)  ! MPI Request handles
  Type(pframe_t), Pointer :: block

  Call MPI_Bcast(ntot, nblocks, MPI_INTEGER, 0, COMM_EXAMD, ierr)

  num_blocks_x = nblocks_xyz(1)
  num_blocks_y = nblocks_xyz(2)
  num_blocks_z = nblocks_xyz(3)

  ! Root process to distribute work
  If (pid.Eq.0) Then
     Allocate(req(noprocs-1))

     Allocate(datasize(0:noprocs-1))
     datasize(:) = 0
     Do k = 1, num_blocks_z
        Do j = 1, num_blocks_y
           Do i = 1, num_blocks_x
              block => grid(i,j,k)
              If (block%host>-1) Then 
                 tgt = block%host
                 b = i + (j-1)*num_blocks_x + (k-1)*(num_blocks_x*num_blocks_y)
                 datasize(tgt) = datasize(tgt)+ntot(b)  
              Endif
           End Do
        End Do
     End Do
     Allocate(sendbuffer(3,0:Maxval(datasize),0:noprocs-1))

     datasize(:) = 0
     cnt = 0
     Do k = 1, num_blocks_z
        Do j = 1, num_blocks_y
           Do i = 1, num_blocks_x
              block => grid(i,j,k)
              tgt = block%host
              b = i + (j-1)*num_blocks_x + (k-1)*(num_blocks_x*num_blocks_y)
              block%natms = ntot(b)  
              If (tgt.Gt.0) Then
                 If (ntot(b).Gt.0) Then  ! i.e. don't waste time sending nothing
                    ! Build up a per processor set of data to send
                    particle_index = first_particle(b)
                    Do loop = 1, ntot(b)
                       sendbuffer(:,datasize(tgt),tgt) = particle_positions(:,particle_index)
                       datasize(tgt) = datasize(tgt) + 1
                       particle_index = next_particle(particle_index)
                    End Do

                    If (block%partype.Eq.3) Then ! i.e. we'll need this block back
                       Allocate(block%ltype(ntot(b)))
                    Endif
                 Endif
              Else If (tgt.Eq.0) Then  ! exclude the null case
                 Allocate(block%posn(1:3,ntot(b)))
                 Allocate(block%forces(1:3,1:ntot(b)))
                 Allocate(block%ltype(ntot(b)))

                 particle_index = first_particle(b)
                 Do loop = 1, ntot(b)
                    block%posn(:,loop) = particle_positions(:,particle_index)
                    particle_index = next_particle(particle_index)
                 End Do

                 block%forces(:,:) = zero
              Endif
           End Do
        End Do
     End Do

     ! Note using cnt rather than nblocks*4 as the count to ignore local blocks

     If (noprocs.Gt.1) Then
        Do p = 1, noprocs-1
           Call MPI_Isend(sendbuffer(1,0,p), 3*datasize(p), MYFLOAT, p, 100, COMM_EXAMD, req(p), ierr)
        End Do        
        Call MPI_Waitall(noprocs-1, req, MPI_STATUSES_IGNORE, ierr)
     End If

     Deallocate(req)
     Deallocate(sendbuffer)

  Else ! i.e. all processes apart from 0

     ! Other processors only to work through the blocks they're expecting in order

     Allocate(datasize(0:0))
     datasize(:) = 0
     Do b = 1, nblocks_loc
        ii = locblocks(1,b)
        jj = locblocks(2,b)
        kk = locblocks(3,b)
        block => grid(ii,jj,kk)
        tag = ii + (jj-1)*num_blocks_x + (kk-1)*(num_blocks_x*num_blocks_y)
        datasize(0) = datasize(0)+ntot(tag)  
     End Do

     Allocate(recvbuffer(3,0:datasize(0)-1,0:0))
     Call MPI_Recv(recvbuffer, 3*datasize(0), MYFLOAT, 0, 100, COMM_EXAMD, MPI_STATUS_IGNORE, ierr)

     datasize(0) = 0

     cnt = 0
     Do b = 1, nblocks_loc
        ii = locblocks(1,b)
        jj = locblocks(2,b)
        kk = locblocks(3,b)
        block => grid(ii,jj,kk)
        tag = ii + (jj-1)*num_blocks_x + (kk-1)*(num_blocks_x*num_blocks_y)
        block%natms = ntot(tag)
        Allocate(block%posn(1:3,ntot(tag)))
        Allocate(block%forces(1:3,1:ntot(tag)))
        Allocate(block%ltype(1:ntot(tag)))

        block%forces(:,:) = zero

        If (ntot(tag).Gt.0) Then
           block%posn(1:3,1:ntot(tag)) = recvbuffer(1:3,datasize(0):datasize(0)+ntot(tag)-1,0)
           datasize(0) = datasize(0)+ntot(tag)
        Endif
     End Do
     
     ! Whilst waiting for messages, let's fill the number of atoms in on the neighbouring (exterior) blocks across the domain.
     Do b = 1, n_exterior
        ii = exterior_blks(1,b)
        jj = exterior_blks(2,b)
        kk = exterior_blks(3,b)
        block => grid(ii,jj,kk)
        tgt = ii + (jj-1)*num_blocks_x + (kk-1)*(num_blocks_x*num_blocks_y)
        block%natms = ntot(tgt)  
        Allocate(block%ltype(block%natms))
     End Do

     Deallocate(recvbuffer)
  End If

  ! All processes execute the rest of the subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
  ! Set up arrays for the communication
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

  Allocate(datasize_s(0:noprocs-1),datasize_r(0:noprocs-1))
  datasize_s(:)=0
  datasize_r(:)=0
  Allocate(datasize_fs(0:noprocs-1),datasize_fr(0:noprocs-1))
  datasize_fr(:)=0
  
  ! Calculate size of messages to receive
  Do b = 1, n_exterior
     ii = exterior_blks(1,b)
     jj = exterior_blks(2,b)
     kk = exterior_blks(3,b)
     block => grid(ii,jj,kk)
     p = block%host
     If (block%host>-1) Then
        datasize_r(p) = datasize_r(p) + Max(block%natms,1)  

        ! Do some initialisation whilst we're here
        Allocate(block%posn(3,block%natms))
        Allocate(block%forces(1:3,1:block%natms))

        block%forces(:,:) = zero
     Endif
  End Do

  ! Count the number of processors we are sending to, and make a reference array.
  ! This will be the same set we want data from.
  n_extneigh = 0
  Do p = 0, noprocs-1
     If (p.Eq.pid) Cycle
     If (datasize_r(p).Gt.0) Then
        n_extneigh = n_extneigh + 1
     Endif
  End Do 
  Allocate(ext_neigh(n_extneigh))
  ! Add arrays for storing target and running totals of boundary blocks' neighbours being sent to each processor 
  Allocate(n_frcsendtgt(n_extneigh))
  Allocate(n_frcsendcnt(n_extneigh))
  n_frcsendtgt(:) = 0 ! Initialise
  n_frcsendcnt(:) = 0 ! Initialise
  n_extneigh = 0
  Do p = 0, noprocs-1
     If (p.Eq.pid) Cycle
     If (datasize_r(p).Gt.0) Then
        n_extneigh = n_extneigh + 1
        ext_neigh(n_extneigh) = p
     Endif
  End Do 

  ! Calculate size of messages to send
  Do b = 1, n_bdry
     ii = bdry_blks(1,b)
     jj = bdry_blks(2,b)
     kk = bdry_blks(3,b)
     block => grid(ii,jj,kk)
     p = block%host
     If (block%host>-1) Then
        datasize_s(p) = datasize_s(p) + block%natms  
        Do i = 1, block%nsendtoremote
           p = block%sendtoremote(i)
           If (p.Gt.-1) Then
              datasize_s(p) = datasize_s(p) + block%natms  

              If (block%sendtoremote(i).Ne.pid) Then
                 Do pp = 1, n_extneigh
                    If (p.Eq.ext_neigh(pp)) Exit
                 End Do
                 If (pp.Eq.n_extneigh+1) Print*,'ERROR in partition work reverse lookup'
                 n_frcsendtgt(pp) = n_frcsendtgt(pp) + 1
              Endif
           Endif
        End Do
     End If
  End Do

  datasize_fs(:) = datasize_s(:)

  ! Now make lists of blocks on the boundary related to each neighbouring processor
  Call exaMD_partition_neighbours()

End Subroutine exaMD_partition_work

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  exaMD_partition_neighbours()
!!!!!
!!!!!  * This routine forms lists of boundary and exterior blocks related to
!!!!!    which process needs them.  A boundary block, for example, may be in
!!!!!    multiple lists.  This is used in the full overlap case.
!!!!!
!!!!!  Returns:
!!!!!  * (none)
!!!!!  
!!!!!  Global structure data set
!!!!!  * pbdryblk_lst - list of boundary blocks needed by each neighbour
!!!!!  * n_pbdryblk_lst - number of boundary blocks needed by each neighbour
!!!!!  * pextblk_lst - list of exterior blocks needed by each neighbour
!!!!!  * n_pextblk_lst - number of exterior blocks needed by each neighbour
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine exaMD_partition_neighbours()

  Use exaMD_mod

  Implicit None

  Integer :: i, j, k, ii, jj, kk, b, cnt, tgt, tag, p, pp
  Integer, Allocatable :: counts(:)
  Type(pframe_t), Pointer :: block

  If (n_bdry.Eq.0) Return

  ! Allocate memory for counter arrays
  Allocate(counts(n_extneigh),n_pbdryblk_lst(n_extneigh+1),n_pextblk_lst(n_extneigh+1))
  counts(:) = 0

  ! Boundary lists first

  ! Count how many boundary blocks have neighbours on each process
  Do b = 1, n_bdry
     ii = bdry_blks(1,b)
     jj = bdry_blks(2,b)
     kk = bdry_blks(3,b)
     block => grid(ii,jj,kk)

     Do j = 1, block%nsendtoremote
        p = block%sendtoremote(j)
        Do pp = 1, n_extneigh
           If (p.Eq.ext_neigh(pp)) Exit
        End Do
        counts(pp) = counts(pp) + 1
     End Do
  End Do

  ! Fill the n_pbdryblk_lst
  n_pbdryblk_lst(1) = 1
  Do i = 2, n_extneigh+1
     n_pbdryblk_lst(i) = n_pbdryblk_lst(i-1) + counts(i-1)
  End Do

  ! Allocate memory for the lists
  Allocate(pbdryblk_lst(Sum(counts(:))))

  ! Reset counter
  counts(:) = 0

  ! Count how many boundary blocks have neighbours on each process
  Do b = 1, n_bdry
     ii = bdry_blks(1,b)
     jj = bdry_blks(2,b)
     kk = bdry_blks(3,b)
     block => grid(ii,jj,kk)

     Do j = 1, block%nsendtoremote
        p = block%sendtoremote(j)
        Do pp = 1, n_extneigh
           If (p.Eq.ext_neigh(pp)) Exit
        End Do
        pbdryblk_lst(n_pbdryblk_lst(pp)+counts(pp)) = b
        counts(pp) = counts(pp) + 1
     End Do
  End Do

  ! Exterior lists

  ! Reset counter
  counts(:) = 0

  ! Count how many exterior blocks are on each process
  Do b = 1, n_exterior
     ii = exterior_blks(1,b)
     jj = exterior_blks(2,b)
     kk = exterior_blks(3,b)
     block => grid(ii,jj,kk)

     p = block%host
     If (p.Eq.-1) Cycle
     Do pp = 1, n_extneigh
        If (p.Eq.ext_neigh(pp)) Exit
     End Do
     counts(pp) = counts(pp) + 1
  End Do

  ! Fill the n_pextblk_lst
  n_pextblk_lst(1) = 1
  Do i = 2, n_extneigh+1
     n_pextblk_lst(i) = n_pextblk_lst(i-1) + counts(i-1)
  End Do

  ! Allocate memory for the lists
  Allocate(pextblk_lst(Sum(counts(:))))

  ! Reset counter
  counts(:) = 0

  ! Assign exterior blocks with neighbours on each process to lists
  Do b = 1, n_exterior
     ii = exterior_blks(1,b)
     jj = exterior_blks(2,b)
     kk = exterior_blks(3,b)
     block => grid(ii,jj,kk)

     p = block%host
     If (p.Eq.-1) Cycle
     Do pp = 1, n_extneigh
        If (p.Eq.ext_neigh(pp)) Exit
     End Do
     pextblk_lst(n_pextblk_lst(pp)+counts(pp)) = b
     counts(pp) = counts(pp) + 1
  End Do

  Deallocate(counts)

End Subroutine exaMD_partition_neighbours

End Module exaMD_partitioning

