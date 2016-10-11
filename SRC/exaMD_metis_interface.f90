!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   exaMD_metis_interface.f90
!!!!!
!!!!!   This file provides interface functions for calling the necessary METIS 
!!!!!   functions from Fortran.  No work done, simply an interface to the 
!!!!!   external C functions.
!!!!!
!!!!!   The routines contained in this module are:
!!!!!    * METIS_SetDefaultOptions
!!!!!    * METIS_PartGraphKway
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module exaMD_METIS_interface

  Use iso_c_binding

  Implicit None

  Interface
    Subroutine METIS_SetDefaultOptions(options)
      Use, Intrinsic :: iso_c_binding

      Implicit None

      Integer (c_int), Dimension(*), Intent(in) :: options

    End Subroutine METIS_SetDefaultOptions
    

    Subroutine METIS_PartGraphKway(   &
                     nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt,nparts,tpwgts,ubvec,options,objval,part) bind(c)
      Use, Intrinsic :: iso_c_binding

      Implicit None

      Integer (c_int), Dimension(*), Intent(in) :: xadj,adjncy,vwgt,vsize,adjwgt
      Integer (c_int), Dimension(*), Intent(in) :: tpwgts,ubvec,options,part
      Integer (c_int),               Intent(in) :: nvtxs,ncon, nparts, objval

    End Subroutine METIS_PartGraphKway
  

!    Subroutine ParMETIS_V3_PartKway( &
!                      vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag,   &
!                      numflag, ncon, nparts, tpwgts, ubvec, options,  &
!                      edgecut, part, comm) bind(c)                    
!      use, intrinsic        :: iso_c_binding
!
!      implicit none
!
!      integer (c_int), dimension(*),    intent(in)  :: vtxdist, xadj, adjncy, vwgt, adjwgt
!      integer (c_int),                  intent(in)  :: wgtflag, numflag, ncon, nparts
!      integer (c_int), dimension(*),    intent(in)  :: tpwgts, ubvec, options
!      integer (c_int),                  intent(in)  :: edgecut
!      integer (c_int), dimension(*),    intent(in)  :: part
!      integer (c_int),                  intent(in)  :: comm
!
!    end subroutine ParMETIS_V3_PartKway
  End Interface
End Module exaMD_METIS_interface
