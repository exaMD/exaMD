!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   user_func.f90
!!!!!
!!!!!   This file contains an example particle interaction function for use with
!!!!!   the exaMD proto-app library for a molecular dynamics calculation.
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module user_func_mod

  Use exaMD_mod

  Implicit None

  ! User level data
  Integer, Allocatable :: kktype(:,:)        ! table of chemical species
  Integer, Allocatable :: ltpvdw(:)          ! table of forces id
  Real(kind=PRC), Allocatable :: prmvdw(:,:) ! parameters of vdw forces

  Public :: lennardjones_12

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  lennardjones_12()
!!!!!
!!!!!  * Example user function for computing forces between particles
!!!!!
!!!!!  Returns:
!!!!!  * f12 - array that updates the forces data structure per block
!!!!!  * e12 - array that is unused
!!!!!  
!!!!!  Global structure data set
!!!!!  * (none)
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!DEC$ ATTRIBUTES FORCEINLINE :: lennardjones_12
Subroutine lennardjones_12(dx, dy, dz, d2, f12, e12, retval)

  Implicit None

  Real(kind=PRC), Intent(in) :: dx,dy,dz,d2
  Real(kind=PRC), Intent(out) :: f12(3),e12
  Integer, Intent(out) :: retval

  Real(kind=PRC) :: rc2,eps4,eps24,sig2,d6,d12,d2i,dg,c2i,c6,c12

  retval = 0

  If (d2 > prmvdw(1,5)**2) Then  ! test cut-off distance
!     f12(1:3) = zero
!     e12 = zero
     retval = 1    ! Set to not update the forces vector
  Else

     eps24 = twentyfour*prmvdw(1,1)
     sig2 = prmvdw(1,2)**2

     d2i = one/d2
     d6 = (sig2*d2i)**3
     d12 = d6*d6

     dg = - eps24*(two*d12 - d6)*d2i

     f12(1) = dg*dx
     f12(2) = dg*dy
     f12(3) = dg*dz

     e12 = zero

  End If

End Subroutine lennardjones_12

End Module user_func_mod
