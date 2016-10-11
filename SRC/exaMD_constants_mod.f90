! Store numerical constants

Module exaMD_constants_mod

  Implicit None

#if PRC==4
  Real(kind=PRC), Parameter :: zero  = 0.e0, &
                               one   = 1.e0, &
                               two   = 2.e0, &
                               three = 3.e0, &
                               four  = 4.e0, &
                               eight = 8.e0, &
                               twentyfour = 24.e0, &
                               half = 0.5e0, &
                               third = 0.3333333e0, &
                               fourth = 0.25e0, &
                               eighti = 0.125e0, &
                               sqrt_twelve = 3.464102e0, &
                               pi = 3.141593e0
#else
  Real(kind=PRC), Parameter :: zero  = 0.d0, &
                               one   = 1.d0, &
                               two   = 2.d0, &
                               three = 3.d0, &
                               four  = 4.d0, &
                               eight = 8.d0, &
                               twentyfour = 24.d0, &
                               half = 0.5d0, &
                               third = 0.333333333333333d0, &
                               fourth = 0.25d0, &
                               eighti = 0.125d0, &
                               sqrt_twelve = 3.46410161513775d0, &
                               pi = 3.1415926535897932d0
#endif

End Module exaMD_constants_mod
