
       parameter (pi=3.141592653589793d0)
       parameter (hbarc=197.3269680000d0)        ! hbar*c [MeV*fm]
       parameter (mp=938.2723100000000d0)        ! proton mass
       parameter (mn=939.5656300000000d0)        ! neutron mass
       parameter (charge2=1.4398840000d0)        ! electric charge to the power of 2

       INTEGER :: AZ, AN        ! proton and neutron number
       INTEGER :: noscmax       ! the maximal oscillator shell (defines the configuration space: #shells = noscmax + 1) for HFB
       INTEGER :: nosc_TDA       ! the maximal oscillator shell (defines the configuration space: #shells = noscmax + 1) for TDA
       INTEGER :: min_p,max_p,min_n,max_n

       DOUBLE PRECISION :: precis  ! parameter to stop the iteration step
       DOUBLE PRECISION :: hbarom  ! hbar*omega [MeV]

       common /spacetda/ min_p,max_p,min_n,max_n
       common /hb/ hbarom
       common /pre/ precis


!       b = sqrt(m*omega/hbar) = sqrt(m*c^2*hbar*omega)/(hbar*c) = sqrt(938.5*41*A^(-1/3) MeV^2)/(197.327 MeV*fm)
