!   Used types

      TYPE level_type
c      SEQUENCE
        INTEGER           :: index        ! the number of level
        INTEGER           :: ipar         ! parity
	INTEGER           :: N            ! the quantum number N=2*n+l
	INTEGER           :: nn           ! the quantum number n
	INTEGER           :: l            ! the quantum number l
        INTEGER           :: j2           ! the quantum number 2*j
        DOUBLE PRECISION  :: spenrg       ! s.p.-energy e_i = hbar*omega*(N+3/2)
	DOUBLE PRECISION  :: ei           ! s.p. HF energy
        DOUBLE PRECISION  :: qei          ! quasi-s.p. HFB energy
        DOUBLE PRECISION  :: ui           ! Bogolyubov coefficient u_i 
        DOUBLE PRECISION  :: vi           ! Bogolyubov coefficient v_i 
      END TYPE level_type

      TYPE phonon_type
        INTEGER           :: noscmax      ! show the space in which the TDA phonon was calculated
	INTEGER           :: ipar         ! parity of the phonon
	INTEGER           :: jproj        ! j-projection of the phonon
	INTEGER           :: dimp         ! the dimension of the fields cph_p
	INTEGER           :: dimn         ! the dimension of the fields cph_n
	DOUBLE PRECISION  :: ener         ! energy of the phonon
	DOUBLE PRECISION  :: cph_p(112,112)   ! the amplitudes c_ph for protons
	DOUBLE PRECISION  :: cph_n(112,112)   ! the amplitudes c_ph for neutrons
      END TYPE phonon_type

      TYPE ph_type
        INTEGER           :: par          ! parity
        INTEGER           :: Jph          ! J
        INTEGER           :: lam          ! index
        DOUBLE PRECISION  :: el           !energy
      END TYPE ph_type

      TYPE cph_type
        INTEGER           :: par          !particle state
        INTEGER           :: hol          !hole state
        INTEGER           :: tz           ! -1 for protons, +1 for neutrons
        DOUBLE PRECISION  :: en           ! energy of phonon
        DOUBLE PRECISION  :: cph          ! Cph^lam value
      END TYPE cph_type
 
      TYPE dimen_type
       integer         :: parity
       integer         :: jph
       integer         :: dimen
      END TYPE dimen_type
     
      TYPE dens_type 
       INTEGER            :: lam
       INTEGER            :: aa    !<lam'|(aa^+ x bb)^sig|lam>
       INTEGER            :: bb
       INTEGER            :: sig
       DOUBLE PRECISION   :: val 
      END TYPE dens_type  
