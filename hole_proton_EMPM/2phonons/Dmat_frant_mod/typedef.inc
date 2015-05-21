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

      TYPE ph_type
        INTEGER           :: par          ! parity
        INTEGER           :: Jph          ! J
        INTEGER           :: lam          ! index
        DOUBLE PRECISION  :: el           !energy
      END TYPE ph_type

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

      TYPE indbas_type
       INTEGER            :: ilam
       INTEGER            :: ir
      END TYPE indbas_type

      TYPE ph2_type
        INTEGER           :: par          ! parity
        INTEGER           :: J            ! J
        INTEGER           :: alfa         ! alfa index
        DOUBLE PRECISION  :: el           !energy
      END TYPE ph2_type

