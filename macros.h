#define LogAlloc(ERR,NAME,LEN,SIZE,TAG) CALL LogMemAlloc(NAME,LEN,SIZE,this_routine,TAG)
#define LogDealloc(TAG) CALL LogMemDealloc(this_routine,TAG)
#define IsNullDet(nI) (nI(1).eq.0)

! Is the specified orbital occupied or not?
#define IsOcc(ilut,orb) btest(ilut((orb-1)/32), mod(orb-1,32))
#define IsNotOcc(ilut,orb) (.not.IsOcc(ilut,orb))

! Is the specified orbital alpha or beta? Generate the appropriate pair.
#define is_beta(orb) btest(orb, 0)
#define is_alpha(orb) (.not.is_beta(orb))
#define is_one_alpha_beta(orb1,orb2) (btest(orb1,0) .xor. btest(orb2,0)) 
#define ab_pair(orb) (ieor(orb-1,1)+1)
#define get_beta(orb) (ibclr(orb-1,0)+1)
#define get_alpha(orb) (ibset(orb-1,0)+1)

! Is the specified orbital part of a doubly occupied pair?
#define IsDoub(ilut,orb) (IsOcc(ilut,orb) .and. IsOcc(ilut,ab_pair(orb)))

! Are the two orbitals specified (may be the same orbital) from the same
! spatial orbital?
#define is_in_pair(orb1,orb2) (ibclr(orb1-1,0) == ibclr(orb2-1,0))

