! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
module gnd_work_type

use SystemData, only: BasisFn, BasisFNSize

! Formerly gndwork.inc.
! Probably only Alex knows what it's for...
!
! You've only to ask!
! Structure holding working data passed recursively into GenNextDet_.  All rather complicated I'm afraid.
!  AJWT 20110121

type GNDWork
    integer        NSWORK(4)
    type(BasisFN)  IMax(2)
    type(BasisFN)  ISym
    integer        nElec
    integer        niWork(1)
    ! (nEl)
    ! INTEGER        nIndJ(nEl)
    ! nIndJ is niWork(1+nEl)
end type

integer, parameter :: GNDWorkSize=6+3*BasisFNSize

end module gnd_work_type
