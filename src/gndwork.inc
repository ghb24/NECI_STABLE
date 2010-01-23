         TYPE GNDWork
            INTEGER        NSWORK(4)
            TYPE(BasisFN)  IMax(2)
            TYPE(BasisFN)  ISym
            INTEGER        nElec
            INTEGER        niWork(1)
!           (nEl)
!            INTEGER        nIndJ(nEl)

!           nIndJ is niWork(1+nEl)
         END TYPE
         integer, PARAMETER :: GNDWorkSize=5+3*BasisFNSize
