      INCLUDE 'basis.inc'
         TYPE GNDWork
            INTEGER        NSWORK(4)
            TYPE(BasisFN)  IMax(2)
            TYPE(BasisFN)  ISym
            INTEGER        nElec
            INTEGER        niWork(1)
C           (nEl)
C            INTEGER        nIndJ(nEl)

C           nIndJ is niWork(1+nEl)
         END TYPE
         PARAMETER GNDWorkSize=5+3*BasisFNSize
