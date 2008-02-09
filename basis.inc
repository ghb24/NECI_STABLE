      TYPE Symmetry
         SEQUENCE
         INTEGER*8 S
      END TYPE
   
      integer, PARAMETER :: SymmetrySize=2
      TYPE BasisFN
         SEQUENCE
         INTEGER k(3)
         INTEGER Ms
         TYPE(Symmetry) sym
      END TYPE

      integer, PARAMETER :: BasisFNSize=SymmetrySize+4


