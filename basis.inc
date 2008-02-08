      TYPE Symmetry
         INTEGER*8 S
      END TYPE
   
      integer, PARAMETER :: SymmetrySize=2
      TYPE BasisFN
         INTEGER k(3)
         INTEGER Ms
         TYPE(Symmetry) sym
      END TYPE

      integer, PARAMETER :: BasisFNSize=SymmetrySize+4


