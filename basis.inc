      TYPE Symmetry
         INTEGER*8 S
      END TYPE
   
      PARAMETER SymmetrySize=2
      TYPE BasisFN
         INTEGER k(3)
         INTEGER Ms
         TYPE(Symmetry) sym
      END TYPE

      PARAMETER BasisFNSize=SymmetrySize+4


