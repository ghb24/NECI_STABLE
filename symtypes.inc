         INCLUDE 'basis.inc'
         TYPE SymPairProd 
            TYPE(Symmetry) Sym
            INTEGER nPairs
            INTEGER nIndex
         ENDTYPE
         integer, PARAMETER  :: SymPairProdSize=SymmetrySize+2

         TYPE SymClass 
            TYPE(Symmetry) SymRem
            INTEGER        SymLab
            INTEGER        spacer
         ENDTYPE
         integer, PARAMETER :: SymClassSize=2+SymmetrySize
!  The spacer is there to make sure we have a structure which is a multiple of 8-bytes for 64-bit machines.
