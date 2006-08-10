         INCLUDE 'basis.inc'
         TYPE SymPairProd 
            TYPE(Symmetry) Sym
            INTEGER nPairs
            INTEGER nIndex
         ENDTYPE
         PARAMETER SymPairProdSize=SymmetrySize+2

         TYPE SymClass 
            TYPE(Symmetry) SymRem
            INTEGER        SymLab
         ENDTYPE
         PARAMETER SymClassSize=1+SymmetrySize
