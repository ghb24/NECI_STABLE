C.. A test subroutine to trial a new symmetry excitation generation algorithm
C.. Algorithm: First list all the different symmetry classes in NI.
C.. For each pair of symmetry classes, determine its symmetry product.
C.. Classify and store each possible pair of orbitals under its symmetry product.
C.. With each symmetry product (denoted []), calculate []'' which is the 
C.. remainder after having removed 
      SUBROUTINE SYMGENEXCITS(NI,NEL)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER NEL,NI(NEL)
C.. List the symmetry classes in NI.
                  
      END
