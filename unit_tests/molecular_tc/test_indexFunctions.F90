#include "macros.h"
program test_index_functions
  use fiveIndexLMat
  use LMat_mod
  use LMat_indexing
  use procedure_pointers, only: lMatInd
  use LMat_aux
  use fruit
  implicit none

  abstract interface
     pure function six_index_func_t(a,b,c,i,j,k) result(index)
       use constants
       integer(int64), value :: a,b,c,i,j,k
       integer(int64) :: index
     end function six_index_func_t

     pure function five_index_func_t(a,b,c,i,j,k) result(index)
       use constants
       integer(int64), value :: a,b,c,i,j,k
       integer(int64) :: index
     end function five_index_func_t  
  end interface 

  call init_fruit()

  call index_functions_test_driver()
  
  call fruit_summary()
  call fruit_finalize()
  
  contains

    subroutine index_functions_test_driver()
      implicit none
      call setup_tests()
      ! test if the indexing functions are 1-to-1 and have the intended symmetry
      call lMatFiveInd_test()
      
    end subroutine index_functions_test_driver

    !------------------------------------------------------------------------------------------!

    subroutine setup_tests()
      implicit none

      ! examplary case: 14 orbitals
      nBI = 14
      call initFiveIndexAccess()

      lMatInd => lMatIndSym
    end subroutine setup_tests

    !------------------------------------------------------------------------------------------!

    subroutine lMatFiveInd_test()
      implicit none
      integer, parameter :: a = 4, b = 8, i = 3, j = 7, n = 10
      ! first test the symmetry with some arbitrary indices
      call assert_equals(LMatFiveInd(a,b,n,i,j,n),LMatFiveInd(i,b,n,a,j,n))
      call assert_equals(LMatFiveInd(a,b,n,i,j,n),LMatFiveInd(a,n,j,i,n,b))
      call assert_equals(LMatFiveInd(a,b,n,j,n,i),LMatFiveInd(j,n,n,a,b,i))
      call assert_equals(LMatFiveInd(a,b,n,j,n,i),LMatFiveInd(n,j,b,i,a,n))
      call assert_equals(LMatFiveInd(a,b,n,n,j,i),LMatFiveInd(n,i,b,a,n,j))

      ! then, test if the indexing function is 1-to-1
      call test_five_index_function(LMatFiveInd)
    end subroutine lMatFiveInd_test

    !------------------------------------------------------------------------------------------!

    subroutine test_six_index_function(indFunc)
      implicit none
      procedure(six_index_func_t) :: indFunc

      integer(int64) :: a,b,c,i,j,k
      integer(int64), allocatable :: indTable(:)

      allocate(indTable(indFunc(nBI,nBI,nBI,nBI,nBI,nBI)))


    end subroutine test_six_index_function

    !------------------------------------------------------------------------------------------!

    subroutine test_five_index_function(indFunc)
      ! generate all possible indices and check that they map the range 1-maxIndex 1-to-1
      implicit none
      procedure(five_index_func_t) :: indFunc

      integer(int64) :: a,b,i,j,n
      integer(int64) :: maxInd, cInd
      integer(int64), allocatable :: indTable(:,:)
      integer :: errSum

      maxInd = indFunc(nBI,nBI,nBI,nBI,nBI,nBI)
      ! the table of the computed indices
      allocate(indTable(maxInd,5), source = 0_int64)

      ! number of fails:
      errSum = 0

      do n = 1, nBI
         do a = 1, nBI
            do b = 1, a
               do i = 1, min(a,n)
                  do j = 1, minval((/b,n,i/))
                     cInd = indFunc(a,b,n,i,j,n)
                     if(sum(indTable(cInd,:)).eq.0) then
                        indTable(cInd,:) = (/a,b,i,j,n/)
                     else
                        print *, "Double entry with indices", indTable(cInd,:), "vs", a,b,i,j,n
                        stop
                        errSum = errSum + 1
                     endif
                  end do
               end do
            end do
         end do
      end do

      call assert_equals(errSum,0)

    end subroutine test_five_index_function


end program test_index_functions
