!A simple module with some interfaces to avoid compilation warnings. Ignore.
module neci_intfce
    implicit none
    interface
        subroutine gensymexcitit2(ni, nel, g1, nbasis, tsetup, nmem, nj, ic, store, ilevel)
            use systemdata, only: basisfn
            implicit none
            integer nel, ni(nel), nbasis
            type(basisfn) g1(nbasis)
            integer store(6)
            integer, target :: nmem(*)
            integer nj(nel), ic
            logical tsetup
            integer ilevel
        end subroutine gensymexcitit2
        subroutine gensymexcitit3par(ni, tsetup, nmem, nj, ic, store, ilevel, iminelec1, imaxelec1)
            use systemdata, only: nel
            implicit none
            integer ni(nel)
            integer, pointer :: dstore(:)
            integer store(6)
            integer, target ::  nmem(*)
            integer nj(nel), ic
            logical tsetup
            integer ilevel
            integer iminelec1, imaxelec1
        end subroutine gensymexcitit3par
        subroutine genexcitprob(ni, nj, nel, niexcitor, g1, nbasismax, arr, nbasis, pgen)
            use constants, only: dp
            use systemdata, only: basisfn
            implicit none
            integer nel, ni(nel), nj(nel), nbasis, nbasismax(*)
            integer, target :: niexcitor(*)
            type(basisfn) g1(nbasis)
            real(dp) pgen
            real(dp) arr(nbasis, 2)
        end subroutine genexcitprob

        subroutine setbasislim_hub(nbasismax, nmaxx, nmaxy, nmaxz, len, &
                                   tpbc, treal)
            integer nbasismax(5, *), nmaxx, nmaxy, nmaxz, len
            logical tpbc, treal
        end subroutine

        subroutine setbasislim_hubtilt(nbasismax, nmaxx, nmaxy, nmaxz, len, &
                                       tpbc, itiltx, itilty)
            use constants, only: sizeof_int, dp
            implicit none
            integer nbasismax(5, *), nmaxx, nmaxy, nmaxz, len
            logical tpbc
            integer itiltx, itilty
        end subroutine

        subroutine calctmathub(nbasis, nbasismax, bhub, ttilt, g1, treal, tpbc)
            use constants, only: dp
            use systemdata, only: basisfn, t_open_bc_x, t_open_bc_y
            use oneeints, only: tmat2d, tmatsym, setuptmat
            use parallel_neci, only: iprocindex
            use sym_mod, only: mompbcsym
            implicit none
            integer nbasis, nbasismax(5, *)
            type(basisfn) g1(nbasis)
            real(dp) bhub
            integer isize
            logical ttilt, treal, tpbc
        end subroutine

    end interface
end module neci_intfce
