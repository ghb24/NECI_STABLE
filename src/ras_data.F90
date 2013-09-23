#include "macros.h"

module ras_data

    use constants

    implicit none

    ! The five parameters defining a RAS space.
    type ras_parameters
        ! The number of (spatial) orbitals in each of the RAS spaces.
        integer(sp) :: size_1, size_2, size_3
        ! The minimum number of electrons (both alpha and beta) in RAS1,
        ! and the maximum number in RAS3.
        integer(sp) :: min_1, max_3
        ! The lower and upper bounds on the number of electrons that occupy
        ! be in ras1 orbitals for each string.
        integer(sp) :: lower_ras1, upper_ras1
        ! The total number of ras classes.
        integer(sp) :: num_classes
        ! cum_clasees(i) holds the cumulative number of strings in classes
        ! up to (and *not* including) class i.
        integer(sp), allocatable, dimension(:) :: cum_classes
        ! The total number of strings (there are the same number of alpha
        ! and beta strings for now).
        integer(sp) :: num_strings
        ! If you input the number of electrons in RAS1 as the first index
        ! and the number of electrons in RAS3 as the second index then
        ! this array will give you the label of corresponding class.
        integer(sp), allocatable, dimension(:,:) :: class_label
    end type

    ! A RAS class refers to a collection of strings with an fixed number
    ! of electrons in RAS1, RAS2 and RAS3. Thus, the collection of all
    ! RAS strings will, in general, be formed from many RAS classes.
    type ras_class_data
        ! The number of electrons in each of RAS1, RAS2 and RAS3.
        integer(sp) :: nelec_1, nelec_2, nelec_3
        ! The total number of strings in this class.
        integer(sp) :: class_size
        ! The vertex weights used in the addressing scheme.
        integer(sp), allocatable, dimension(:,:) :: vertex_weights
        ! The number of classes which can be combined with this one in
        ! a full determinant to give the correct overall RAS parameters.
        integer(sp) :: num_comb
        ! The labels of the classes which can be combined with this one.
        integer(sp), allocatable, dimension(:) :: allowed_combns
        ! A one-to-one map from the address of a string, as obtained by the
        ! function get_address, to one where states are sorted by their symmetry.
        integer(sp), allocatable, dimension(:) :: address_map ! (class_size)
        ! The number of strings in this class with each of the symmetry labels.
        integer(sp) :: num_sym(0:7)
        ! cum_sym(i) holds the cumulative number of strings in this class with
        ! symmetry labels up to (and *not* including) i.
        integer(sp) :: cum_sym(0:7)
    end type

    type ras_vector
        ! The elements of a single block of a vector.
        real(dp), allocatable, dimension(:,:) :: elements
    end type

    type ras_factors
        ! The elements of a single block of a vector.
        real(dp), allocatable, dimension(:) :: elements
    end type

    type direct_ci_excit
        ! The addresses of the excitations.
        integer(sp), allocatable, dimension(:) :: excit_ind
        ! The corresponding parities for the excitations.
        integer(sp), allocatable, dimension(:) :: par
        ! orb(1:2,i) holds the orbitals involved in excitation i.
        integer(sp), allocatable, dimension(:,:) :: orbs
        ! The total number of excitations.
        integer(sp) :: nexcit
    end type direct_ci_excit

    ! The number of electrons occupying one alpha or beta string. As only
    ! Ms=0 is implemented, this is just nOccAlpha=nOccBeta=nEl/2
    integer(sp) :: tot_nelec
    ! The number of spatial orbitals.
    integer(sp) :: tot_norbs
    ! The total symmetry of the Hartree-Fock state (always 0 for the only case
    ! that can be treated so far...)
    integer(sp) :: HFSym_sp

    ! For a semi-stochastic ras core space.
    type(ras_parameters) :: core_ras

end module ras_data
