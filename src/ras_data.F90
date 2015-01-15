! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
#include "macros.h"

module ras_data

    use constants

    implicit none

    ! The five parameters defining a RAS space.
    type ras_parameters
        ! The number of (spatial) orbitals in each of the RAS spaces.
        integer :: size_1, size_2, size_3
        ! The minimum number of electrons (both alpha and beta) in RAS1,
        ! and the maximum number in RAS3.
        integer :: min_1, max_3
        ! The lower and upper bounds on the number of electrons that occupy
        ! be in ras1 orbitals for each string.
        integer :: lower_ras1, upper_ras1
        ! The total number of ras classes.
        integer :: num_classes
        ! cum_clasees(i) holds the cumulative number of strings in classes
        ! up to (and *not* including) class i.
        integer, allocatable, dimension(:) :: cum_classes
        ! The total number of strings (there are the same number of alpha
        ! and beta strings for now).
        integer :: num_strings
        ! If you input the number of electrons in RAS1 as the first index
        ! and the number of electrons in RAS3 as the second index then
        ! this array will give you the label of corresponding class.
        integer, allocatable, dimension(:,:) :: class_label
    end type

    ! A RAS class refers to a collection of strings with an fixed number
    ! of electrons in RAS1, RAS2 and RAS3. Thus, the collection of all
    ! RAS strings will, in general, be formed from many RAS classes.
    type ras_class_data
        ! The number of electrons in each of RAS1, RAS2 and RAS3.
        integer :: nelec_1, nelec_2, nelec_3
        ! The total number of strings in this class.
        integer :: class_size
        ! The vertex weights used in the addressing scheme.
        integer, allocatable, dimension(:,:) :: vertex_weights
        ! The number of classes which can be combined with this one in
        ! a full determinant to give the correct overall RAS parameters.
        integer :: num_comb
        ! The labels of the classes which can be combined with this one.
        integer, allocatable, dimension(:) :: allowed_combns
        ! A one-to-one map from the address of a string, as obtained by the
        ! function get_address, to one where states are sorted by their symmetry.
        integer, allocatable, dimension(:) :: address_map ! (class_size)
        ! The number of strings in this class with each of the symmetry labels.
        integer :: num_sym(0:7)
        ! cum_sym(i) holds the cumulative number of strings in this class with
        ! symmetry labels up to (and *not* including) i.
        integer :: cum_sym(0:7)
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
        integer, allocatable, dimension(:) :: excit_ind
        ! The corresponding parities for the excitations.
        integer, allocatable, dimension(:) :: par
        ! orb(1:2,i) holds the orbitals involved in excitation i.
        integer, allocatable, dimension(:,:) :: orbs
        ! The total number of excitations.
        integer :: nexcit
    end type direct_ci_excit

    ! The number of electrons occupying one alpha or beta string. As only
    ! Ms=0 is implemented, this is just nOccAlpha=nOccBeta=nEl/2
    integer :: tot_nelec
    ! The number of spatial orbitals.
    integer :: tot_norbs
    ! The total symmetry of the Hartree-Fock state (always 0 for the only case
    ! that can be treated so far...)
    integer :: HFSym_ras

    ! For a semi-stochastic ras space.
    type(ras_parameters) :: core_ras
    ! For a trial-wavefunction ras space.
    type(ras_parameters) :: trial_ras

end module ras_data
