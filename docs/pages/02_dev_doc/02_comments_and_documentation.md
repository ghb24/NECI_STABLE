---
title: Comments and documentation
---

## Comments and documentation

**Comments**<br>

Comments should not explain what the code does, but why.
```Fortran
! Not recommended
    ! increment i
    i = i + 1
```
Code should be written to require as least comments as possible.
This is for example achieved by using self explaining variable names
```Fortran
! Not recommended

    ! n is the number of spatial orbitals
    n = ...

! Recommended
    n_spat_orbs = ...
```
and encoding assumptions in the type system.
```Fortran
! Not recommended

    elemental real function get_area(r)
        real, intent(in) :: r
            !! r encodes the radius of a circle.
        get_area = PI * r**2
    end function

! Recommended

    elemental real function get_area(circle)
        type(Circle_t), intent(in) :: circle
        get_area = PI * circle%r**2
    end function
```

**Documentation**<br>
The documentation of NECI is automatically generated from the
code using [ford](https://github.com/Fortran-FOSS-Programmers/ford).

From its documentation we copy the most relevant part.

In modern (post 1990) Fortran, comments are indicated by an exclamation mark (`!`).
FORD will ignore a normal comment like this.
However, comments with two exclamation marks (`!!`)
or exclamation mark and  greater sign (`!>`) are interpreted
as documentation and will be captured for inclusion in the output.

Documentating with `!!` is preferred and comes after whatever you are
documenting.
```Fortran
subroutine feed_pets(cats, dogs, food, angry)
    !! Feeds your cats and dogs, if enough food is available.
    !!
    !! If not enough food is available, some of your pets will get angry.
    integer, intent(in)  :: cats
        !! The number of cats to keep track of.
    integer, intent(in)  :: dogs
        !! The number of dogs to keep track of.
    real, intent(inout)  :: food
        !! The ammount of pet food (in kilograms) which you have on hand.
    integer, intent(out) :: angry
	    !! The number of pets angry because they weren't fed.

    !...
end subroutine feed_pets
```

To keep compatibility to existing documentation in the Doxygen format
is possible to use `!>` to document **before** whatever you are documenting.
```Fortran
!> Feeds your cats and dogs, if enough food is available.
!>
!> If not enough food is available, some of your pets will get angry.
subroutine feed_pets(cats, dogs, food, angry)
    !> The number of cats to keep track of.
    integer, intent(in)  :: cats
    !> The number of dogs to keep track of.
    integer, intent(in)  :: dogs
    !> The ammount of pet food (in kilograms) which you have on hand.
    real, intent(inout)  :: food
	!> The number of pets angry because they weren't fed.
    integer, intent(out) :: angry

    !...
end subroutine feed_pets
```

Note that unlike Doxygen `intent`, `allocatable`, `target` etc. attributes
are automatically
parsed by `ford` and **shall not be** specified in the documentation redundantly.

Ford has some useful notes you may recognise from using Doxygen, like

@note
notes
@endnote

@warning
warning tags
@endwarning

@todo
todo tags
@endtodo

@bug
bug tags
@endbug
