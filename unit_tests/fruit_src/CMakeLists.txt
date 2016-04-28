
# Libfruit provides the machinery for unit testing in Fortran

neci_add_library(
    TARGET fruit
    TYPE STATIC
    SOURCES fruit.f90
    LINKER_LANGUAGE Fortran )
