#!/bin/bash

for f in *; do
    if [ -d "${f}" ]; then
        cd "${f}"

        rm -f test.*
        rm -f FCIMCStats*
        rm -f INITIATORStats
        rm -f HAMIL
        rm -f TMAT
        rm -f UMAT
        rm -f NodeFile*
        rm -f Blocks_*
        rm -f ENERGIES
        rm -f MCPATHS
        rm -f RDM*
        rm -f Two*
        rm -f CLASS*
        rm -f InitPops*
        rm -f ORBOCC*
        rm -f hamil.*
        rm -f overlap.*
        rm -f lowdin.*
        rm -f gram_schmidt.*
        rm -f POPS_NORM
        rm -f DETS
        rm -f SpaceMCStats
        rm -f SPECTRAL_DATA
        rm -f FTLM_EIGV
        rm -f fciqmc_data*
        rm -f EIGV_DATA

        cd ..
    fi
done
