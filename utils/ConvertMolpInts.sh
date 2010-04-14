#!/bin/bash -x

wd=$(pwd)

grep -i final Molp.out | sed -e 's/\(^.*:\)\(.*\)/\2/' > MolpStats

`./ConvertMolpFCID.x > Convert.out`
