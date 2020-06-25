#!/usr/bin/env bash

# Prettify f90 source files

fprettify $@
sed -i 's/write\s*(\(.*,.*\))/write\(\1\)/i' $@
sed -i 's/read\s*(\(.*,.*\))/read\(\1\)/i' $@
sed -i 's/allocate\s*(\(.*\))/allocate\(\1\)/i' $@
sed -i 's/associate\s*(\(.*\))/associate\(\1\)/i' $@
sed -i 's/elseif/else if/i' $@
# Ignore preprocessor #endif
sed -i 's/[^#]endif/ end if/i' $@
sed -i 's/enddo/end do/i' $@
