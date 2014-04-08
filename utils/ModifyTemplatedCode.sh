#!/bin/bash

sed -i -e '/INTEGER iC,i,j,l,iC2/s/$/,k/' \
       -e '/realtemp_sign=transfer(temp_sign, realtemp_sign)/s/\( *\).*/\1do k=1,lenof_sign\n\1 realtemp_sign(k)=temp_sign(k)\n\1enddo/' \
       -e '/use AmpList/s/$/\n   use neci_conversion/' \
       -e '/temp_sign(1)=iSgn$/s/=.*/=convert_type(iSgn,temp_sign(1))/' \
       -e '/i=i\*E/s/=i\*\(.*\)/=i*convert_type(\1,i)/' \
       -e '/nparts=nparts+/s/=nparts+\(.*\)/=nparts+convert_type(\1,nparts)/' \
       -e '/amp=GetAmpl/s/=\(.*\)/=convert_type(\1,amp)/' \
       -e '/accum=accum+/s/=accum+\(.*\)/=accum+convert_type(\1,accum)/' \
       -e '/dCur=/s/=\(.*\)/=convert_type(\1,dCur)/' \
       CCMCClusterList.F90

var[1]="integer"; kind[1]="int32"; conv[1]="int"
var[2]="integer"; kind[2]="int64"; conv[2]="int"
var[3]="real";    kind[3]="sp";    conv[3]="real"
var[4]="real";    kind[4]="dp";    conv[4]="real"
file="src/neci/neci_conversion.F90"
cat << _EOF_ > ${file}
#if defined(MOLPRO_NECI)
module neci_conversion

 use constants, only : int32, int64, sp, dp, lenof_sign

 interface convert_type
_EOF_
for i in ${!var[*]}; do
 for j in ${!var[*]}; do
  echo "  module procedure ${kind[${i}]}_to_${kind[${j}]}" >> ${file}
  echo "  module procedure lenof_sign_${kind[${i}]}_to_${kind[${j}]}" >> ${file}
 done
done
cat << _EOF_ >> ${file}
 end interface convert_type

 contains

_EOF_
for i in ${!var[*]}; do
 for j in ${!var[*]}; do
cat << _EOF_ >> ${file}
! convert ${var[${i}]}(${kind[${i}]}) to ${var[${j}]}(${kind[${j}]})
  function ${kind[${i}]}_to_${kind[${j}]}(var_in,var_out)
  implicit none
  ${var[${j}]}(${kind[${j}]}) :: ${kind[${i}]}_to_${kind[${j}]}
  ${var[${i}]}(${kind[${i}]}) :: var_in
  ${var[${j}]}(${kind[${j}]}) :: var_out ! used for procedure matching only
  ${kind[${i}]}_to_${kind[${j}]}=${conv[${j}]}(var_in,${kind[${j}]})
  end function ${kind[${i}]}_to_${kind[${j}]}

! convert ${var[${i}]}(${kind[${i}]})(1:lenof_sign) to ${var[${j}]}(${kind[${j}]})(1:lenof_sign)
  function lenof_sign_${kind[${i}]}_to_${kind[${j}]}(var_in,var_out)
  implicit none
  ${var[${j}]}(${kind[${j}]}), dimension(lenof_sign) :: lenof_sign_${kind[${i}]}_to_${kind[${j}]}
  ${var[${i}]}(${kind[${i}]}), dimension(lenof_sign) :: var_in
  ${var[${j}]}(${kind[${j}]}), dimension(lenof_sign) :: var_out ! used for procedure matching only
  lenof_sign_${kind[${i}]}_to_${kind[${j}]}=${conv[${j}]}(var_in,${kind[${j}]})
  end function lenof_sign_${kind[${i}]}_to_${kind[${j}]}

_EOF_
 done
done
cat << _EOF_ >> ${file}
end module neci_conversion
#else
subroutine neci_conversion_dummy()
 implicit none
 return
end subroutine neci_conversion_dummy
#endif
_EOF_

