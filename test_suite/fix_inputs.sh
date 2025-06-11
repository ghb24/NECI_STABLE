while read p; do
  cd $p
  NORB=$(head -30 {FCIDUMP,*FciDmp}  2>/dev/null | grep -i norb | grep -oP 'NORB=\s*\K\d+')
  pwd
  sed -i "/endsys/I i\
    NORB ${NORB}
" neci.inp
  echo $NORB
  cd $OLDPWD
done < tmpfailures
