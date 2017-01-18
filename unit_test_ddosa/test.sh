TMPDIR=$PWD/tmp

(
    mkdir -pv $TMPDIR
    cd $TMPDIR
    rundda.py ibis_isgr_energy -m ddosa  -m /onlybright/v4.2 -a 'ddosa.ScWData(input_scwid="023900320010.001")'
)
