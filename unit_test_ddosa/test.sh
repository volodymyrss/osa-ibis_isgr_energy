TMPDIR=$PWD/tmp

(
    mkdir -pv $TMPDIR
    cd $TMPDIR
    rundda.py ibis_isgr_energy -m ddosa  -m /onlybright/v4.2 -a 'ddosa.ScWData(input_scwid="023900320010.001")'
 #   rundda.py BinBackgroundSpectrum -m ddosa -m /eddosa/GG5  -m /onlybright/v4.2 -a 'ddosa.ScWData(input_scwid="023900320010.001")' -f BinBackgroundSpectrum
)
