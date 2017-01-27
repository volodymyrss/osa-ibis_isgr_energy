scw=${1:?}

TMPDIR=$PWD/tmp

add2python $PWD

(
    mkdir -pv $TMPDIR
    cd $TMPDIR
    rundda.py ibis_isgr_energy -m ddosa  -m ddosa11 -m /onlybright/v4.2 -a 'ddosa.ScWData(input_scwid="'$scw'")'
 #   rundda.py BinBackgroundSpectrum -m ddosa -m /eddosa/GG5  -m /onlybright/v4.2 -a 'ddosa.ScWData(input_scwid="'$scw'")' -f BinBackgroundSpectrum
)

exit

(
    mkdir -pv $TMPDIR/osa10
    cd $TMPDIR/osa10

    source /sps/integral//data/scripts/init_integral_generic.sh 

    rundda.py ibis_isgr_energy -m ddosa  -m /onlybright/v4.2 -a 'ddosa.ScWData(input_scwid="'$scw'")'
 #   rundda.py BinBackgroundSpectrum -m ddosa -m /eddosa/GG5  -m /onlybright/v4.2 -a 'ddosa.ScWData(input_scwid="'$scw'")' -f BinBackgroundSpectrum
)
