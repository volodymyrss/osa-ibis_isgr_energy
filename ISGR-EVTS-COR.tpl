################################################################################
#
#DATANAME IBIS-ISGRI - EVENTS - CORRECTED DATA
#
#DESCRIPTION
# Contains the corrected photon-by-photon events from the
# IBIS-ISGRI detector.
#
#TEMPLATE ISGR-EVTS-COR
#
#CONTENT BINARY TABLE
#
#RP obs/scw
#DB No
#VR
#
#CHANGES
# 1.4   : Template creation
# 2.1   : Added ISDCLEVL keyword
#       : Added the SELECT_FLAG column (previously it was in ISGR-EVTS-PRP)
# 2.2   : Removed EXPID keyword
#       : Added SWBOUND keyword
# 2.3   : Removed ERTFIRST, ERTLAST, BCPPID and PREVSWID keywords
#       : Template put under Configuration Control
# 7.0   : Replaced ISGRI_PI column by ISGRI_RT (SCREW-01384)
# 7.1   : Renamed ISGRI_RT column as ISGRI_PI (SCREW-01430)
# 7.2   : Including original PHA and RT
#
################################################################################
XTENSION	BINTABLE	/ Binary table extension
EXTNAME		ISGR-EVTS-COR	/ Extension name
EXTREL		'7.2'		/ ISDC release number
BASETYPE	DAL_TABLE	/ Data Access Layer base type
TELESCOP	INTEGRAL	/ Telescope or mission name
ORIGIN		ISDC		/ Origin of FITS file
INSTRUME	IBIS		/ Instrument name
DETNAM		ISGRI		/ Name of the detector layer
ISDCLEVL	COR		/ ISDC level of data processing
CREATOR		String		/ Program that created this FITS file
CONFIGUR	String		/ Software configuration
DATE		UTC_format	/ FITS file creation date
REVOL		Integer		/ Revolution number
SWID		String		/ Science Window identifier
SW_TYPE		String		/ Type of the Science Window
SWBOUND		String		/ Reason for Science Window ending
OBTSTART	OBT_format	/ OBT of the start of the Science Window
OBTEND		OBT_format	/ OBT of the end of the Science Window
	TTYPE#	ISGRI_PI	/ Corrected rise time for ISGRI
	TFORM#	1B		/ Format of column ISGRI_PI
	TTYPE#	ISGRI_ENERGY	/ Deposited energy in the ISGRI layer
	TFORM#	1E		/ Format of column ISGRI_ENERGY
	TUNIT#	keV		/ Unit of column ISGRI_ENERGY
	TTYPE#	SELECT_FLAG	/ Selection flag
	TFORM#	1B		/ Format of column SELECT_FLAG
    TTYPE#  ISGRI_PHA2   / Pulse height in the ISGRI layer, after LUT1
    TFORM#  1U      / Format of column ISGRI_PHA2

