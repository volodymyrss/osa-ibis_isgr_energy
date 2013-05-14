################################################################################
#
#DATANAME IBIS-ISGRI - EVENTS - RAW DATA
#
#DESCRIPTION
# Contains the raw photon-by-photon events from the IBIS-ISGRI detector.
#
#TEMPLATE ISGR-EVTS-RAW
#
#CONTENT BINARY TABLE
#
#RP scw
#DB No
#V1
#
#CHANGES
# 1.4   : Template creation
# 2.0   : Added ISDCLEVL keyword
# 2.1.3 : Added SWBOUND keyword
# 2.3   : Changed APID value to 1348_or_1476 (previously 1348)
#       : Template put under Configuration Control
# 3.1   : Removed PLOBTSTR and PLOBTEND keywords (SPR-00570)
# 3.4   : Changed APID value to Integer (SCREW-00412)
#
################################################################################
XTENSION	BINTABLE	/ Binary table extension
EXTNAME		ISGR-EVTS-RAW	/ Extension name
EXTREL		'3.4'		/ ISDC release number
BASETYPE	DAL_TABLE	/ Data Access Layer base type
TELESCOP	INTEGRAL	/ Telescope or mission name
ORIGIN		ISDC		/ Origin of FITS file
INSTRUME	IBIS		/ Instrument name
DETNAM		ISGRI		/ Name of the detector layer
ISDCLEVL	RAW		/ ISDC level of data processing
CREATOR		String		/ Program that created this FITS file
CONFIGUR	String		/ Software configuration
DATE		UTC_format	/ FITS file creation date
ERTFIRST	UTC_format	/ Earth received time of the first packet
ERTLAST		UTC_format	/ Earth received time of the last packet
REVOL		Integer		/ Revolution number
SWID		String		/ Science Window identifier
SW_TYPE		String		/ Type of the Science Window
SWBOUND		String		/ Reason for Science Window ending
BCPPID		String		/ Broadcast packet pointing ID at ScW start
PREVSWID	String		/ Identification of the previous Science Window
PCKSTART	Integer		/ Packet time of the start of the Science Window
PCKEND		Integer		/ Packet time of the end of the Science Window
APID		Integer		/ Application process identifier (1348 or 1476)
SSC_BEG		Integer		/ Source sequence count start
SSC_END		Integer		/ Source sequence count end
SSC_GAP		Integer		/ Number of missing packets
	TTYPE#	DELTA_TIME	/ Delta time to previous event
	TFORM#	1B		/ Format of column DELTA_TIME
	TTYPE#	RISE_TIME	/ Event rise time
	TFORM#	1B		/ Format of column RISE_TIME
	TTYPE#	ISGRI_PHA	/ Pulse height in the ISGRI layer
	TFORM#	1U		/ Format of column ISGRI_PHA
	TTYPE#	ISGRI_Y		/ Y location in the ISGRI layer
	TFORM#	1B		/ Format of column ISGRI_Y
	TTYPE#	ISGRI_Z		/ Z location in the ISGRI layer
	TFORM#	1B		/ Format of column ISGRI_Z
