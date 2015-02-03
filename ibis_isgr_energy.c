/************************************************************************
 * FILE:        ibis_isgr_energy.c
 * VERSION:     8.2
 * COMPONENT:   ibis_isgr_energy
 * AUTHOR:      P. Laurent,     SAp-CEA & APC
 * DESCRIPTION: component subroutine source code
 * HISTORY:
 *   NB: previous versions are made by A.Sauvageon
 *   MD,     28/03/1999, template version
 *   MD,     09/04/1999, demi-module version
 *   PL,     20/05/1999, first module version
 *   PL,     09/06/1999, add lecture of the gain-offset file
 *   PL,     09/06/1999, add a column in the output file for corrected energy
 *   MD,       /11/1999, fixed bug of index calculation in ibis_isgr_energy_corr
 *   PL,     14/02/2000, add gain-offset for the rise-time values
 *   PL,     07/03/2000, add STAMP for corrected data 
 *   PL,     23/06/2000, implement the modification of the RISE-MOD file,
 *                       (6 columns instead of 260), add error codes,
 *                       try to exit correctly when no input are present
 *   PL,     18/12/2000, try again to exit correctly when no input are present,
 *                       and to add error codes
 *   PL,     18/01/2001, replace RETURN command by EXIT in accordance to CTS
 *   PL,     21/03/2001, try again and again to exit correctly when no input
 *                       are present,
 *   PL,     27/03/2002, change energy correction algorithm, global rise-time
 *                       correction for different pixel type, add pixel type
 *                       in GO calibration file
 *   PL,     10/04/2002, correct bug in correction algorithm
 *   PL,     22/04/2002, add new groups for IC files
 *   PL,     31/07/2002, correct list of copied KEYWORDS
 *   PL,3.6, 19/09/2002, SINT32BITPR 1739 (modify the random generator initialisation)
 *   PL,4.0, 30/07/2003, remove random generator no longer needed, change
 *                       correction algorithm for new version IC files (SCREW 1229)
 *  ASA,4.1, 25/08/2003, SPR 3137 corrected: errors in formula
 *                       if rise-time or pha out of range (0-255, 0-2047)
 *                       last available LUT2 is applied (instead of LUT1 energy)
 *   PL,4.2, 06/10/2003, SCREW 1283: get ISGRI mean temperature from HK to
 *                       include temperature variation for gain and offset
 *
 *   code into C, ASA, december 2003
 *   all calculations are now made in double, finally written into float
 *   (ener is double to make cast to int identical Solaris/Linux)
 *
 *  ASA,5.0, 30/03/2004, SPR 3418 (partly) SCREW 1398 (new parameters)
 *                       SCREW 1384 (column ISGRI_RT)
 *  ASA,5.1, 08/04/2004, SPR 3418 (header) and new PHA gain laws
 *                       SCREW 1430 (column back to ISGRI_PI)
 *  ASA,5.2, 25/05/2004, SPR 3560, 3561 (misleading error message)
 *  ASA,5.3, 16/06/2004, SPR 3686 (don't calculate means on whole ScW)
 *
 *  ASA,5.4, 31/03/2005, SCREW 1499 (read parameters from IC file)
 *  ASA,5.5, 26/04/2005, SPR 4112 (offset decrease, better gain decrease)
 *  ASA, 5.5.1, 09/2005, SPR 4263 (replace Warn by Log) SPR 4138
 *
 *  New version 5.6.0 made by Christophe Couvreur (CC)
 *  CC,5.6.1,28/07/2006,   SPR 4537 (ibis_isgr_energy uses too much memory)
 *  CC, 5.6.2, 02/08/2006, SPR 4549 (in function ibis_isgr_energyTransform
 *                         index for tables s_gh and s_oh is not checking and 
 *                         it can be out of bound) 
 *  CC, 5.6.3, 19/09/2006, SPR 4579 (ibis_isgri_energy does not interpolate switch on time)
 *  CC, 5.6.4, 19/09/2006, SPR 4579 (ibis_isgri_energy does not interpolate switch on time)
 *                         correction into routine ibis_isgr_switchDetection:  
 *                         if( (*switchOnTime) - tstart > 3.0 ) modified into 
 *                         if( (tstart - (*switchOnTime)) > 3.0 ).
 *  CC, 5.7
 * ASA, 5.7                SPR 4641 + cosmetic changes
 *  CC, 5.7.1              SPR 4664 (correction of energy drift around 23kev)
 *  CC, 6.0, 12/07/2007,   SPR 4695 (calibration is based on 2 laws to avoid increase 
 *                         of Crab flux in 20-40kev). 
 *                         Below 50kev we apply the first law and above, 
 *                         we apply the second one. In the first law, Ph gain and Offset 
 *                         don't depend on RT contrary to the second law.
 *                         - Switch on Time routine is no more considerated.
 *                         - The Table for rt effect has disapeared and included into
 *                           PH gain and offset for 2nd law.
 *  CC, 6.1, 19/07/2007,   SPR 4696 (correction bug for E=50kev)
 * ASA, 6.2, 20/12/2007,   SPR 4773 4698 (+ cosmetic changes)
 *  NP, 6.3, 24/02/2009,   SPR 4830 (interpolation for wiggles)
 *  NP, 6.4, 27/02/2009,   SCREW 2117 extrapolation for end of validity of dose
 * ASA, 6.5  16/03/2009,   SPR 4838 (+ cosmetic changes = ibis_comp_energy)
 *  NP, 6.6  14/04/2009,   SPR 4855, revert SPR 4830
 *  NP, 7.0  17/12/2009,   SCREW 2178
 *  PL, 8.0  02/02/2012,   remove IREM counters, Temperature correction by MDU
 *  NP, 8.1  21/02/2012,   completely remove IREM counters
 *  PL, 8.2  02/04/2012,   modify coefficients PAR1_..._corrPH1 to correct <50 keV behavior
 ******************************************************************************/

#include "ibis_isgr_energy.h"

static const double DtempH1[8] = {0.43, -0.39, -0.77, 0.84, -0.78, 1.09, -0.08, -0.31};
             /* delta from ISGRI mean temperature, to check probe temp2 is OK */

/************************************************************************
 * FUNCTION:  ibis_isgr_energyWork
 * DESCRIPTION:
 *  Does the work, i.e. inputs data, reads calibration data, 
 *  computes corrected energies. Returns ISDC_OK if everything is fine,
 *  else returns an error code transmitted from the called function.
 * ERROR CODES:
 *  DAL error codes
 *  I_ISGR_ERR_MEMORY             Memory allocation error
 *  ibis_isgr_energyCheckOut()    error codes
 *  ibis_isgr_energyCheckInNEW()  error codes
 *
 * PARAMETERS:
 *  workGRP   dal_element *     in  DOL of the working group
 *  gti              int        in  0 to allow PRP with 0 rows (no OBT), 1 otherwise
 *  erase            int        in  0 to replace output columns, 1 to erase rows
 *  chatter          int        in  verbosity level
 *  acorName         char *     in  acorDOL string
 *  riseName         char *     in  riseDOL string
 *  phGainDOLstr     char *     in  phGainDOL string for 2nd calibration law
 *  phOffDOLstr      char *     in  phOffsetDOL string for 2nd calibration law
 * RETURN:            int     current status
 ************************************************************************/
int ibis_isgr_energyWork(dal_element *workGRP,
                         int          gti,
                         int          erase,
                         int          chatter,
                         char        *acorName,
                         char        *riseName,
                         char        *phGainDOLstr,
                         char        *phOffDOLstr,
                         int          corGainDrift)
{
  int    i, status = ISDC_OK,
         freeStatus= ISDC_OK;
  long   numEvents = 0l,      /* number of S1 events   */
         buffSize,            /* size of output buffers in bytes */
         revol;               /* to calculate the PHA gain decrease */
  char  *logString=NULL,
         hkName[11];
  double meanT[8],               /* mean of the 8 mce temperatures */
    meanBias[8] ;

  /* input data */
  DAL3_Word *isgriPha = NULL; /* ISGRI_PHA */
  DAL3_Byte *riseTime = NULL; /* RISE_TIME */
  DAL3_Byte *isgriY   = NULL; /* ISGRI_Y   */
  DAL3_Byte *isgriZ   = NULL; /* ISGRI_Z   */

  /* output data */
  DAL3_Word *isgrPHA2 = NULL;   /* ISGRI_PI, now RT corr */
  DAL3_Byte *isgrPi = NULL;   /* ISGRI_PI, now RT corr */
  DAL3_Word *isgrPHA1 = NULL;   /* ISGRI_PI, now RT corr */
  DAL3_Byte *isgrRT1 = NULL;   /* ISGRI_PI, now RT corr */
  float *isgrEnergy = NULL;   /* ISGRI_ENERGY  */
  OBTime obtStart=DAL3_NO_OBTIME,
         obtEnd=DAL3_NO_OBTIME;

  /* data from ISGRI_GO */
  /* needs to be double otherwise 1E-6 differences */
  double **goTab = NULL;
  int     *pixTab= NULL;

  /* data from ISGRI-3DL2-MOD */
  short *isgriRtTab=NULL;
  /* new LUT2 3D, SPR 4537,before: float ***isgriRtTab=NULL; */

  /*double *s_gh=NULL,*s_oh=NULL;*/ /*coefficients for rt effect*/

  /*input parameters*/
  dal_element *in_table = NULL,
              *outTable = NULL;
  dal_element *isgrOffsTabPtr = NULL;
  dal_element *isgrRiseTabPtr = NULL;
  dal_element *isgrHK1_Ptr = NULL;
  dal_element *supGainTabPtr=NULL;
  dal_element *supOffsTabPtr=NULL;

  /* coefficients for 2nd calibration method */
  double *supCoeffg=NULL, *supCoeffo=NULL;

  do {

    if ((logString=(char *)calloc(DAL_BIG_STRING, sizeof(char))) == NULL) {
      status=I_ISGR_ERR_MEMORY;
      break;
    }
    /*#################################################################*/
    /* Locate ISGR-EVTS-RAW bintable, for info: to know if ALL is used */
    /*#################################################################*/
    if (chatter > 3) {
      status=DALobjectFindElement(workGRP, DS_ISGR_RAW, &in_table, status);
      status=DALtableGetNumRows(in_table, &numEvents, status);
      if (status != ISDC_OK) {
        RILlogMessage(NULL, Log_0,
                      "%13s bintable NOT found. Reverting to ISDC_OK", DS_ISGR_RAW);
        status=ISDC_OK;
        in_table=NULL;
      }
      else
        RILlogMessage(NULL, Log_0, "%13s bintable found: %9ld rows.",
                                  DS_ISGR_RAW, numEvents);
    }

    /*#################################################################*/
    /* Read S1 RAW data */
    /*#################################################################*/
    status=ibis_isgr_energyReadData(workGRP, gti, chatter,
                                    &numEvents, &revol,
                                    &obtStart, &obtEnd,
                                    &isgriPha,
                                    &riseTime,
                                    &isgriY, &isgriZ);
    if ((status != ISDC_OK) || (numEvents == 0l))
      break;
    if ((chatter > 2) && (gti))
      RILlogMessage(NULL, Log_0, "OBT range: %020lld , %020lld", obtStart, obtEnd);
    if ((obtStart < 0) || (obtEnd < 0)) {
      if (gti) {
        RILlogMessage(NULL, Warning_1, "At least one OBT limit is negative.");
        RILlogMessage(NULL, Warning_1, "Using all ScW to calculate mean bias and temperature.");
      }
      obtStart=DAL3_NO_OBTIME;
      obtEnd=DAL3_NO_OBTIME;
    }

    /*#################################################################*/
    /* Check if output structures exist */
    /*#################################################################*/
    strncpy(logString, DS_ISGR_RAW, 10);
    logString[10]='\0';
    strcat(logString, "COR");
    status=ibis_isgr_energyCheckOut(workGRP, logString, chatter,
                                    numEvents, erase, &outTable);
    if (status != ISDC_OK) break;

    /*#################################################################*/
    /* Check input data */
    /*#################################################################*/
    if (chatter > 3) RILlogMessage(NULL, Log_0, "Opening 4 calibration tables...");
    status=ibis_isgr_energyCheckIn(acorName, riseName,
                                      phGainDOLstr, phOffDOLstr,
                                      &isgrOffsTabPtr,
                                      &isgrRiseTabPtr,
                                      &supGainTabPtr,
                                      &supOffsTabPtr);
    if (status != ISDC_OK) break;

    /*#################################################################*/
    /*Allocate memory for correction tables */
    /*#################################################################*/
    status=I_ISGR_ERR_MEMORY;

    /*------------- memory allocation for lut1----------------------*/
    goTab=(double **)calloc(ISGRI_GO_N_COL-1,sizeof(double *));
    if (goTab == NULL) break;
    for (i=0; i<ISGRI_GO_N_COL-1; i++) goTab[i]=NULL;
    /* to correctly free memory if next loop is broken */
    for (i=0; i<ISGRI_GO_N_COL-1; i++) {
      goTab[i]=(double *)calloc(ISGRI_N_PIX,sizeof(double));
      if (goTab[i] == NULL) break;
    }
    if (i < ISGRI_GO_N_COL-1) break;

    pixTab=(int *)calloc(ISGRI_N_PIX,sizeof(int));
    if (pixTab == NULL) break;

    /*------------- memory allocation for new lut2----------------------*/
    buffSize=ISGRI_RT_N_ENER_SCALED*ISGRI_RT_N_DATA*ISGRI_RT_N_RANDOM_DIM;
    isgriRtTab=(short *)calloc(buffSize,sizeof(short));
    if (isgriRtTab == NULL) break;

    /*----------- Allocate memory for 2nd correction table -----------*/
    supCoeffg=(double*)calloc(ISGRI_RT_N_DATA*ISGRI_PHG2_N_COL,sizeof(double));
    if (supCoeffg == NULL) break;
    supCoeffo=(double*)calloc(ISGRI_RT_N_DATA*ISGRI_PHO2_N_COL,sizeof(double));
    if (supCoeffo == NULL) break;
    status=ISDC_OK;

    /*#################################################################*/
    /* Allocate memory for the output vectors */
    /*#################################################################*/
    buffSize=numEvents*sizeof(DAL3_Word);
    status=DALallocateDataBuffer((void **)&isgrPHA2, buffSize, status);
    buffSize=numEvents*sizeof(DAL3_Word);
    status=DALallocateDataBuffer((void **)&isgrPHA1, buffSize, status);
    buffSize=numEvents*sizeof(DAL3_Byte);
    status=DALallocateDataBuffer((void **)&isgrRT1, buffSize, status);
    buffSize=numEvents*sizeof(DAL3_Byte);
    status=DALallocateDataBuffer((void **)&isgrPi, buffSize, status);
    buffSize=numEvents*sizeof(float);
    status=DALallocateDataBuffer((void **)&isgrEnergy, buffSize, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Cannot allocate buffers for output data");
      break;
    }

    /*#################################################################*/
    /* Read ISGRI gain, offset and rise time correction tables */
    /*#################################################################*/
    status=ibis_isgr_energyReadCal(isgrOffsTabPtr,
                                   isgrRiseTabPtr,
                                   supGainTabPtr,
                                   supOffsTabPtr,
                                   chatter,
                                   goTab, pixTab,
                                   isgriRtTab,
                                   phGainDOLstr, phOffDOLstr,
                                   supCoeffg, supCoeffo);
    if (status != ISDC_OK) break;

    /*#################################################################*/
    /* Locate bintable for ISGRI temperature and bias */
    /*#################################################################*/
    status=DALobjectFindElement(workGRP, DS_ISGR_HK, &isgrHK1_Ptr, status);
    if (status != ISDC_OK)
      RILlogMessage(NULL, Warning_2, "%13s bintable NOT found.", DS_ISGR_HK);
    else if (chatter > 2)
      RILlogMessage(NULL, Log_0, "%13s bintable found.", DS_ISGR_HK);

    for (i=0; i<8; i++) 
      {
	meanBias[i]=KEY_DEF_BIAS;
	meanT[i]   =KEY_DEF_TEMP;
      }

    /*#################################################################*/
    /* Calculation of Temperature*/
    /*#################################################################*/
    status=ibis_energyIsgrHkCal(workGRP, obtStart, obtEnd,
                                meanT, meanBias, chatter, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Warning_2, "Reverting from status=%d to ISDC_OK",
                                    status);
      RILlogMessage(NULL, Warning_2, "Using constant ISGRI temperature and bias (%+6.2f %+6.1f)",
                                    KEY_DEF_TEMP, KEY_DEF_BIAS);
    }
    else if (chatter > 3) {
      RILlogMessage(NULL, Log_0, "Mean ISGRI module bias (V):");
      strcpy(logString, "");
      for (i=0; i<8; i++) {
        sprintf(hkName, " %+6.1f", meanBias[i]);
        strcat(logString, hkName);
      }
      RILlogMessage(NULL, Log_0, logString);
      RILlogMessage(NULL, Log_0, "Mean ISGRI module Temperature (C):");
      strcpy(logString, "");
      for (i=0; i<8; i++) {
        sprintf(hkName, " %+6.1f", meanT[i]);
        strcat(logString, hkName);
      }
      RILlogMessage(NULL, Log_0, logString);
    }
    status=ISDC_OK;

    /*#################################################################*/
    /* Perform energy correction */
    /*#################################################################*/
    status=ibis_isgr_energyTransform(outTable, numEvents, revol,
                                 goTab, pixTab,isgriRtTab,
                                 supCoeffg, supCoeffo,
                                 chatter,
                                 meanT,
				 meanBias,    
                                 isgriPha, riseTime,
                                 isgriY, isgriZ,
                                 isgrPi, isgrPHA2, isgrRT1, isgrPHA1, isgrEnergy, corGainDrift);
  } while(0);

  /*#################################################################*/
  /* Close structures opened by  DALobjectOpen() */
  /*#################################################################*/
  /* if error while closing, status is only changed if it was ISDC_OK */

  if (isgrOffsTabPtr != NULL) {
    i=DALobjectClose(isgrOffsTabPtr, DAL_SAVE, ISDC_OK);
    if (i != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Error %d when closing structure: %13s",
                                  i, DS_ISGR_GO);
      if (status == ISDC_OK) status=i;
      else RILlogMessage(NULL, Error_1,
                              "Error not taken into account for exit value");
    }
  }
  if (isgrRiseTabPtr != NULL) {
    i=DALobjectClose(isgrRiseTabPtr, DAL_SAVE, ISDC_OK);
    if (i != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Error %d when closing structure: %13s",
                                  i, DS_ISGR_3DL2_MOD);
      if (status == ISDC_OK) status=i;
      else RILlogMessage(NULL, Error_1,
                              "Error not taken into account for exit value");
    }
  }

  if (supGainTabPtr != NULL) {
    i=DALobjectClose(supGainTabPtr, DAL_SAVE, ISDC_OK);
    if (i != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Error %d when closing structure: %13s",
		    i, DS_PHG2);
      if (status == ISDC_OK) status=i;
      else RILlogMessage(NULL, Error_1,
			 "Error not taken into account for exit value");
    }
  }
  
  if (supOffsTabPtr != NULL) {
    i=DALobjectClose(supOffsTabPtr, DAL_SAVE, ISDC_OK);
    if (i != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Error %d when closing structure: %13s",
		    i, DS_PHO2);
      if (status == ISDC_OK) status=i;
      else RILlogMessage(NULL, Error_1,
			 "Error not taken into account for exit value");
    }
  }

  /*#################################################################*/
  /* Release data buffers */
  /*#################################################################*/
  if (logString != NULL) free(logString);
  if (goTab != NULL) {
    for (i=0; i<ISGRI_GO_N_COL-1; i++) if (goTab[i] != NULL) free(goTab[i]);
    free(goTab);
  }
  if (pixTab != NULL) free(pixTab);
  if (isgriRtTab != NULL) free(isgriRtTab);
  if(supCoeffg != NULL) free(supCoeffg);
  if(supCoeffo != NULL) free(supCoeffo);

  if (isgriPha != NULL)
    freeStatus=DALfreeDataBuffer((void *)isgriPha, ISDC_OK);
  if (riseTime != NULL)
    freeStatus=DALfreeDataBuffer((void *)riseTime, ISDC_OK);
  if (isgriY != NULL)
    freeStatus=DALfreeDataBuffer((void *)isgriY, ISDC_OK);
  if (isgriZ != NULL)
    freeStatus=DALfreeDataBuffer((void *)isgriZ, ISDC_OK);

  if (isgrPi != NULL)
    freeStatus=DALfreeDataBuffer((void *)isgrPi, ISDC_OK);
  if (isgrEnergy != NULL)
    freeStatus=DALfreeDataBuffer((void *)isgrEnergy, ISDC_OK);

  if (freeStatus != ISDC_OK) {
    RILlogMessage(NULL, Error_1, "Error when freeing data buffers: %d",
                                freeStatus);
    RILlogMessage(NULL, Error_1, "Error not taken into account for exit value");
  }
  else if (chatter > 2) RILlogMessage(NULL, Log_0, "Data buffers released.");

  return status;
}


/************************************************************************
 * FUNCTION:  ibis_isgr_energyCheckIn
 * DESCRIPTION:
 *  Checks if IC tables exist and have a correct format.
 *  Returns ISDC_OK if everything is fine, else returns an error code.
 * ERROR CODES:
 *  DAL error codes
 *  I_ISGR_ERR_ISGR_OFFS_BAD       ISGR-OFFS-MOD table wrong size
 *  I_ISGR_ERR_ISGR_RISE_BAD       ISGR-3DL2-MOD table wrong size 
 *  I_ISGR_ERR_IBIS_IREM_BAD       IBIS-IREM-CAL table wrong size
 *  I_ISGR_ERR_ISGR_PHGO2_BAD      other IC tables wrong size
 *
 PARAMETERS:
 *  acorName              char *  in   acorDOL string
 *  riseName              char *  in   riseDOL string
 *  phGainDOLstr          char *  in   Ph gain DOL string for 2nd calibration law
 *  phOffDOLstr           char *  in   Ph offset DOL string for 2nd calibration law
 *  isgrOffsTabPtr dal_element ** out  indicates ISGR-OFFS-MOD
 *  isgrRiseTabPtr dal_element ** out    "       ISGR-3DL2-MOD
 *  supGainTabPtr  dal_element ** out  Coeff for gain for the 2nd law
 *  supOffsTabPtr  dal_element ** out  Coeff for offset for the 2nd law
 * RETURN:            int     current status
 ************************************************************************/
int ibis_isgr_energyCheckIn(
                         char         *acorName,
                         char         *riseName,
                         char         *phGainDOLstr,
                         char         *phOffDOLstr,
                         dal_element **isgrOffsTabPtr,
                         dal_element **isgrRiseTabPtr,
                         dal_element **supGainTabPtr,
                         dal_element **supOffsTabPtr)
{
  int  status = ISDC_OK,
       numCol,
       numAxes; /*for new lut2*/
  long numRow,
       dimAxes[ISGRI_DIM_LUT2_3D]; /*for new lut2*/
  char keyVal[DAL_MAX_STRING];
  dal_dataType type;

  do {

    /*#################################################################*/
    /* Locate the bintable of the ISGRI gain-offset data */
    /*#################################################################*/
    status=DALobjectOpen(acorName, isgrOffsTabPtr, status);
    status=DALelementGetName(*isgrOffsTabPtr, keyVal, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "%13s bintable cannot be opened. Status=%d",
                                  DS_ISGR_GO, status);
      break;
    }
    if (strcmp(keyVal, DS_ISGR_GO)) {
      RILlogMessage(NULL, Error_2, "File (%s) should be a %13s",
                                  acorName, DS_ISGR_GO);
      status=I_ISGR_ERR_BAD_INPUT;
      break;
    }
    status=DALtableGetNumRows(*isgrOffsTabPtr, &numRow, status);
    status=DALtableGetNumCols(*isgrOffsTabPtr, &numCol, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Cannot get size of table %13s. Status=%d",
                                  DS_ISGR_GO, status);
      break;
    }
    if (numRow != ISGRI_N_PIX) {
      status=I_ISGR_ERR_ISGR_OFFS_BAD;
      RILlogMessage(NULL, Error_2, "Wrong number of rows (%ld) in %13s.",
                                  numRow, DS_ISGR_GO);
      break;
    }
    if (numCol != ISGRI_GO_N_COL) {
      status=I_ISGR_ERR_ISGR_OFFS_BAD;
      RILlogMessage(NULL, Error_2, "Wrong number of columns (%d) in %13s.",
                                  numCol, DS_ISGR_GO);
      break;
    }

    /*#################################################################*/
    /* Locate the bintable of the ISGRI rise-time calibration data */
    /*#################################################################*/
    status=DALobjectOpen(riseName, isgrRiseTabPtr, status);
    status=DALelementGetName(*isgrRiseTabPtr, keyVal, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "%13s image cannot be opened. Status=%d",
                                  DS_ISGR_3DL2_MOD, status);
      break;
    }
    if (strcmp(keyVal, DS_ISGR_3DL2_MOD)) {
      RILlogMessage(NULL, Error_2, "File (%s) should be a %13s",
                                  riseName, DS_ISGR_3DL2_MOD);
      status=I_ISGR_ERR_BAD_INPUT;
      break;
    }
    status=DALarrayGetStruct(*isgrRiseTabPtr, &type, &numAxes, dimAxes, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Cannot get the 3 sizes of array %13s. Status=%d",
                                  DS_ISGR_3DL2_MOD, status);
      break;
    }
    if (numAxes != ISGRI_DIM_LUT2_3D) {
      status= I_ISGR_ERR_ISGR_RISE_BAD;
      RILlogMessage(NULL, Error_2, "%13s image must be a 3D array.", DS_ISGR_3DL2_MOD);
      break;
    }
    if (  (dimAxes[0]!=ISGRI_RT_N_ENER_SCALED) || (dimAxes[1]!=ISGRI_RT_N_DATA)
        ||(dimAxes[2]!=ISGRI_RT_N_RANDOM_DIM)) {
      status=I_ISGR_ERR_ISGR_RISE_BAD;
      RILlogMessage(NULL, Error_2, "%13s array dimensions must be: %d*%d*%d",
                                  DS_ISGR_3DL2_MOD, ISGRI_RT_N_ENER_SCALED,
                                  ISGRI_RT_N_DATA, ISGRI_RT_N_RANDOM_DIM);
      break;
    }

    /*#################################################################*/
    /* Locate the bintable of the correction coefficients for gain and offset */
    /*#################################################################*/

    status=DALobjectOpen(phGainDOLstr, supGainTabPtr, status);
    //    status=DALelementGetName(*coeffGainTabPtr, keyVal, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "%13s bintable for 2nd method gain cannot be opened. Status=%d",
		    DS_PHG2, status);
      break;
    }
    status=DALtableGetNumRows(*supGainTabPtr, &numRow, status);
    status=DALtableGetNumCols(*supGainTabPtr, &numCol, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Cannot get size of table for 2nd method gain. Status=%d",
		    status);
      break;
    }
    if (numRow != ISGRI_RT_N_DATA) {
      status=I_ISGR_ERR_ISGR_PHGO2_BAD;
      RILlogMessage(NULL, Error_2, "Wrong number of rows (%ld) in gain table %13s.",
		    numRow, DS_PHG2);
      break;
    }
    if (numCol != ISGRI_PHG2_N_COL) {
      status=I_ISGR_ERR_ISGR_PHGO2_BAD;
      RILlogMessage(NULL, Error_2, "Wrong number of columns (%d) in gain table %13s.",
		    numCol, DS_PHG2);
      break;
    }
    
    status=DALobjectOpen(phOffDOLstr, supOffsTabPtr, status);
//    status=DALelementGetName(*supOffsTabPtr, keyVal, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "%13s bintable for 2nd method offset cannot be opened. Status=%d",
		    DS_PHO2, status);
      break;
    }
    status=DALtableGetNumRows(*supOffsTabPtr, &numRow, status);
    status=DALtableGetNumCols(*supOffsTabPtr, &numCol, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Cannot get size of table for 2nd method offset. Status=%d",
		    status);
      break;
    }
    if (numRow != ISGRI_RT_N_DATA) {
      status=I_ISGR_ERR_ISGR_PHGO2_BAD;
      RILlogMessage(NULL, Error_2, "Wrong number of rows (%ld) in offset table %13s.",
		    numRow, DS_PHO2);
      break;
    }
    if (numCol != ISGRI_PHO2_N_COL) {
      status=I_ISGR_ERR_ISGR_PHGO2_BAD;
      RILlogMessage(NULL, Error_2, "Wrong number of columns (%d) in offset table %13s.",
		    numCol, DS_PHO2);
      break;
    }
    
  } while(0);
  return status;
}


/************************************************************************
 * FUNCTION:  ibis_isgr_energyReadData
 * DESCRIPTION:
 *  Reads the S1 RAW events bintables.
 *  Returns ISDC_OK if everything is fine, else returns an error code.
 * ERROR CODES:
 *  DAL error codes
 *
 * PARAMETERS:
 *  workGRP  dal_element *    in   DOL of the working group
 *  gti              int      in   0 to allow PRP with 0 rows (no OBT), 1 otherwise
 *  chatter          int      in   verbosity level
 *  numEvents      long  *   out   length of the event list
 *  revol          long  *   out   revolution number
 *  obtStart     OBTime  *  in/out start of time range
 *  obtEnd       OBTime  *  in/out end   of time range
 *  isPha     DAL3_Word **   out   ISGRI_PHA 
 *  riseT     DAL3_Byte **   out   RISE_TIME 
 *  isY       DAL3_Byte **   out   ISGRI_Y 
 *  isZ       DAL3_Byte **   out   ISGRI_Z 
 * RETURN:            int     current status
 ************************************************************************/
int ibis_isgr_energyReadData(dal_element *workGRP,
                         int         gti,
                         int         chatter,
                         long       *numEvents,
                         long       *revol,
                         OBTime     *obtStart,
                         OBTime     *obtEnd,
                         DAL3_Word **isPha,
                         DAL3_Byte **riseT,
                         DAL3_Byte **isY,
                         DAL3_Byte **isZ)
{
  int    status = ISDC_OK;
  short  selected=0;
  long   buffSize;
  dal_dataType type;
  ISDCLevel    myLevel;

  do {

    status=DALattributeGetInt(workGRP, "REVOL", revol, NULL, NULL, status);

    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Cannot get attribute REVOL in input group.");
      break;
    }

  /*#################################################################*/
  /* Select and get RAW events */
  /*#################################################################*/
    selected=1;
    if (gti) myLevel=PRP;  else myLevel=RAW;
    status=DAL3IBISselectEvents(workGRP, ISGRI_EVTS, myLevel, gti,
                                obtStart, obtEnd, NULL, status);
    status=DAL3IBISgetNumEvents(numEvents, status);
    /* dal3gen>5.0.0 do not error any more with: -2504=DAL_TABLE_HAS_NO_ROWS */
    if ((status == DAL3IBIS_NO_IBIS_EVENTS) || (status == DAL_TABLE_HAS_NO_ROWS)) {
      RILlogMessage(NULL, Warning_1, "Reverting from status=%d to ISDC_OK", status);
      RILlogMessage(NULL, Warning_1, "NO event selected. Execution stopped.");
      status=ISDC_OK;
      *numEvents=0l;
      break;
    }
    else if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Events selection failed. Status=%d",
                                  status);
      break;
    }
    /* test useful when: dal3gen>5.0.0 or RAW <> PRP (0 rows, useGTI=yes) */
    if (*numEvents == 0l) {
      RILlogMessage(NULL, Warning_1, "NO event selected. Execution stopped.");
      break;
    }
    if (chatter > 0)
      RILlogMessage(NULL, Log_2, "Number of selected events: %9ld", *numEvents);

  /*#################################################################*/
  /* Allocate memory buffers */
  /*#################################################################*/
    buffSize= *numEvents * sizeof(DAL3_Word);
    status=DALallocateDataBuffer((void **)isPha, buffSize, status);
    buffSize= *numEvents * sizeof(DAL3_Byte);
    status=DALallocateDataBuffer((void **)riseT, buffSize, status);
    status=DALallocateDataBuffer((void **)isY,   buffSize, status);
    status=DALallocateDataBuffer((void **)isZ,   buffSize, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Cannot allocate buffers for input data");
      break;
    }

  /*#################################################################*/
  /* Read RAW data */
  /*#################################################################*/
    type=DAL_USHORT;
    status=DAL3IBISgetEvents(ISGRI_PHA, &type, (void *)*isPha, status);
    type=DAL_BYTE;
    status=DAL3IBISgetEvents(RISE_TIME, &type, (void *)*riseT, status);
    type=DAL_BYTE;
    status=DAL3IBISgetEvents(ISGRI_Y,   &type, (void *)*isY,   status);
    type=DAL_BYTE;
    status=DAL3IBISgetEvents(ISGRI_Z,   &type, (void *)*isZ,   status);
    if (gti) {
      type=DAL3_OBT;
      status=DAL3IBISgetEventsBins(OB_TIME, &type, 1,1, obtStart, status);
      buffSize= *numEvents;
      status=DAL3IBISgetEventsBins(OB_TIME, &type, buffSize,buffSize,
                                   obtEnd, status);
    }
    else {
      *obtStart=DAL3_NO_OBTIME;
      *obtEnd=DAL3_NO_OBTIME;
    }
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Cannot get input data");
      break;
    }

  } while(0);
  if (selected > 0) status=DAL3IBIScloseEvents(status);
  return status;
}


/************************************************************************
 * FUNCTION:  ibis_isgr_energyCheckOut
 * DESCRIPTION:
 *  Checks for the presence and size of the output columns, add rows if
 * necessary. Copy the necessary attributes.
 *  Returns ISDC_OK if everything is fine, else returns an error code.
 * ERROR CODES:
 *  DAL error codes
 *  I_ISGR_ERR_ISGR_OUT_COR   if number of output rows inconsistent
 *
 * PARAMETERS:
 *  workGRP  dal_element *    in    DOL of the working group
 *  outName         char *    in    bintable name of the output data
 *  chatter          int      in    verbosity level
 *  numEvents       long      in    length of the event list
 *  erase            int      in    0 to replace output columns, 1 to erase rows
 *  outTable dal_element ** in/out  DOL of the output table
 * RETURN:            int     current status
 ************************************************************************/
int ibis_isgr_energyCheckOut(dal_element *workGRP,
                         char         *outName,
                         int           chatter,
                         long          numEvents,
                         int           erase,
                         dal_element **outTable)
{
  int   status = ISDC_OK;
  long  outRow;

  do {

  /*#################################################################*/
  /* Locate the ISGR-EVTS-COR table */
  /*#################################################################*/
    status=DALobjectFindElement(workGRP, outName, outTable, status);
    status=DALtableGetNumRows(*outTable, &outRow, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "%13s bintable NOT found.", outName);
      break;
    }
    if (chatter > 2)
      RILlogMessage(NULL, Log_0, "%13s output table: %9ld rows.", outName, outRow);

  /*#################################################################*/
  /* Check the size of the corrected energy vectors (0 or numEvents) */
  /*#################################################################*/
    if (outRow > 0l) {

      if (erase) {
        DALtableDelRows(*outTable, 1l, outRow, status);
        if (status != ISDC_OK) {
          RILlogMessage(NULL, Error_2, "Cannot delete all rows. Status=%d",
                                      status);
          break;
        }
        outRow=0;
        if (chatter > 1)
          RILlogMessage(NULL, Log_1, "Output table: deleted OLD values.");
      }
      else if (outRow != numEvents) {
        RILlogMessage(NULL, Error_2, "%13s has wrong length (%ld rows).",
                                    outName, outRow);
        status=I_ISGR_ERR_ISGR_OUT_COR;
        break;
      }

    }
    if (outRow == 0l) {
      status=DALtableAddRows(*outTable, 0l, numEvents, status);
      if (status != ISDC_OK) {
        RILlogMessage(NULL, Error_2, "Cannot ADD rows. Status=%d", status);
        break;
      }
    }

  /*#################################################################*/
  /* Copy the required attributes (SPR  2958)                        */
  /*#################################################################*/
    status=DAL3GENattributeCopy(workGRP, *outTable,
                                "REVOL,SWID,SW_TYPE,SWBOUND,OBTSTART,OBTEND", status);
    if (status !=ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Could not copy ScW attributes");
      break;
    }

  } while(0);
  return status;
}


/************************************************************************
 * FUNCTION:  ibis_isgr_energyReadCal
 * DESCRIPTION:
 *  Reads the calibration data for the ISGRI layer.
 *  Returns ISDC_OK if everything is fine, else returns an error code.
 * ERROR CODES:
 *  DAL error codes
 *
 * PARAMETERS:
 *  isgrOffsTabPtr  dal_element *    in  ISGR-OFFS-MOD
 *  isgrRiseTabPtr  dal_element *    in  ISGR-3DL2-MOD
 *  supGainTabPtr   dal_element *    in  DOL for gain of 2nd law
 *  supOffsTabPtr   dal_element *    in  DOL for offset of 2nd law
 *  chatter                int       in  verbosity level
 *  isgriGoTab          double **   out  ISGRI gain-offset table
 *  isgriPixTab            int  *   out  ISGRI pixel type table
 *  isgriRtTab           short  *   out  ISGRI rise-time calibration table (3D)
 *  supCoeffg           double  *   out  Coefficients of gain for the 2nd law
 *  supCoeffo           double  *   out  Coefficients of offset for the 2nd law
 * RETURN:            int     current status
 ************************************************************************/
int ibis_isgr_energyReadCal(dal_element *isgrOffsTabPtr,
                         dal_element    *isgrRiseTabPtr,
                         dal_element    *supGainTabPtr,
                         dal_element    *supOffsTabPtr,
                         int             chatter,
                         double        **isgriGoTab,
                         int            *isgriPixTab,
                         short          *isgriRtTab,
                         char         *phGainDOLstr,
                         char         *phOffDOLstr,
                         double         *supCoeffg,
                         double         *supCoeffo)
{
  int    status = ISDC_OK,
         i, j;
  double gt,gh,ot,oh;
  double *tempo_tab=NULL;
  dal_dataType type;
  /*------------ variable for Lut3D ----------------------*/
  long my_numValues,
       my_startValues[3]={1,1,1},
       my_endValues[3]={ISGRI_RT_N_ENER_SCALED,ISGRI_RT_N_DATA,ISGRI_RT_N_RANDOM_DIM};
  /*{1024,256,500};*/
  /*------------------------------------------------*/

  do {

    /*#################################################################*/
    /* ISGRI gain-offset calibration data */
    /*#################################################################*/
    if (chatter > 3) RILlogMessage(NULL, Log_0, "ISGRI gain-offset LUT1 reading...");
    for (j=1; j<ISGRI_GO_N_COL; j++) {
      type=DAL_DOUBLE;
      status=DALtableGetCol(isgrOffsTabPtr, NULL, j, &type, NULL,
                            (void *)(isgriGoTab[j-1]), status);
    }
    type=DAL_INT;
    status=DALtableGetCol(isgrOffsTabPtr, NULL, j, &type, NULL,
                          (void *)(isgriPixTab), status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Cannot read LUT1 columns. Status=%d", status);
      break;
    }
    
    for (i=0; i<ISGRI_N_PIX; i++) {
      /* basic coefficients */
      gh = isgriGoTab[0][i]/10. * 2.;
      oh = isgriGoTab[1][i] * 2.;
      gt = isgriGoTab[2][i]/100.* 30.;
      ot = (isgriGoTab[3][i]+2) * 20.;
      isgriGoTab[0][i]=gh;
      isgriGoTab[1][i]=oh;
      isgriGoTab[2][i]=gt;
      isgriGoTab[3][i]=ot;
    }

    /*#################################################################*/
    /*  ISGRI rise-time calibration data */
    /*#################################################################*/
    if (chatter > 3) RILlogMessage(NULL, Log_0, "ISGRI rise-time LUT2 reading...");

    /*--------- data reading from new 3D LUT2 --------------------------*/
    type=DAL_SHORT;/*before DAL_FLOAT, SPR 4537*/
    status=DALarrayGetSection(isgrRiseTabPtr, ISGRI_DIM_LUT2_3D, my_startValues,
                              my_endValues, &type, &my_numValues,
                              (void *)isgriRtTab, status);
    /* NB:
       isgriRtTab[i][j][k]
       <=>
       isgri_RtTab[i+ISGRI_RT_N_ENER_SCALED*j+ISGRI_RT_N_ENER_SCALED*ISGRI_RT_N_DATA*k];
    */
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "Cannot read LUT2 array. Status=%d", status);
      break;
    }

    /*#################################################################*/
    /*  Calibration data for Method2 (events above 50kev) */
    /*#################################################################*/

    tempo_tab=(double*)calloc(ISGRI_RT_N_DATA,sizeof(double));
    if (tempo_tab == NULL) break;
    
    if (chatter > 3) RILlogMessage(NULL, Log_0, "ISGRI gain and offset corrections reading...");

    type=DAL_DOUBLE;
    for (j=0; j<ISGRI_PHG2_N_COL; j++) {
      status=DALtableGetCol(supGainTabPtr,NULL,j+1,&type,NULL,(void *)tempo_tab,status);
      if (status != ISDC_OK) {
        RILlogMessage(NULL, Error_2, "Cannot read column %d in gain table %13s. Status=%d",
		      j, DS_PHG2, status);
        break;
      }
      for (i=0; i<ISGRI_RT_N_DATA; i++) supCoeffg[i+j*ISGRI_RT_N_DATA]=tempo_tab[i];
    }  
    if (status != ISDC_OK) break;
    
    for (j=0; j<ISGRI_PHO2_N_COL; j++) {
      status=DALtableGetCol(supOffsTabPtr,NULL,j+1,&type,NULL,(void *)tempo_tab,status);
      if (status != ISDC_OK) {
        RILlogMessage(NULL, Error_2, "Cannot read column %d in offset table %13s. Status=%d",
		      j+1, DS_PHO2, status);
        break;
      }
      for (i=0; i<ISGRI_RT_N_DATA; i++) supCoeffo[i+j*ISGRI_RT_N_DATA]=tempo_tab[i];
    }
    
  /*-----------------------------------------------------*/


    if (status != ISDC_OK) break;

  } while(0);

  /*------ free memory -----*/
  if (tempo_tab != NULL) free(tempo_tab);
  return status;
}


/************************************************************************
 * FUNCTION:  ibis_energyIsgrHkCal
 * DESCRIPTION:
 *  Reads the converted HK1 needed for ISGRI LUT1.
 *  Returns ISDC_OK if everything is fine, else returns an error code.
 * ERROR CODES:
 *  DAL error codes
 *  I_ISGR_ERR_MEMORY             Memory allocation error
 *
 * PARAMETERS:
 *  workGRP       dal_element *      in  DOL of the working group
 *  obtStart           OBTime        in  start of time range
 *  obtEnd             OBTime        in  end   of time range
 *  meanT[]              double  in/out  calculated ISGRI mean temperature per MCE
 *  meanBias[]         double    in/out  calculated mean bias per MCE
 *  chatter               int        in  verbosity level
 *  status                int        in  input status
 * RETURN:            int     current status
 ************************************************************************/
int ibis_energyIsgrHkCal(dal_element *workGRP,
                         OBTime       obtStart,
                         OBTime       obtEnd,
                         double       meanT[8],
                         double       meanBias[8],
                         int          chatter,
                         int          status)
{
  int     j, totMCE;
  char    hkName[20], num[3];
  long    i, nValues,
          totVal[8];
  double  myMean, myTot,
    //meanTscw[8]={-51.0,-51.0,-51.0 -51.0, -51.0,-51.0,-51.0,-51.0},
         *hkBuff=NULL;
  OBTime *obtime,
          startTime=DAL3_NO_OBTIME,
          endTime = DAL3_NO_OBTIME;
  dal_dataType dataType;

  if (status != ISDC_OK) return status;
  /* SPR 3686: if OBT limits are valid, use S1 PRP OBT limits */
  if (obtStart != DAL3_NO_OBTIME) {

    if ( (obtime=(OBTime *)calloc(1, sizeof(OBTime))) == NULL)
      return(I_ISGR_ERR_MEMORY);
    endTime=obtEnd;
    /* check if OBT are not too close */
  do {
    status=DAL3GENelapsedOBT(obtEnd, obtStart, &startTime, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Warning_1, "Error calculating elapsed OBT");
      RILlogMessage(NULL, Warning_1, "Reverting from status=%d to ISDC_OK",
                                    status);
      status=ISDC_OK;
      break;
    }
    *obtime=(DAL3_OBT_SECOND) * (SEC_DELTA_MIN);
    if (startTime < *obtime) {

      RILlogMessage(NULL, Warning_1, "Last OBT (%020lld) too close from first OBT (%020lld)",
                                    obtEnd, obtStart);
      status=DAL3GENskipOBT(obtStart, *obtime, &startTime, status);
      if (status != ISDC_OK) {
        RILlogMessage(NULL, Warning_1, "Error adding %d seconds to first OBT", SEC_DELTA_MIN);
        RILlogMessage(NULL, Warning_1, "Reverting from status=%d to ISDC_OK",
                                      status);
        status=ISDC_OK;
        break;
      }
      endTime=startTime;
      RILlogMessage(NULL, Warning_1, "Adding %d seconds to first OBT, last OBT = %020lld",
                                    SEC_DELTA_MIN, endTime);
    }
  } while (0);
    free(obtime);
    startTime=obtStart;

  }
  obtime=NULL;
  totMCE=0;
  for (j=0; j<8; j++) {

    strcpy(hkName, KEY_MCE_TEMP);
    sprintf(num, "%d", j);
    strcat(hkName, num);
    /* get the information for the requested data */
    status=DAL3HKgetValueInfo(workGRP, hkName, IBIS, DAL3HK_CONVERTED,
                              startTime, endTime,  &nValues, &dataType, status);
    if (status != ISDC_OK) break;
    if (nValues < 1) {
      RILlogMessage(NULL, Warning_1, "%13s has NO valid row.", DS_ISGR_HK);
      status=I_ISGR_ERR_BAD_INPUT;
      break;
    }
    /* now allocate the required memory and get the data values */
    status=DALallocateDataBuffer((void **)&hkBuff, nValues*sizeof(double), status);
    status=DALallocateDataBuffer((void **)&obtime, nValues*sizeof(OBTime), status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Cannot allocate buffers for HK data: %s",
                                    hkName);
      break;
    }
    /* retrieve the housekeeping data from the fits file */
    dataType=DAL_DOUBLE;
    status=DAL3HKgetValues(workGRP, hkName, IBIS, DAL3HK_CONVERTED,
                           startTime, endTime, obtime, hkBuff, dataType, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Cannot get buffers for HK data: %s",
                                    hkName);
      break;
    }
    totVal[j]=0;
    myMean=0.0;
    for (i=0; i<nValues; i++) {
      if (hkBuff[i]  >  KEY_MIN_TEMP) { myMean+=hkBuff[i]; totVal[j]++; }
    }
    if (totVal[j] < 1) {
      RILlogMessage(NULL, Warning_1, "Column %19s has NO valid data.", hkName);
      RILlogMessage(NULL, Warning_1, "Continue to next HK data");
    }
    else {
      meanT[j]=myMean/totVal[j];
      totMCE++;
    }

  }
  if (status != ISDC_OK) {
    /* big problem: don't even try to get next HK data */
    DAL3HKfreeData(status);
    if (hkBuff != NULL) DALfreeDataBuffer((void *)hkBuff, ISDC_OK);
    if (obtime != NULL) DALfreeDataBuffer((void *)obtime, ISDC_OK);
    return status;
  }
  if (totMCE == 0) {
    myMean=KEY_DEF_TEMP;
    RILlogMessage(NULL, Warning_1, "Using default ISGRI mean temperature: %+6.2f degC",
                                  myMean);
  }
  else {
    myMean=0.0; myTot=0.0;
    for (j=0; j<8; j++)
      if (totVal[j] > 0) { myMean+=meanT[j]; myTot+=totVal[j]; }
    myMean/=totMCE;
    if (chatter > 1)
      RILlogMessage(NULL, Log_1, "Mean temp. (%05.1f values) on %d MCEs: %+6.2f degC",
                                myTot/totMCE, totMCE, myMean);
    /* Check probe is OK, otherwise re-computes the mean with valid values */
    i=totMCE;
    totMCE=0;
    for (j=0; j<8; j++)
      if (totVal[j] > 0) {
        if (fabs(meanT[j]-myMean-DtempH1[j]) > KEY_RMS_TEMP) {
          RILlogMessage(NULL, Warning_2,
                             "REJECTING mean temp. on MDU%d: %+6.2f degC",
                             j, meanT[j]);
          totVal[j]=0;
        }
        else totMCE++;
      }
    if (i != totMCE) {
      if (totMCE == 0) {
/*        myMean=KEY_DEF_TEMP; keep previous calculation SPR 4838*/
        RILlogMessage(NULL, Warning_2,
                           "NO new mean temp., CHANGE DtempH1 array");
      }
      else {
        myMean=0.0; myTot=0.0;
        for (j=0; j<8; j++)
          if (totVal[j] > 0) { myMean+=meanT[j]; myTot+=totVal[j]; }
        myMean/=totMCE;
        if (chatter > 1)
          RILlogMessage(NULL, Log_1, "NEW  mean  (%05.1f values) on %d MCEs: %+6.2f degC",
                                    myTot/totMCE, totMCE, myMean);
      }
    }
  }
if (chatter > 3) {
  totMCE=0;
  for (j=0; j<8; j++) {

    strcpy(hkName, KEY_MCE_BIAS);
    sprintf(num, "%d", j);
    strcat(hkName, num);
    /* get the information for the requested data */
    status=DAL3HKgetValueInfo(workGRP, hkName, IBIS, DAL3HK_CONVERTED,
                              startTime, endTime,  &nValues, &dataType, status);
    if (status != ISDC_OK) break;
    if (nValues < 1) {
      RILlogMessage(NULL, Warning_1, "%13s has NO valid row.", DS_ISGR_HK);
      status=I_ISGR_ERR_BAD_INPUT;
      break;
    }
    /* now allocate the required memory and get the data values */
    status=DALallocateDataBuffer((void **)&hkBuff, nValues*sizeof(double), status);
    status=DALallocateDataBuffer((void **)&obtime, nValues*sizeof(OBTime), status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Cannot allocate buffers for HK data: %s",
                                    hkName);
      break;
    }
    /* retrieve the housekeeping data from the fits file */
    dataType=DAL_DOUBLE;
    status=DAL3HKgetValues(workGRP, hkName, IBIS, DAL3HK_CONVERTED,
                           startTime, endTime, obtime, hkBuff, dataType, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Cannot get buffers for HK data: %s",
                                    hkName);
      break;
    }
    totVal[j]=0;
    myMean=0.0;
    for (i=0; i<nValues; i++) {
      if (hkBuff[i]  <  KEY_MAX_BIAS) { myMean+=hkBuff[i]; totVal[j]++; }
    }
    if (totVal[j] < 1) {
      RILlogMessage(NULL, Warning_1, "Column %19s has NO valid data.", hkName);
      RILlogMessage(NULL, Warning_1, "Continue to next HK data");
      /* meanBias[j] not changed, already contains default KEY_DEF_BIAS */
    }
    else {
      meanBias[j]=myMean/totVal[j];
      totMCE++;
    }

  }
  if (totMCE) {
    myMean=0.0; myTot=0.0;
    for (j=0; j<8; j++)
      if (totVal[j] > 0) { myMean+=meanBias[j]; myTot+=totVal[j]; }
    myMean/=totMCE;
    if (chatter > 1)
      RILlogMessage(NULL, Log_1, "Mean bias (%05.1f values) on %d MCEs: %+6.1f V",
                                myTot/totMCE, totMCE, myMean);
  }
  else
    RILlogMessage(NULL, Warning_1, "Using default ISGRI mean bias: %+6.1f V",
                                  KEY_DEF_BIAS);
}
  /* in any case must call this function to free internal buffers */
  DAL3HKfreeData(status);
  if (hkBuff != NULL) DALfreeDataBuffer((void *)hkBuff, ISDC_OK);
  if (obtime != NULL) DALfreeDataBuffer((void *)obtime, ISDC_OK);
  return status;
}
/************************************************************************
 * FUNCTION:  ibis_isgr_energyTransform
 * DESCRIPTION:
 *  Computes corrected energies and saves the results.
 *  Returns ISDC_OK if everything is fine, else returns an error code.
 * ERROR CODES:
 *  DAL error codes
 *  I_ISGR_ERR_MEMORY             Memory allocation error
 *
 * PARAMETERS:
 *  outTable  dal_element*    in      DOL of the output table
 *  numEvents        long     in      number of events
 *  revol            long     in      revolution number
 *  isgriGoTab[]   double*    in      ISGRI gain-offset table
 *  isgriPixTab       int*    in      ISGRI pixel type table 
 *  isgriRtTab      short*    in      ISGRI rise-time calibration table (3D)
 *  supCoeffg      double*    in      Coefficients of gain for the 2nd law
 *  supCoeffo      double*    in      Coefficients of offset for the 2nd law
 *  chatter           int     in      verbosity level
 *  meanTemp    double[8]     in      MDU temperature
 *  meanBias    double[8]     in      MDU bias
 *  isgriPha    DAL3_Word*    in      ISGRI raw energy (pulse height)
 *  riseTime    DAL3_Byte*    in      ISGRI rise time
 *  isgriY      DAL3_Byte*    in      ISGRI Y position
 *  isgriZ      DAL3_Byte*    in      ISGRI Z position
 *  isgrPi      DAL3_Byte*  in/out    ISGRI corrected rise-time
 *  isgrEnergy      float*  in/out    ISGRI corrected energy vector (keV)
 * RETURN:            int     current status
 ************************************************************************/
int ibis_isgr_energyTransform   (dal_element *outTable,
				 long         numEvents,
				 long         revol,	 
				 double     **isgriGoTab,
				 int         *isgriPixTab,
				 short       *isgriRtTab,
				 double      *supCoeffg,
				 double      *supCoeffo,
				 int          chatter,
				 double       meanTemp[8],
				 double       meanBias[8],
				 DAL3_Word   *isgriPha,
				 DAL3_Byte   *riseTime,
				 DAL3_Byte   *isgriY,
				 DAL3_Byte   *isgriZ,
				 DAL3_Byte   *isgrPi,
				 DAL3_Word   *isgrPHA2,
				 DAL3_Byte   *isgrRT1,
				 DAL3_Word   *isgrPHA1,
				 float       *isgrEnergy,
                 int          corGainDrift)
{
  int    status = ISDC_OK,
    pixelNo,mce,
         irt, ipha,ipha2;
  long   i,j,              /* loop index */
         index_cc,
         infoEvt[7]={0, 0, 0, 0, 0, -99, 0};
  double index_cc_real,deltin,delte;
  double gt[ISGRI_N_PIX],
         gh[ISGRI_N_PIX],
         ot[ISGRI_N_PIX],
         oh[ISGRI_N_PIX] ;  /*lut1 coefficients*/
  double gain_corrPH1,offset_corrPH1, /*correction coefficient for the first law*/
         rt, pha,
         ener, corr;
  double random_num;

  double slopeMDU[8]={-1.8,-2.0,-2.3,-2.7,-0.5,-2.4,-0.8,-0.5} ;
  
  /*------------ variables for new method -----------*/
  /*correction coefficients for the second law */
  double *gain_corrPH2=NULL,
         *offset_corrPH2=NULL,
	  off_scale=0.,g_scale=1.,
          Ec=50.,*Chc=NULL;/*Threshold in keV and in channel to 
                             apply both calibration method*/
do {
  status=I_ISGR_ERR_MEMORY;
  gain_corrPH2=(double*)calloc(ISGRI_RT_N_DATA ,sizeof(double));
  if (gain_corrPH2 == NULL) break;
  offset_corrPH2=(double*)calloc(ISGRI_RT_N_DATA ,sizeof(double));
  if (offset_corrPH2 == NULL) break;

  Chc=(double*)calloc(ISGRI_RT_N_DATA ,sizeof(double));
  if (Chc == NULL) break;
  status=ISDC_OK;

  for (j=0;j<8;j++) 
    {
      meanTemp[j]=(meanTemp[j]+273.0)/273.0; /* to scale temperature in ratio to minimum Kelvin */
      meanBias[j]=-meanBias[j]/100. ;
      meanBias[j]=1.2 ;
    }

  for (i=0;i<128;i++)
    for (j=0;j<128;j++)
      {
	pixelNo = 128*i+j;
	mce     = 7 - j/32 - 4*(i/64);
	gh[pixelNo] = isgriGoTab[0][pixelNo]*pow(meanTemp[mce],-1.11)*pow(meanBias[mce],-0.0832);
	oh[pixelNo] = isgriGoTab[1][pixelNo]*pow(meanTemp[mce],slopeMDU[mce])*pow(meanBias[mce],0.0288);
	gt[pixelNo] = isgriGoTab[2][pixelNo]*pow(meanTemp[mce],0.518)*pow(meanBias[mce],0.583);
	ot[pixelNo] = isgriGoTab[3][pixelNo]*pow(meanTemp[mce],0.625)*pow(meanBias[mce],0.530);
      }
  
  /*Law 1 correction: offset = law 2(revol=0,RT=0) fixed, gain calibration with W line
  par[12]_[gain,offset]_corrPH1 must also be defined in ii_shadow_build_types.h*/
  gain_corrPH1   = PAR1_GAIN_corrPH1+PAR2_GAIN_corrPH1*revol;
  offset_corrPH1 = PAR1_OFFSET_corrPH1;
  
  g_scale = G_SCALE0+G_SCALE1*revol;
  off_scale = OFF_SCALE0;

  for (irt=0; irt<256; irt++) {
    gain_corrPH2[irt]   = supCoeffg[irt+0*ISGRI_RT_N_DATA] + supCoeffg[irt+1*ISGRI_RT_N_DATA]*revol;
    offset_corrPH2[irt] = supCoeffo[irt+0*ISGRI_RT_N_DATA] + supCoeffo[irt+1*ISGRI_RT_N_DATA]*revol + supCoeffo[irt+2*ISGRI_RT_N_DATA]*revol*revol;
    Chc[irt]=0.; 
    if (gain_corrPH1!=gain_corrPH2[irt]) {
      Ec=(offset_corrPH2[irt]-offset_corrPH1)/(gain_corrPH1-gain_corrPH2[irt]);
      Chc[irt]=Ec*gain_corrPH1+offset_corrPH1;
    }
  }

  if (chatter > 2) {
    RILlogMessage(NULL, Log_0, "GAIN1, OFFSET1 for PHA        : %8.6f  %8.4f",
                              gain_corrPH1, offset_corrPH1);
    RILlogMessage(NULL, Log_0, "GAIN2, OFFSET2 for PHA (RT=35): %8.6f  %8.4f",
                              gain_corrPH2[35], offset_corrPH2[35]);
    RILlogMessage(NULL, Log_0, "Equal energy for laws  (RT=35): %5.2f keV",
                              (Chc[35]-offset_corrPH1)/gain_corrPH1 );
  }
 
  /*#################################################################*/
  /* Compute corrected energies  */
  /*#################################################################*/
    for (j=0L; j < numEvents; j++)  {
      
      /*--- LUT1 corrections ---*/
      pixelNo = 128*(int)isgriY[j] + (int)isgriZ[j];

    rt = 2.*riseTime[j]/2.4+5.0;  /* 256 channels for LUT2 calibration scaled */
    rt = rt*gt[pixelNo] + ot[pixelNo];

    pha = isgriPha[j]*gh[pixelNo] + oh[pixelNo];

    isgrRT1[j]=(DAL3_Byte)rt;
    isgrPHA1[j]=(DAL3_Word)pha;

    /*--- drift correction ---*/
    /*rt=(rt-offset_corrRT)/gain_corrRT; Rt effect included into PH gain2 and offset2*/
    corr=rt*1000.;
    if (corr > infoEvt[5]) infoEvt[5]=corr;
    else if (corr < infoEvt[4]) infoEvt[4]=corr;

    if ((rt-floor(rt)) <0.5) irt = (long)floor(rt) ;
        else irt = (long)ceil(rt);
    /* if rise-time out of range take LUT2 limit, SPR 4549 */
    if (irt < 0)        {irt=0;   infoEvt[0]++;}
    else if (irt > 255) {irt=255; infoEvt[1]++;}

    if (corGainDrift) {
        random_num=DAL3GENrandomDoubleX1()-0.5;
        if (pha < Chc[irt]) pha=(pha+random_num-offset_corrPH1)/gain_corrPH1;
        else           pha=(pha+random_num-offset_corrPH2[irt])/gain_corrPH2[irt];
        pha=  2*(pha-off_scale)/g_scale;
    } else {
        pha=  (pha-off_scale)/g_scale;
    };

    ipha= (int)pha;

    /*random_num=500.0*rand()/(1+(double)RAND_MAX);*/
    random_num=500.0*DAL3GENrandomDoubleX1();


   if ((ipha >= 0) && (ipha < 2048)) {
      /*ener = (double) 2.0*isgriRtTab[(int)(pha/2.0)][irt][(int)random_num]/30.0;*/
      if ((pha/2.0-floor(pha/2.0)) <0.5) ipha2 = (long)floor(pha/2.0) ;
        else ipha2 = (long)ceil(pha/2.0);
      index_cc= ipha2
	+ISGRI_RT_N_ENER_SCALED*irt
                    +ISGRI_RT_N_ENER_SCALED*ISGRI_RT_N_DATA*(long)random_num;
      /*NP Piotr start*/
      index_cc_real=pha/2.0
	+ISGRI_RT_N_ENER_SCALED*(long)irt
	+ISGRI_RT_N_ENER_SCALED*ISGRI_RT_N_DATA*(long)random_num;
      deltin = index_cc_real-(double)index_cc;
      if (index_cc+1<ISGRI_RT_N_ENER_SCALED*ISGRI_RT_N_DATA*ISGRI_RT_N_RANDOM_DIM){
        delte = deltin*((double)isgriRtTab[index_cc+1]-(double)isgriRtTab[index_cc]);
      }
      else{
        delte=0;
      }
      ener=FINALOFFS+((double)isgriRtTab[index_cc]+delte)/30.0;
      /*NP Piotr end*/
    }
    else if (ipha >= 2048) {
      /*ener =(double) 2.0*isgriRtTab[(int)(ISGRI_RT_N_ENER/2.0-1)][irt][(int)random_num]/30.0;*/
      index_cc= ISGRI_RT_N_ENER_SCALED-1
               +ISGRI_RT_N_ENER_SCALED*(long)irt
               +ISGRI_RT_N_ENER_SCALED*ISGRI_RT_N_DATA*(long)random_num;
      ener=FINALOFFS+(double)isgriRtTab[index_cc]/30.0;
      infoEvt[3]++;
    }
    else {ener=0.0; infoEvt[2]++;}

    if (ener <= FINALOFFS) { isgrEnergy[j] = 0.0;  /* if, by mistake, isgriRtTab<0*/
                             infoEvt[6]++;}        /* replace 0 by 1.3664, SPR 4664 */
    else isgrEnergy[j] = (float)ener;

    isgrPi[j]=(DAL3_Byte)irt;
    isgrPHA2[j]=(DAL3_Word)ipha2;

  }
  if (chatter > 1) {
    RILlogMessage(NULL, Log_1, "Total COR rise-time <= -1: %9ld", infoEvt[0]);
    RILlogMessage(NULL, Log_1, "Total COR rise-time >=256: %9ld", infoEvt[1]);
    ener=(double)infoEvt[4]/1000.;
    corr=(double)infoEvt[5]/1000.;
    RILlogMessage(NULL, Log_1, "COR rise-time interval: %4.1f to %5.1f", ener,corr);
    RILlogMessage(NULL, Log_1, "Total COR amplitude <= -1: %9ld", infoEvt[2]);
    RILlogMessage(NULL, Log_1, "Total COR amplitude>=2048: %9ld", infoEvt[3]);
    RILlogMessage(NULL, Log_1, "Total with LUT2 coef. <=0: %9ld",
                              infoEvt[6]-infoEvt[2]);
  }

  /*#################################################################*/
  /* Put computed columns into the output table */
  /*#################################################################*/
  status=DALtablePutCol(outTable, "ISGRI_PHA2", 0, DAL_USHORT, numEvents,
                        isgrPHA2, status);
  status=DALtablePutCol(outTable, "ISGRI_RT1", 0, DAL_BYTE, numEvents,
                        isgrRT1, status);
  status=DALtablePutCol(outTable, "ISGRI_PHA1", 0, DAL_USHORT, numEvents,
                        isgrPHA1, status);
  status=DALtablePutCol(outTable, KEY_COL_OUT, 0, DAL_BYTE, numEvents,
                        isgrPi, status);
  status=DALtablePutCol(outTable, "ISGRI_ENERGY", 0, DAL_FLOAT, numEvents,
                        isgrEnergy, status);
  if (status != ISDC_OK)
    RILlogMessage(NULL, Error_2, "Cannot write output data. Status=%d", status);
  status=CommonStampObject(outTable, "Energy correction.", status);

} while (0);
  /*--------- free memory ---------------------*/
  if (gain_corrPH2 !=NULL) free(gain_corrPH2);
  if (offset_corrPH2 !=NULL) free(offset_corrPH2);
  if (Chc != NULL) free(Chc);
  return status;
}
