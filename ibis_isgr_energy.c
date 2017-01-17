/************************************************************************
 * FILE:        ibis_isgr_energy.c
 * VERSION:     8.3
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
 *  CF, 8.3  05/11/2015,   SCREW 2624 path for low energy correction implemented with Paris agreement
 *  VS, 9.0  17/01/2017,   substrantially rewritten, a lot of functionality moved to dal3ibis_calib
 ******************************************************************************/

#include "ibis_isgr_energy.h"

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
 * RETURN:            int     current status
 ************************************************************************/

int ibis_isgr_energyWork(dal_element *workGRP,
                         ibis_isgr_energy_settings *ptr_ibis_isgr_energy_settings,
                         ibis_isgr_energy_settings *ptr_ISGRI_energy_caldb_dols,
                         int chatter,
                         int status)
{
  int    i,
         freeStatus= ISDC_OK;

  long   numEvents = 0l,      /* number of S1 events   */
         buffSize,            /* size of output buffers in bytes */
         revol;               /* to calculate the PHA gain decrease */

  char  *logString=NULL,

  /* input data */
  DAL3_Word *isgriPha = NULL; /* ISGRI_PHA */
  DAL3_Byte *riseTime = NULL; /* RISE_TIME */
  DAL3_Byte *isgriY   = NULL; /* ISGRI_Y   */
  DAL3_Byte *isgriZ   = NULL; /* ISGRI_Z   */

  /* output data */
  DAL3_Byte *isgrPi = NULL;  
  float *isgrEnergy = NULL;   /* ISGRI_ENERGY  */

  OBTime obtStart=DAL3_NO_OBTIME,
         obtEnd=DAL3_NO_OBTIME;

  ISGRI_energy_calibration_struct ISGRI_energy_calibration;

  dal_element *in_table = NULL,
              *outTable = NULL;

  dal_element *isgrOffsTabPtr = NULL;
  dal_element *isgrRiseTabPtr = NULL;
  dal_element *isgrHK1_Ptr = NULL;
  dal_element *supGainTabPtr=NULL;
  dal_element *supOffsTabPtr=NULL;

  do {

    if ((logString=(char *)calloc(DAL_BIG_STRING, sizeof(char))) == NULL) {
      status=I_ISGR_ERR_MEMORY;
      break;
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

    statis=DAL3IBISupdateTBiasMDUCorrection(status);


    /*#################################################################*/
    /* Perform energy correction */
    /*#################################################################*/

    

    status=DAL3IBISReconstructISGRIEnergy(workGRP,status);

    status=ibis_isgr_energyTransform(outTable, numEvents, revol,
                                 goTab, pixTab,isgriRtTab,
                                 supCoeffg, supCoeffo,
                                 chatter,
                                 meanT,
	                			 meanBias,    
                                 isgriPha, riseTime,
                                 isgriY, isgriZ,
                                 isgrPi, isgrEnergy);
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


int DAL3IBIS_read_ISGRI_events(dal_element *workGRP,
                                 int gti,
                                 ISGRI_events_struct *ISGRI_events,
                                 int chatter,
                                 int status)

{
  short  selected=0;
  long   buffSize;

  dal_dataType type;
  ISDCLevel    myLevel;

  do {

    selected=1;
    if (gti) myLevel=PRP;  else myLevel=RAW;
    status=DAL3IBISselectEvents(workGRP, ISGRI_EVTS, myLevel, gti,
                                obtStart, obtEnd, NULL, status);
    status=DAL3IBISgetNumEvents(numEvents, status);

    if ((status == DAL3IBIS_NO_IBIS_EVENTS) || (status == DAL_TABLE_HAS_NO_ROWS)) {
      RILlogMessage(NULL, Warning_1, "Reverting from status=%d to ISDC_OK", status);
      RILlogMessage(NULL, Warning_1, "NO event selected. Execution stopped.");
    //  status=ISDC_OK;
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
				 float       *isgrEnergy)
{
    int    status = ISDC_OK,
    long   j,              /* loop index */
           infoEvt[7]={0, 0, 0, 0, 0, -99, 0};

    do {
        for (j=0L; j < numEvents; j++)  {

            DAL3IBISReconstructISGRIEnergy(isgriPha,
                    riseTime,
                    isgriY,
                    isgriZ,
                    ISGRI_energy_correction,
                    &isgri_energy,
                    &isgri_pi,
                    &infoEvt,
                    int status);

        };

        if (chatter > 1) {
            RILlogMessage(NULL, Log_1, "Total COR rise-time <= -1: %9ld", infoEvt[0]);
            RILlogMessage(NULL, Log_1, "Total COR rise-time >=256: %9ld", infoEvt[1]);
          //  ener=(double)infoEvt[4]/1000.;
         //   corr=(double)infoEvt[5]/1000.;
         //   RILlogMessage(NULL, Log_1, "COR rise-time interval: %4.1f to %5.1f", ener,corr);
            RILlogMessage(NULL, Log_1, "Total COR amplitude <= -1: %9ld", infoEvt[2]);
            RILlogMessage(NULL, Log_1, "Total COR amplitude>=2048: %9ld", infoEvt[3]);
            RILlogMessage(NULL, Log_1, "Total with LUT2 coef. <=0: %9ld",
                    infoEvt[6]-infoEvt[2]);
        }

        /*#################################################################*/
        /* Put computed columns into the output table */
        /*#################################################################*/
        status=DALtablePutCol(outTable, "ISGRI_PI", 0, DAL_BYTE, numEvents,
                isgrPi, status);
        status=DALtablePutCol(outTable, "ISGRI_ENERGY", 0, DAL_FLOAT, numEvents,
                isgrEnergy, status);
        if (status != ISDC_OK)
            RILlogMessage(NULL, Error_2, "Cannot write output data. Status=%d", status);
        status=CommonStampObject(outTable, "Energy correction.", status);

    } while (0);
    return status;
}
