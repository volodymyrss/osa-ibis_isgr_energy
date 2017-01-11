/************************************************************************
 * FILE:        ibis_isgr_energy_main.c
 * VERSION:     8.3
 * COMPONENT:   ibis_isgr_energy
 * AUTHOR:      C. Couvreur &
 *              A. Sauvageon (SAp-CEA,  asauvageon@cea.fr)
 * DESCRIPTION: main program of the component
 * HISTORY:
 *   MD,     28/03/1999, template version
 *   MD,     09/04/1999, demi-module version
 *   PL,     20/05/1999, first module version
 *   PL,     09/06/1999, add lecture of the gain-offset file, 
 *                       add a column for the corrected energy in the cor 
 *   PL,     09/06/1999, output FITS file
 *   PL,     17/11/1999, add COMMON_PREPARE_PARS and COMMON_CLOSE_SWG
 *   PL,     07/03/2000, add STAMPS
 *   PL,     24/04/2002, SPR 1052 1259 1369 (correct reading of IC files)
 *   PL,3.5,   /09/2002, SPR 1413 1511 1636 (filling header attributes)
 *   PL,4.0, 30/07/2003, SCREW 1229 (new LUT2)
 *  ASA,4.1, 25/08/2003, SPR 3137 (bad formula to use LUT2)
 *   PL,4.2, 07/10/2003, SCREW 1283 (read bias and temperature in HK)
 *   PL,4.3, 17/10/2003, SPR 3257 (mean HK correctly calculated)
 *
 *   code into C, ASA, december 2003
 *
 *  ASA,5.0, 30/03/2004, SPR 3418 (partly) SCREW 1398 (new parameters)
 *                       SCREW 1384 (column ISGRI_RT)
 *  ASA,5.1, 08/04/2004, SPR 3418 (header) and new PHA gain laws
 *                       SCREW 1430 (column back to ISGRI_PI)
 *  ASA,5.2, 25/05/2004, SPR 3560, 3561
 *  ASA,5.3, 16/06/2004, SPR 3686
 *
 *  ASA,5.4, 31/03/2005, SCREW 1499 (read parameters from IC file)
 *  CC, 5.6, 28/07/2006, new version (new parameters)
 *  CC, 5.7
 *  ASA,5.7, 14/12/2006, SPR 4641 + cosmetic changes
 *  CC,5.7.1,05/03/2007, SPR 4664 (correction of energy drift around 23kev, 
 *                       new input parameter)
 *  CC, 6.0, 11/07/2007, SPR 4695 new parameters: supGDOL and supODOL
 *                       parameter switDOL and rtcDOL deleted.
 * ASA, 6.2, 20/12/2007, SPR 4773 4698 (+ cosmetic changes)
 *  PL, 8.0  02/02/2012,   remove IREM counters, Temperature correction by MDU
 *  PL, 8.2  02/04/2012,   modify coefficients PAR1_..._corrPH1 to correct <50 keV behavior
 ************************************************************************/

#include "ibis_isgr_energy.h" 


int main (int argc, char *argv[])
{
  int   status = ISDC_OK,
        makeUnique = 1,
        clobber,
        gti,
        erase, chatter;

  char *riseDOLstr = NULL,
       *acorDOLstr = NULL,
   /*  *switDOLstr = NULL,*/
   /*  *rteffectDOLstr = NULL,*/
       *phGainDOLstr=NULL,
    *phOffsDOLstr=NULL ;
  char  randName[DAL_BIG_STRING];
  unsigned long seed;
  dal_element  *workGRP = NULL;

  /* initialize the common library stuff */
  status=CommonInit(COMPONENT_NAME, COMPONENT_VERSION, argc, argv);
  if (status != ISDC_SINGLE_MODE) {
    RILlogMessage(NULL, Warning_2, "CommonInit status = %d", status);
    RILlogMessage(NULL, Warning_2, "number of command line arguments = %d", argc);
    RILlogMessage(NULL, Warning_2, "program name : %s", argv[0]);
    RILlogMessage(NULL, Warning_2, "Program aborted : could not initialize.");
    CommonExit(status);
  }
  do {
    status=I_ISGR_ERR_MEMORY;
    if ( (riseDOLstr=(char *)calloc(DAL_FILE_NAME_STRING, sizeof(char))) == NULL)
      break;
    if ( (acorDOLstr=(char *)calloc(DAL_FILE_NAME_STRING, sizeof(char))) == NULL)
      break;
    /*if ( (switDOLstr=(char *)calloc(DAL_FILE_NAME_STRING, sizeof(char))) == NULL)
      break;*/
    /*if ( (rteffectDOLstr=(char *)calloc(DAL_FILE_NAME_STRING, sizeof(char))) == NULL)
      break;*/
    if ( (phGainDOLstr=(char *)calloc(DAL_FILE_NAME_STRING, sizeof(char))) == NULL)
      break;
    if ( (phOffsDOLstr=(char *)calloc(DAL_FILE_NAME_STRING, sizeof(char))) == NULL)
      break;
    status=ISDC_OK;
  /*#################################################################*/
  /* get the parameters */
  /*#################################################################*/
    status=PILGetInt("chatter", &chatter);
    if (status != ISDC_OK) break;
    RILlogMessage(NULL, Log_2, "Verbosity level = %d", chatter);

    status=PILGetBool("useGTI", &gti);
    if (status != ISDC_OK) break;
    if (chatter > 0) {
      if (gti)
        RILlogMessage(NULL, Log_2,"Number of GTI: 1 (PRP must contain OBTs)");
      else
        RILlogMessage(NULL, Log_2,"Number of GTI: 0 (PRP can have 0 OBT)");
    }
    status=PILGetBool("eraseALL", &erase);
    if (status != ISDC_OK) break;
    if (chatter > 0) {
      if (erase)
        RILlogMessage(NULL, Log_2,"Erase output rows");
      else
        RILlogMessage(NULL, Log_2,"Replace output columns");
    }
    /*SPR 4664 ----*/
    /*status=PILGetInt("RTdriftCor", &RTdriftCor);
    if (status != ISDC_OK) break;
    switch (RTdriftCor){
    case 0: 
      RILlogMessage(NULL, Log_2,"RT corrections are not taken into account.");
      break;
    case 1: 
      RILlogMessage(NULL, Log_2,"RT gain and offset are taken into account.");
      break;
    case 2: 
      RILlogMessage(NULL, Log_2,"All RT corrections are taken into account.");
      break;
    default : 
      status=I_ISGR_ERR_BAD_INPUT;
      RILlogMessage(NULL, Error_2,"Bad Value for the parameter RTdriftCor");
    }
    if (status != ISDC_OK) break;*/
    /*--------*/

    status=PILGetString("randSeed", randName);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "The parameter 'randSeed' is not found.");
      break;
    }
    if (strlen(randName) > 0) {
      seed=strtoul(randName, (char **)NULL, 10);
      if (seed == ULONG_MAX) {
        RILlogMessage(NULL, Error_2, "Parameter 'randSeed' interval: 0<=  <%lu",
                                    ULONG_MAX);
        status=I_ISGR_ERR_BAD_INPUT;
        break;
      }
      RILlogMessage(NULL, Log_2, "Seed for random number generator: %010lu", seed);
      DAL3GENrandomSeed(seed);
    }

    status=PILGetString("riseDOL", riseDOLstr);
    if (status != ISDC_OK) break;
    if (strlen(riseDOLstr) == 0) {
      RILlogMessage(NULL, Error_2, "The parameter 'riseDOL' is empty");
      status=I_ISGR_ERR_BAD_INPUT;
      break;
    }
    status=PILGetString("GODOL", acorDOLstr);
    if (status != ISDC_OK) break;
    if (strlen(acorDOLstr) == 0) {
      RILlogMessage(NULL, Error_2, "The parameter 'GODOL' is empty");
      status=I_ISGR_ERR_BAD_INPUT;
      break;
    }

    status=PILGetString("supGDOL", phGainDOLstr);
    if (status != ISDC_OK) break;
    if (strlen(phGainDOLstr) == 0) {
      RILlogMessage(NULL, Error_2, "The parameter 'supGDOL' is empty");
      status=I_ISGR_ERR_BAD_INPUT;
      break;
    }
    status=PILGetString("supODOL", phOffsDOLstr);
    if (status != ISDC_OK) break;
    if (strlen(phOffsDOLstr) == 0) {
      RILlogMessage(NULL, Error_2, "The parameter 'supODOL' is empty");
      status=I_ISGR_ERR_BAD_INPUT;
      break;
    } /* hkCnvDOL was erased (so useless !) since 5.6.0 */
    status=CommonPreparePARsStrings("inGRP",
                                    "inRawEvts,hkCnvDOL",
                                    "outGRP",
                                    "outCorEvts",
                                    makeUnique,
                                    &workGRP,
                                    &clobber,
                                    status);
    if (status == GROUP_NOT_FOUND) {
      RILlogMessage(NULL, Error_2, "Program aborted: GROUP not found.");
      break;
    }
    else if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_2, "CommonPreparePARsStrings error. Status=%d",
                                  status);
      break;
    }

  /*#################################################################*/
  /* do the work */
  /*#################################################################*/
    status=ibis_isgr_energyWork(workGRP, gti, erase, chatter,
                                acorDOLstr, riseDOLstr,
                                phGainDOLstr, phOffsDOLstr);
    if (status == I_ISGR_ERR_MEMORY) {
      RILlogMessage(NULL, Error_2, "Program aborted: memory allocation error.");
    }
    if (workGRP != NULL) status=CommonCloseSWG(workGRP, status);
    /* if error while closing, status is only changed if it was ISDC_OK */

  } while(0);

  if (riseDOLstr != NULL) free(riseDOLstr);
  if (acorDOLstr != NULL) free(acorDOLstr);
  /*if(switDOLstr !=NULL) free(switDOLstr);*/
  /*if(rteffectDOLstr !=NULL) free(rteffectDOLstr);*/
  if (phGainDOLstr != NULL) free(phGainDOLstr);
  if (phOffsDOLstr != NULL) free(phOffsDOLstr);

  CommonExit(status);
  return(status);   /* to make lint happy */
}
