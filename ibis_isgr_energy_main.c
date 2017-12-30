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
 *  VS, 9.0. 27/01/2017, most of the functionality moved to DAL3IBIS
 *  see also https://github.com/volodymyrss/osa-ibis_isgr_energy
 ************************************************************************/

#include "ibis_isgr_energy.h" 


/************************************************************************
 * FUNCTION:  get_all_PIL
 * DESCRIPTION:
 *  gets the settings from PIL
 * ERROR CODES:
 *  DAL error codes
 *  I_ISGR_ERR_ISGR_OUT_COR   if number of output rows inconsistent
 *  I_ISGR_ERR_BAD_INPUT      bad input somehow
 *  see other parameter by their names
 ************************************************************************/

// returns status
int get_all_PIL(dal_element **ptr_workGRP,
                ibis_isgr_energy_settings_struct *ptr_ibis_isgr_energy_settings,
                ISGRI_energy_caldb_dols_struct *ptr_ISGRI_energy_caldb_dols,
                int *ptr_chatter,
                int status) {
    char  randName[DAL_BIG_STRING];
    unsigned long seed;
    int chatter;

    TRY_BLOCK_BEGIN

        TRY( PILGetInt("chatter", ptr_chatter), status, "reading chatter parameter" );
        chatter=*ptr_chatter;

        RILlogMessage(NULL, Log_2, "Verbosity level = %d", chatter);

        TRY( PILGetBool("useGTI", &ptr_ibis_isgr_energy_settings->gti), status, "unable to read gti parameter" );

        if (chatter > 0) {
            if (ptr_ibis_isgr_energy_settings->gti)
                RILlogMessage(NULL, Log_2,"Number of GTI: 1 (PRP must contain OBTs)");
            else
                RILlogMessage(NULL, Log_2,"Number of GTI: 0 (PRP can have 0 OBT)");
        }

        TRY( PILGetBool("eraseALL", &ptr_ibis_isgr_energy_settings->erase), status, "unable to read eraseALL parameter" );

        if (chatter > 0) {
            if (ptr_ibis_isgr_energy_settings->erase)
                RILlogMessage(NULL, Log_2,"Erase output rows");
            else
                RILlogMessage(NULL, Log_2,"Replace output columns");
        }

        TRY( PILGetString("randSeed", randName), status, "The parameter 'randSeed' is not found.");
        if (strlen(randName) > 0) {
            seed=strtoul(randName, (char **)NULL, 10);
            if (seed == ULONG_MAX) {
                FAIL(I_ISGR_ERR_BAD_INPUT, "Parameter 'randSeed' interval: 0<=  <%lu",ULONG_MAX);
            }
            RILlogMessage(NULL, Log_2, "Seed for random number generator: %010lu", seed);
            DAL3GENrandomSeed(seed);
        }
        
        TRY( PILGetString("GODOL", ptr_ISGRI_energy_caldb_dols->lut1_DOL), status, "reading GODOL parameter" );
        if (strlen(ptr_ISGRI_energy_caldb_dols->lut1_DOL) == 0) FAIL(I_ISGR_ERR_BAD_INPUT,"The parameter 'GODOL' is empty");

        TRY( PILGetString("riseDOL", ptr_ISGRI_energy_caldb_dols->lut2_DOL), status, "reading riseDOL parameter");
        if (strlen(ptr_ISGRI_energy_caldb_dols->lut2_DOL) == 0) FAIL(I_ISGR_ERR_BAD_INPUT,"The parameter 'riseDOL' is empty");

        TRY( PILGetString("mcecDOL", ptr_ISGRI_energy_caldb_dols->mcec_DOL), status, "reading mcecDOL parameter");
        if (strlen(ptr_ISGRI_energy_caldb_dols->mcec_DOL) == 0) FAIL(I_ISGR_ERR_BAD_INPUT,"The parameter 'mcecDOL' is empty");
        
        TRY( PILGetString("l2reDOL", ptr_ISGRI_energy_caldb_dols->l2re_DOL), status, "reading mcecDOL parameter");
        if (strlen(ptr_ISGRI_energy_caldb_dols->l2re_DOL) == 0) FAIL(I_ISGR_ERR_BAD_INPUT,"The parameter 'l2reDOL' is empty");

        TRY( CommonPreparePARsStrings("inGRP",
                "inRawEvts,hkCnvDOL",
                "outGRP",
                "outCorEvts",
                ptr_ibis_isgr_energy_settings->makeUnique,
                ptr_workGRP,
                &ptr_ibis_isgr_energy_settings->clobber,
                status),
                    status,"CommonPreparePARsStrings"
                );

    TRY_BLOCK_END
    return status;
}
  

int main (int argc, char *argv[])
{
  int  status = ISDC_OK;
  int chatter;

  ibis_isgr_energy_settings_struct ibis_isgr_energy_settings;
  ISGRI_energy_caldb_dols_struct ISGRI_energy_caldb_dols;
  dal_element  *workGRP = NULL;
  
  ibis_isgr_energy_settings.makeUnique = 1;

  TRY_BLOCK_BEGIN

      status = CommonInit(COMPONENT_NAME, COMPONENT_VERSION, argc, argv);
      if (status != ISDC_SINGLE_MODE) {
        RILlogMessage(NULL, Warning_2, "CommonInit status = %d", status);
        RILlogMessage(NULL, Warning_2, "number of command line arguments = %d", argc);
        RILlogMessage(NULL, Warning_2, "program name : %s", argv[0]);
        RILlogMessage(NULL, Warning_2, "Program aborted : could not initialize.");
        CommonExit(status);
      }

      TRY( get_all_PIL(&workGRP,&ibis_isgr_energy_settings,&ISGRI_energy_caldb_dols,&chatter,status), status, "get_all_PIL" );

      TRY( ibis_isgr_energyWork(workGRP, &ibis_isgr_energy_settings, &ISGRI_energy_caldb_dols,chatter,status), status, "ibis_isgr_energyWork" );

      if (workGRP != NULL)
          TRY( CommonCloseSWG(workGRP, status), status, "CommonCloseSWG");

      CommonExit(status);

  TRY_BLOCK_END

  DAL_GC_free_all(chatter,status); //GC status is discarded

  return(status); 
}
