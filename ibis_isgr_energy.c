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
                         ibis_isgr_energy_settings_struct *ptr_ibis_isgr_energy_settings,
                         ISGRI_energy_caldb_dols_struct *ptr_ISGRI_energy_caldb_dols,
                         int chatter,
                         int status)
{
    int    i,
           freeStatus= ISDC_OK;
    char  logString[DAL_BIG_STRING];

    IBIS_events_struct IBIS_events;
    ISGRI_energy_calibration_struct ISGRI_energy_calibration;

    TRY_BLOCK_BEGIN

        TRY( DAL3IBIS_print_all_events(workGRP,status), status, "showing events" );

        TRY( DAL3IBIS_read_IBIS_events(workGRP,ISGRI_EVTS,&IBIS_events,ptr_ibis_isgr_energy_settings->gti,chatter,status), -1, "reading events" );

        TRY( DAL3IBIS_init_ISGRI_energy_calibration(&ISGRI_energy_calibration,status), status, "initializing ISGRI energy calibration");

        TRY( DAL3IBIS_populate_newest_DS(&IBIS_events, &ISGRI_energy_calibration, DS_ISGR_LUT1, &DAL3IBIS_open_LUT1, &DAL3IBIS_read_LUT1,chatter,status), status, "reading LUT1" );
        TRY( DAL3IBIS_correct_LUT1_for_temperature_bias(workGRP,&ISGRI_energy_calibration,&IBIS_events,chatter,status), status, "correcting for LUT1 temperature bias");

        TRY( DAL3IBIS_populate_newest_DS(&IBIS_events, &ISGRI_energy_calibration, DS_ISGR_MCEC, &DAL3IBIS_open_MCEC, &DAL3IBIS_read_MCEC, chatter,status), status, "loading MCE evolution correction");
        TRY( DAL3IBIS_populate_newest_DS(&IBIS_events, &ISGRI_energy_calibration, DS_ISGR_LUT2, &DAL3IBIS_open_LUT2, &DAL3IBIS_read_LUT2, chatter,status), status, "loading LUT2" );
        TRY( DAL3IBIS_populate_newest_DS(&IBIS_events, &ISGRI_energy_calibration, DS_ISGR_L2RE, &DAL3IBIS_open_L2RE, &DAL3IBIS_read_L2RE, chatter,status), status, "loading LUT2 rapid evolution" );

    TRY_BLOCK_END

    status=DAL3IBIS_reconstruct_ISGRI_energies(&ISGRI_energy_calibration,&IBIS_events,chatter,status);

    //status=DAL3IBIS_dealocate(&ISGRI_energy_calibration,status);

    ibis_isgr_energyCheckOut(&IBIS_events,workGRP,"ISGR-EVTS-COR",ptr_ibis_isgr_energy_settings,chatter,status);

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
 * RETURN:            int     current status
 ************************************************************************/
int ibis_isgr_energyCheckOut(
        IBIS_events_struct *ptr_IBIS_events,
        dal_element *workGRP,
        char         *outName,
        ibis_isgr_energy_settings_struct *ptr_ibis_isgr_energy_settings,
        int           chatter,
        int           status
        )
{
    long  outRow;
    dal_element *c_outTable;
    dal_element **outTable = &c_outTable;


    do {
        status=DALobjectFindElement(workGRP, outName, outTable, status);
        status=DALtableGetNumRows(*outTable, &outRow, status);
        if (status != ISDC_OK) {
            RILlogMessage(NULL, Error_2, "%13s bintable NOT found.", outName);
            break;
        }
        if (chatter > 2)
            RILlogMessage(NULL, Log_0, "%13s output table: %9ld rows.", outName, outRow);

        if (outRow > 0l) {

            if (ptr_ibis_isgr_energy_settings->erase) {
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
            else if (outRow != ptr_IBIS_events->numEvents) {
                RILlogMessage(NULL, Error_2, "%13s has wrong length (%ld rows).",
                        outName, outRow);
                status=I_ISGR_ERR_ISGR_OUT_COR;
                break;
            }

        }

        if (outRow == 0l) {
            status=DALtableAddRows(*outTable, 0l, ptr_IBIS_events->numEvents, status);
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
            
        RILlogMessage(NULL, Log_0, "will write %li events",ptr_IBIS_events->numEvents, status);

        status=DALtablePutCol(*outTable, "ISGRI_PI", 0, DAL_BYTE, ptr_IBIS_events->numEvents,
                ptr_IBIS_events->isgri_pi, status);
        status=DALtablePutCol(*outTable, "ISGRI_ENERGY", 0, DAL_FLOAT, ptr_IBIS_events->numEvents,
                ptr_IBIS_events->isgri_energy, status);

        if (status != ISDC_OK) {
            RILlogMessage(NULL, Error_2, "Cannot write output data. Status=%d", status);
        };

        status=CommonStampObject(*outTable, "Energy correction.", status);

    } while(0);
    return status;
}


