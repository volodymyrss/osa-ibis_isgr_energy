/************************************************************************
 * FILE:        ibis_isgr_energy.h
 * VERSION:     8.3
 * COMPONENT:   ibis_isgr_energy
 * AUTHOR:      A.Sauvageon & C. Couvreur,     SAp-CEA
 * DESCRIPTION: definition of the globals and functions' prototypes
 * HISTORY:
 *   MD,     28/03/1999, template version
 *   MD,     09/04/1999, demi-module version
 *   PL,     20/05/1999, first module version
 *   PL,     09/06/1999, full module version
 *   PL,     17/11/1999, add common_prepare
 *   PL,     27/03/2002, change type of caltable from integer to real*4
 *   PL,     22/04/2002, add groups for IC files
 *   PL,3.6, 19/09/2002, put the random number generator as a parameter
 *   PL,4.0, 30/07/2003, remove call to random number no longer useful 
 *   PL,4.2, 07/10/2003, SCREW 1283 (get ISGRI mean temperature from HK)
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
 *  ASA,5.5, 26/04/2005, SPR 4112 (offset decrease, better gain decrease)
 *
 *  New version 5.6.0 made by Christophe Couvreur (CC)
 *  CC, 5.6.1, 28/07/2006, SPR 4537 (ibis_isgr_energy uses too much memory)
 *  CC, 5.6.2, 02/08/2006, SPR 4549 (in function ibis_isgr_energyTransformNEW
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
 *  CC, 6.0                SPR 4695 (calibration is based on 2 laws to avoid increase 
 *                         of Crab flux in 20-40kev). 
 *                         Below 50kev we apply the first law and above, 
 *                         we apply the second one. In the first law, Ph gain and Offset 
 *                         don't depend on RT contrary to the second law.
 *                         - Switch on Time routine is no more considerated.
 *                         - The Table for rt effect has disapeared and included into
 *                           PH gain and offset for 2nd law. 
 *  CC, 6.1                SPR 4696
 * ASA, 6.2,  20/12/2007,  SPR 4773 4698 (+ cosmetic changes)
 * ASA, 6.5   16/03/2009,   SPR 4838 (+ cosmetic changes = ibis_comp_energy)
 *  NP, 6.6  14/04/2009    SPR 4855, revert SPR 4830
 *  PL, 8.0  02/02/2011,   remove IREM counters, Temperature correction by MDU
 *  PL, 8.2  02/04/2012,   modify coefficients PAR1_..._corrPH1 to correct <50 keV behavior
 *  CF, 8.3  05/11/2015,   SCREW 2624 path for low energy correction implemented with Paris agreement
 ************************************************************************/

#ifndef IBIS_ENERGY_H_INCLUDED
#define IBIS_ENERGY_H_INCLUDED

#include "isdc.h"
#include "dal3ibis.h"
#include "dal3hk.h"

#define COMPONENT_NAME           "ibis_isgr_energy"
#define COMPONENT_VERSION        "9.0"

#define I_ISGR_ERR_MEMORY         -122050
#define I_ISGR_ERR_BAD_INPUT      -122051
#define I_ISGR_ERR_ISGR_OFFS_BAD  -122052
#define I_ISGR_ERR_ISGR_RISE_BAD  -122053
#define I_ISGR_ERR_ISGR_OUT_COR   -122054
#define I_ISGR_ERR_IBIS_IREM_BAD  -122055
#define I_ISGR_ERR_ISGR_PHGO2_BAD -122056

#define ISGRI_N_PIX     16384l
#define ISGRI_GO_N_COL      5
#define ISGRI_DIM_LUT2_3D   3
#define ISGRI_RT_N_ENER_SCALED 1024
#define ISGRI_RT_N_DATA   256
#define ISGRI_RT_N_RANDOM_DIM   500

#define ISGRI_PHG2_N_COL    2
#define ISGRI_PHO2_N_COL    3

/* unused since V 6.0
#define ISGRI_RT_N_COL      3
*/ 
/* unused dimensions since V 5.6.0
#define ISGRI_RT_N_PAR   257
#define ISGRI_RT_N_ENER 2048
#define ISGRI_RT_N_TYPE    1
*/

#define DS_ISGR_RAW       "ISGR-EVTS-ALL"
#define DS_ISGR_GO        "ISGR-OFFS-MOD"
#define DS_ISGR_3DL2_MOD  "ISGR-3DL2-MOD"
#define DS_PHG2           "ISGR-GAIN-MOD"
#define DS_PHO2           "ISGR-OFF2-MOD"
/* unused structures since V 6.0
#define DS_ISGR_RISE_PRO  "ISGR-RISE-PRO"
#define DS_ISGR_SWIT      "IBIS-SWIT-CAL"
*/
/* unused structures since V 5.6.0
#define DS_ENER_MOD       "ISGR-DROP-MOD"
#define DS_ISGR_RT        "ISGR-RISE-MOD"
*/
#define KEY_COL_OUT  "ISGRI_PI"

#define DS_ISGR_HK   "IBIS-DPE.-CNV"
#define KEY_MCE_BIAS "I0E_MCDTE_MBIAS"
#define KEY_DEF_BIAS      -120.0    /* default when HK1 missing */
#define KEY_MAX_BIAS       -60.0    /*  155 to disregard RAW 0  */
#define KEY_MCE_TEMP "I0E_MTEMP2_MMDU"
#define KEY_DEF_TEMP        -8.0    /* default when HK1 missing */
#define KEY_RMS_TEMP         1.2    /* disregard MDU Temp */
#define KEY_MIN_TEMP       -50.5    /* to disregard RAW 0 */

#define SEC_DELTA_MIN       25      /* minimal time range to search HK1 */

/* constant parameters for the energy correction */
#define OFF_SCALE0          -1.997
#define G_SCALE0             1.0184
#define G_SCALE1             0.0000089
#define FINALOFFS            0.
#define PAR1_GAIN_corrPH1    2.047
#define PAR2_GAIN_corrPH1   -0.00061
#define PAR1_OFFSET_corrPH1 -5.655


typedef struct {
    char riseDOLstr[DAL_FILE_NAME_STRING]; 
    char acorDOLstr[DAL_FILE_NAME_STRING];
    char phGainDOLstr[DAL_FILE_NAME_STRING];
    char phOffsDOLstr[DAL_FILE_NAME_STRING];
} ISGRI_energy_caldb_dols_struct;

typedef struct  {
    int makeUnique,
        clobber,
        gti,
        erase, chatter;
} ibis_isgr_energy_settings_struct;

int ibis_isgr_energyWork(dal_element *workGRP,
                         int          gti,
                         int          erase,
                         int          chatter,
                         char        *acorName,
                         char        *riseName,
                         char        *phGainDOLstr,
                         char        *phOffDOLstr);

int ibis_isgr_energyCheckIn(
                         char         *acorName,
                         char         *riseName,
                         char         *phGainDOLstr,
                         char         *phOffDOLstr,
                         dal_element **isgrOffsTabPtr,
                         dal_element **isgrRiseTabPtr,
                         dal_element **supGainTabPtr,
                         dal_element **supOffsTabPtr);

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
                         DAL3_Byte **isZ);

int ibis_isgr_energyCheckOut(dal_element *workGRP,
                         char         *outName,
                         int           chatter,
                         long          numEvents,
                         int           erase,
                         dal_element **outTable);

int ibis_isgr_energyReadCal(dal_element *isgrOffsTabPtr,
                         dal_element    *isgrRiseTabPtr,
                         dal_element    *supGainTabPtr,
                         dal_element    *supOffsTabPtr,
                         int             chatter,
                         double        **isgriGoTab,
                         int            *isgriPixTab,
                         short          *isgriRtTab,
			    char        *phGainDOLstr,
			    char        *phOffDOLstr,
                         double         *supCoeffg,
                         double         *supCoeffo);

int ibis_energyIsgrHkCal(dal_element *workGRP,
                         OBTime       obtStart,
                         OBTime       obtEnd,
                         double       meanT[8],
                         double       meanBias[8],
                         int          chatter,
                         int          status);

int ibis_isgr_energyTransform(dal_element *outTable,
                         long         numEvents,
                         long         revol,
                         double     **isgriGoTab,
                         int         *isgriPixTab,
                         short       *isgriRtTab,
                         double      *supCoeffg,
                         double      *supCoeffo,
                         int          chatter,
                         double       meanT[8],
	       		 double       meanBias[8],
                         DAL3_Word   *isgriPha,
                         DAL3_Byte   *riseTime,
                         DAL3_Byte   *isgriY,
                         DAL3_Byte   *isgriZ,
                         DAL3_Byte   *isgrPi,
                         float       *isgrEnergy);
#endif
