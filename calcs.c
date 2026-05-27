/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "model_constants.h"
#include "model_structs.h"
#include "model_fns.h"
#include "model_globals_ext.h"
#include "binio.h"


	/* function BedsAvailablePerAdUnit
	 *
	 * Purpose: to return the number of beds available at a given time
	 *
	 * Parameters:
	 *	double t: the current time
	 *	int auindex: the index of the admin unit to check
	 *
	 * Returns: an integer number of beds available in that AU at time t
	 *
	 * Author: ggilani, Date: 29/10/14
	 */
int BedsAvailablePerAdUnit(double t, int auindex)
{
	return AdUnits[auindex].totalETUBeds;

	/*if(t>=AdUnits[auindex].timeBedsAvailable)
	{
		return AdUnits[auindex].totalBeds;
	}
	else
	{
		return 0;
	}*/
}


double CalcHouseInf(int j, unsigned short int ts)
{
	return (HOST_ISOLATED(j) ? P.CaseIsolationHouseEffectiveness : 1.0)
		* (HOST_QUARANTINED(j) ? P.HQuarantineHouseEffect : 1.0)
		* P.HouseholdDenomLookup[Households[Hosts[j].hh].nhr - 1] * CalcPersonInf(j, ts);
}

double CalcPlaceInf(int j, int k, unsigned short int ts)
{

	return (HOST_ISOLATED(j) ? P.CaseIsolationEffectiveness : 1.0)
		* (HOST_QUARANTINED(j) ? P.HQuarantinePlaceEffect[k] : 1.0)
		* ((Hosts[j].inf == -2) ? P.SymptPlaceTypeContactRate[k] : 1.0)
		* P.PlaceTypeTrans[k] / P.PlaceTypeGroupSizeParam1[k] * CalcPersonInf(j, ts);
}

double CalcSpatialInf(int j, unsigned short int ts)
{
	return (HOST_ISOLATED(j) ? P.CaseIsolationEffectiveness : 1.0)
		* (HOST_QUARANTINED(j) ? P.HQuarantineSpatialEffect : 1.0)
		* ((Hosts[j].inf == -2) ? P.SymptSpatialContactRate : 1.0)
		* P.RelativeSpatialContact[HOST_AGE_GROUP(j)]
		* CalcPersonInf(j, ts); 		/*	*Hosts[j].spatial_norm */
}

double CalcPersonInf(int j, unsigned short int ts)
{

	return (HOST_TREATED(j) ? P.EvolResistRelTreatInfDrop[Hosts[j].resist] : 1.0)
		* (HOST_VACCED(j) ? P.VaccInfDrop : 1.0)
		* fabs(Hosts[j].infectiousness)
		* Hosts[j].infectiousMult
		* Hosts[j].base_inf_level
		* P.EvolResistRelInf[Hosts[j].resist]
		* P.infectiousness[ts - Hosts[j].latent_time - 1]
		* (Hosts[j].etu ? P.RelativeInfectiousnessETU : 1.0)
		* ((!Hosts[j].etu && Hosts[j].contactTraced) ? P.RelativeInfectiousnessContactTraced : 1.0);
}

double CalcHouseSusc(int ai, unsigned short int ts, int infector, int tn)
{
	return CalcPersonSusc(ai, ts, infector, tn)
		* ((Mcells[Hosts[ai].mcell].socdist == 2) ? ((Hosts[ai].esocdist_comply) ? P.ESocDistHouseholdEffect : P.SocDistHouseholdEffect) : 1.0);
}

double CalcPlaceSusc(int ai, int k, unsigned short int ts, int infector, int tn)
{
	return CalcPersonSusc(ai, ts, infector, tn) * (HOST_QUARANTINED(ai) ? P.HQuarantinePlaceEffect[k] : 1.0)
		* ((Mcells[Hosts[ai].mcell].socdist == 2) ? ((Hosts[ai].esocdist_comply) ? P.ESocDistPlaceEffect[k] : P.SocDistPlaceEffect[k]) : 1.0);
}

double CalcSpatialSusc(int ai, unsigned short int ts, int infector, int tn)
{
	return CalcPersonSusc(ai, ts, infector, tn) * (HOST_QUARANTINED(ai) ? P.HQuarantineSpatialEffect : 1.0)
		* ((Mcells[Hosts[ai].mcell].socdist == 2) ? ((Hosts[ai].esocdist_comply) ? P.ESocDistSpatialEffect : P.SocDistSpatialEffect) : 1.0);
}


double CalcPersonSusc(int ai, unsigned short int ts, int infector, int tn)
{
	double suscdrop;

	if ((Hosts[ai].keyworker) && ((Hosts[ai].vacc_start_time <= 0) || (HOST_TO_BE_VACCED(ai) && Hosts[ai].revacc)))
	{
		suscdrop = P.HCWVaccSuscDrop;
	}
	else
	{
		suscdrop = (HOST_VACCED(ai) ? (HOST_VACCED_SWITCH(ai) ? P.VaccSuscDrop2 : P.VaccSuscDrop) : 1.0);
	}

	return P.WAIFW_Matrix[HOST_AGE_GROUP(ai)][HOST_AGE_GROUP(infector)]
		* P.AgeSusceptibility[HOST_AGE_GROUP(ai)] * Hosts[ai].susc
		* (HOST_TREATED(ai) ? P.EvolResistRelProphSusc[Hosts[infector].resist] : 1.0)
		* suscdrop * CalcPrevalenceDepTransmission(ai)
		* (((Hosts[ai].hcw || Hosts[ai].flw) && P.OutbreakDetected) ? P.RelSuscPPE : 1.0);
}


/* function: CalcPrevalenceDepTransmission(int ai)
 *
 * purpose: to determine whether susceptibility should be scaled due to prevalence dependent transmission and to return the relative reduction in susceptibility
 *
 * parameters:
 *	int ai, the index of the current host
 *
 * returns: relative reduction in susceptibility
 *
 * author: ggilani, date: 22/12/14
 */
double CalcPrevalenceDepTransmission(int ai)
{
	if (P.DoPrevDependTrans)
	{
		if (P.DoPrevDepTransTotalCases)
		{
			if (((double)(Cells[Hosts[ai].pcell].D + Cells[Hosts[ai].pcell].R)) / (double)Cells[Hosts[ai].pcell].n >= P.PrevDependTransThresh)
			{
				return P.PrevDependRelativeSusc;
			}
			else
			{
				return 1.0;
			}
		}
		else if (P.DoPrevDepTransCurrInf)
		{
			if ((double)Cells[Hosts[ai].pcell].I / (double)Cells[Hosts[ai].pcell].n >= P.PrevDependTransThresh)
			{
				return P.PrevDependRelativeSusc;
			}
			else
			{
				return 1.0;
			}
		}
		else
		{
			return 1.0;
		}
	}
	else
	{
		return 1.0;
	}
}

