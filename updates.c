/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "model_constants.h"
#include "model_structs.h"
#include "model_fns.h"
#include "model_globals_ext.h"
#include "binio.h"




void UpdateProbs(int DoPlace)
{
	int j, k, m, same_country, mc, mc_ind, q, p, ci;
	float t;

	if (!DoPlace)
	{
#pragma omp parallel for private(j,k,m,t) schedule(static,500)
		for (j = 0; j < P.NCP; j++)
		{
			CellLookup[j]->tot_prob = 0;
			CellLookup[j]->S0 = CellLookup[j]->S + CellLookup[j]->L + CellLookup[j]->I;
			if ((P.DoDeath) && (!P.DoSIS))
			{
				CellLookup[j]->S0 += CellLookup[j]->n / 5;
				if ((CellLookup[j]->n < 100) || (CellLookup[j]->S0 > CellLookup[j]->n)) CellLookup[j]->S0 = CellLookup[j]->n;
			}
		}
	}
	else
	{
#pragma omp parallel for private(j,k,m,t) schedule(static,500) //added j to private
		for (j = 0; j < P.NCP; j++)
		{
			CellLookup[j]->S0 = CellLookup[j]->S;
			CellLookup[j]->tot_prob = 0;
		}
	}
#pragma omp parallel for private(j,k,m,t,q,p,same_country,mc_ind,mc,ci) schedule(static,500)
	for (j = 0; j < P.NCP; j++)
	{

		if ((P.DoRoadNetwork) && (P.DoRoadPopEffect) && ((CellLookup[j]->road_connection) || (CellLookup[0]->road_connection)))
		{
			CellLookup[j]->cum_trans[0] = ((float)(CellLookup[0]->S0)) * CellLookup[j]->road_access * CellLookup[0]->road_access * CellLookup[j]->max_trans[0]; //scaling distance according to accessibility by road - was accessTo*accessFrom CellLookup[j]->road_access*CellLookup[0]->road_access
			t = ((float)CellLookup[0]->n) * CellLookup[j]->road_access * CellLookup[0]->road_access * CellLookup[j]->max_trans[0];
		}
		else if ((P.DoCapitalCityEffect) && (P.DoCapitalCityPopEffect) && (P.DoAdUnits) && ((CellLookup[j]->capital_city) ^ (CellLookup[0]->capital_city))) //if we're doing the capital city effect and either cell l or cell m (but not both) are in capital cities
		{
			same_country = 0;
			if (CellLookup[j]->capital_city)
			{
				//if cell j contains a capital city district, check to see whether any of cell 0 is in the same country
				ci = CellLookup[0] - Cells;
				mc = (ci / P.nch) * P.NMCL * P.nmch + (ci % P.nch) * P.NMCL;
				for (q = 0; q < P.NMCL; q++)
				{
					for (p = 0; p < P.NMCL; p++)
					{
						//get index of microcell
						mc_ind = mc + p + q * P.nmch;
						if ((AdUnits[Mcells[mc_ind].adunit].id / P.CountryDivisor) == (CellLookup[j]->capital_city / P.CountryDivisor))
						{
							same_country = 1;//P.RoadAccessDistance;
						}
					}
				}
			}
			else if (CellLookup[0]->capital_city)
			{
				//if cell 0 contains a capital city district, check to see whether any of cell j is in the same country
				ci = CellLookup[j] - Cells;
				mc = (ci / P.nch) * P.NMCL * P.nmch + (ci % P.nch) * P.NMCL;
				for (q = 0; q < P.NMCL; q++)
				{
					for (p = 0; p < P.NMCL; p++)
					{
						//get index of microcell
						mc_ind = mc + p + q * P.nmch;
						if ((AdUnits[Mcells[mc_ind].adunit].id / P.CountryDivisor) == (CellLookup[0]->capital_city / P.CountryDivisor))
						{
							same_country = 1;//P.RoadAccessDistance;
						}
					}
				}
			}
			if (same_country)
			{
				CellLookup[j]->cum_trans[0] = ((float)(CellLookup[0]->S0)) * P.CapitalCityPopEffect * CellLookup[j]->max_trans[0];
				t = ((float)CellLookup[0]->n) * P.CapitalCityPopEffect * CellLookup[j]->max_trans[0];
			}
			else
			{
				CellLookup[j]->cum_trans[0] = ((float)(CellLookup[0]->S0)) * CellLookup[j]->max_trans[0];
				t = ((float)CellLookup[0]->n) * CellLookup[j]->max_trans[0];
			}
		}
		else
		{
			CellLookup[j]->cum_trans[0] = ((float)(CellLookup[0]->S0)) * CellLookup[j]->max_trans[0];
			t = ((float)CellLookup[0]->n) * CellLookup[j]->max_trans[0];
		}


		//CellLookup[j]->cum_trans[0]=((float) (CellLookup[0]->S0))*CellLookup[j]->max_trans[0];
		//t=((float) CellLookup[0]->n)*CellLookup[j]->max_trans[0];
		for (m = 1; m < P.NCP; m++)
		{

			if ((P.DoRoadNetwork) && (P.DoRoadPopEffect) && ((CellLookup[j]->road_connection) || (CellLookup[m]->road_connection)))
			{
				CellLookup[j]->cum_trans[m] = CellLookup[j]->cum_trans[m - 1] + ((float)(CellLookup[m]->S0)) * CellLookup[j]->road_access * CellLookup[m]->road_access * CellLookup[j]->max_trans[m]; //scaling distance according to accessibility by road - was accessTo*accessFrom
				t += ((float)CellLookup[m]->n) * CellLookup[j]->road_access * CellLookup[m]->road_access * CellLookup[j]->max_trans[m];
			}
			else if ((P.DoCapitalCityEffect) && (P.DoCapitalCityPopEffect) && (P.DoAdUnits) && ((CellLookup[j]->capital_city) ^ (CellLookup[m]->capital_city))) //if we're doing the capital city effect and either cell l or cell m (but not both) are in capital cities
			{
				same_country = 0;
				if (CellLookup[j]->capital_city)
				{
					//if cell j contains a capital city district, check to see whether any of cell m is in the same country
					ci = CellLookup[m] - Cells;
					mc = (ci / P.nch) * P.NMCL * P.nmch + (ci % P.nch) * P.NMCL;
					for (q = 0; q < P.NMCL; q++)
					{
						for (p = 0; p < P.NMCL; p++)
						{
							//get index of microcell
							mc_ind = mc + p + q * P.nmch;
							if ((AdUnits[Mcells[mc_ind].adunit].id / P.CountryDivisor) == (CellLookup[j]->capital_city / P.CountryDivisor))
							{
								same_country = 1;//P.RoadAccessDistance;
							}
						}
					}
				}
				else if (CellLookup[m]->capital_city)
				{
					//if cell m contains a capital city district, check to see whether any of cell j is in the same country
					ci = CellLookup[j] - Cells;
					mc = (ci / P.nch) * P.NMCL * P.nmch + (ci % P.nch) * P.NMCL;
					for (q = 0; q < P.NMCL; q++)
					{
						for (p = 0; p < P.NMCL; p++)
						{
							//get index of microcell
							mc_ind = mc + p + q * P.nmch;
							if ((AdUnits[Mcells[mc_ind].adunit].id / P.CountryDivisor) == (CellLookup[m]->capital_city / P.CountryDivisor))
							{
								same_country = 1;//P.RoadAccessDistance;
							}
						}
					}
				}
				if (same_country)
				{
					CellLookup[j]->cum_trans[m] = CellLookup[j]->cum_trans[m - 1] + ((float)(CellLookup[m]->S0)) * P.CapitalCityPopEffect * CellLookup[j]->max_trans[m];
					t += ((float)CellLookup[m]->n) * P.CapitalCityPopEffect * CellLookup[j]->max_trans[m];
				}
				else
				{
					CellLookup[j]->cum_trans[m] = CellLookup[j]->cum_trans[m - 1] + ((float)(CellLookup[m]->S0)) * CellLookup[j]->max_trans[m];
					t += ((float)CellLookup[m]->n) * CellLookup[j]->max_trans[m];
				}
			}
			else
			{
				CellLookup[j]->cum_trans[m] = CellLookup[j]->cum_trans[m - 1] + ((float)(CellLookup[m]->S0)) * CellLookup[j]->max_trans[m];
				t += ((float)CellLookup[m]->n) * CellLookup[j]->max_trans[m];
			}


			//CellLookup[j]->cum_trans[m]=CellLookup[j]->cum_trans[m-1]+((float) (CellLookup[m]->S0))*CellLookup[j]->max_trans[m];
			//t+=((float) CellLookup[m]->n)*CellLookup[j]->max_trans[m];
		}
		CellLookup[j]->tot_prob = CellLookup[j]->cum_trans[P.NCP - 1];
		for (m = 0; m < P.NCP; m++)
			CellLookup[j]->cum_trans[m] /= CellLookup[j]->tot_prob;
		CellLookup[j]->tot_prob /= t;
		for (k = m = 0; k <= 1024; k++)
		{
			while (CellLookup[j]->cum_trans[m] * 1024 < ((float)k)) m++;
			CellLookup[j]->InvCDF[k] = m;
		}
	}
}


	/*
	 * Function: UpdateHospitals
	 *
	 * Purpose: updates hospitalisation parameters such as mean time to hospitalisation and number of beds per admin unit at each time step
	 * Parameters: none
	 * Returns: none
	 *
	 * Author: ggilani, 11/03/2017
	 */
void UpdateHospitals(double t)
{
	int i, j;//,nHospBeds,pop,newBeds;
	int nCases, numBedsTotal, numBedsInUse;
	double capacityFlag;
	//pop=P.N;

	//first update current mean time from onset to hospitalisation
	for (i = P.CurrIndMeanTimeToHosp; i < P.NMeanTimeToHosp; i++)
	{
		if (t >= P.ChangePointMeanTimeToHosp[i] + P.ETUTimeStart)
		{
			P.HospitalisationTime = P.MeanTimeToHosp[i];
			P.CurrIndMeanTimeToHosp++;
		}
	}

	//first update current mean time from onset to hospitalisation
	if (P.DoContactTracing)
	{
		for (i = P.CurrIndMeanTimeToHospCT; i < P.NMeanTimeToHospCT; i++)
		{
			if (t >= P.ChangePointMeanTimeToHospCT[i] + P.ETUTimeStart)
			{
				P.HospitalisationTime_contactTrace = P.MeanTimeToHospCT[i];
				P.CurrIndMeanTimeToHospCT++;
			}
		}
	}

	if ((P.DoReactETUBeds) && (t >= P.ETUTimeStart))
	{
		//new code for reactive provisioning of hospital beds - check to see if new beds are planned
		for (i = 0; i < P.NumAdunits; i++)
		{
			//go through each admin unit in term to see if we've passed the first threshold - based on active cases
			nCases = State.cumDC_adunit[i]; //calculated current number of active cases
			if (nCases >= P.InitCasesToETUBeds)
			{
				if (AdUnits[i].ETUbedsActive == 0)
				{
					//this means we have only reached the first threshold - set time to increase number of beds and number of beds
					AdUnits[i].nextTimeToETUBeds = t + P.InitDelayToETUBeds;
					//AdUnits[i].timeBedsAvailable = t + P.InitDelayToBeds;
					AdUnits[i].nextETUBeds = P.InitNumETUBeds;
					AdUnits[i].ETUbedsActive = 1;
				}
				else if ((AdUnits[i].ETUbedsActive) && (AdUnits[i].nextTimeToETUBeds < t))
				{
					//calculate current capacity
					numBedsTotal = BedsAvailablePerAdUnit(t, i);
					numBedsInUse = AdUnits[i].currentETUBeds;
					capacityFlag = (double)(numBedsInUse) / (double)(numBedsTotal);
					//if this has changed (more specifically, if it has got bigger because we've crossed the threshold for adding beds again)
					if (capacityFlag > P.CapacityToMoreETUBeds) //we also can't add more beds until we've added the last set
					{
						//set time for next lot of beds to be added
						AdUnits[i].nextTimeToETUBeds = t + P.SubDelayToETUBeds;
						AdUnits[i].nextETUBeds = P.SubNumETUBeds;
					}
				}
			}
		}
		// check to see if new beds should be allocated
		for (i = 0; i < P.NumAdunits; i++)
		{
			//check to see if new beds should be added
			if (((int)t == (int)AdUnits[i].nextTimeToETUBeds) && (AdUnits[i].ETUbedsActive == 1))
			{
				AdUnits[i].totalETUBeds += AdUnits[i].nextETUBeds;
				State.NumBeds += AdUnits[i].nextETUBeds;
				State.NumBeds_adunits[i] = AdUnits[i].totalETUBeds;
			}
		}
	}
	else
	{
		// check to see if new beds should be allocated
		for (i = 0; i < P.NumAdunits; i++)
		{
			//check to see if new beds should be added
			if ((t >= AdUnits[i].timeETUBedsAvailable) && (AdUnits[i].ETUbedsActive == 0))
			{
				AdUnits[i].totalETUBeds = AdUnits[i].initialETUBeds;
				AdUnits[i].ETUbedsActive = 1;
			}
		}
	}

	// something to get rid of capacity
	for (i = 0; i < P.NumAdunits; i++)
	{
		//check to see if new beds should be added
		if ((AdUnits[i].lastCaseDay > 0) && ((t - AdUnits[i].lastCaseDay) > P.DaysToRemoveCapacity))
		{
			AdUnits[i].totalETUBeds = 0;
			State.NumBeds -= AdUnits[i].totalETUBeds;
			State.NumBeds_adunits[i] = AdUnits[i].totalETUBeds;
		}
	}

	//Commented this out as we're changing to reactive provisioning of beds

	////secondly update the current number of hospital beds per admin unit
	//for(i=P.CurrIndHospBeds;i<P.NHospBeds;i++)
	//{
	//	if(t>=P.ChangePointHospBeds[i])
	//	{
	//		if(P.CurrIndHospBeds==0)
	//		{
	//			nHospBeds=P.HospBeds[i];
	//			for(j=0;j<P.NumAdunits;j++)
	//			{
	//				if(j==(P.NumAdunits-1))
	//				{
	//					AdUnits[j].totalBeds=nHospBeds;
	//				}
	//				else
	//				{
	//					AdUnits[j].totalBeds=ignbin(nHospBeds,((double) P.PopByAdunit[j][0])/((double) pop));
	//					nHospBeds=nHospBeds-AdUnits[j].totalBeds;
	//					pop=pop-P.PopByAdunit[j][0];
	//				}
	//			}
	//		}
	//		else
	//		{
	//			nHospBeds=P.HospBeds[i]-P.HospBeds[i-1]; //only allocate increased number of beds, so that we don't randomly reduce hospital beds in a place
	//			for(j=0;j<P.NumAdunits;j++)
	//			{
	//				if(j==(P.NumAdunits-1))
	//				{
	//					AdUnits[j].totalBeds+=nHospBeds;
	//				}
	//				else
	//				{
	//					newBeds=ignbin(nHospBeds,((double) P.PopByAdunit[j][0])/((double) pop));
	//					AdUnits[j].totalBeds+=newBeds;
	//					nHospBeds=nHospBeds-newBeds;
	//					pop=pop-P.PopByAdunit[j][0];
	//				}
	//			}

	//		}
	//		P.CurrIndHospBeds++;
	//	}
	//}

}


/*
 * Function: UpdateContactTracing
 *
 * Purpose: updates contact tracing capacity
 * Parameters: time t
 * Returns: none
 *
 * Author: ggilani, 11/03/2017
 */
void UpdateContactTracing(double t)
{
	int i, j;
	int nCT, numBedsTotal, numBedsInUse;
	double capacityFlag;
	//pop=P.N;

	if ((P.DoContactTracing) && (t >= P.ContactTracingTimeStart))
	{
		//code to see if contact tracing capacity should be increased
		for (i = 0; i < P.NumAdunits; i++)
		{
			if ((AdUnits[i].contactTraceThresholdCrossed) && (AdUnits[i].nextTimeToCT < t))
			{
				capacityFlag = (double)(AdUnits[i].nct) / (double)(AdUnits[i].contactTraceCapacity);
				//if this has changed (more specifically, if it has got bigger because we've crossed the threshold for adding beds again)
				if (capacityFlag > P.CapacityToMoreCT) //we also can't add more beds until we've added the last set
				{
					//set time for next lot of beds to be added
					AdUnits[i].nextTimeToCT = t + P.DelayToCT;
				}
			}
		}
		// check to see if new beds should be allocated
		for (i = 0; i < P.NumAdunits; i++)
		{
			//check to see if new beds should be added
			if (((int)t == (int)AdUnits[i].nextTimeToCT) && (AdUnits[i].contactTraceCaseThreshold == 1))
			{
				AdUnits[i].contactTraceCapacity += AdUnits[i].contactTraceCapacityInc;
				//State.NumBeds += AdUnits[i].nextETUBeds;
				//State.NumBeds_adunits[i] = AdUnits[i].totalETUBeds;
			}
		}
	}

	// something to get rid of capacity
	for (i = 0; i < P.NumAdunits; i++)
	{
		//check to see if new beds should be added
		if ((AdUnits[i].lastCaseDay > 0) && ((t - AdUnits[i].lastCaseDay) > P.DaysToRemoveCapacity))
		{
			AdUnits[i].contactTraceCapacity = 0;
		}
	}
}

/*
 * Function: UpdateSDB
 *
 * Purpose: updates maximum number of safe and dignified burials per admin unit
 * Parameters: time t
 * Returns: none
 *
 * Author: ggilani, 11/03/2017
 */
void UpdateSDB(double t)
{
	int i, j;
	int nCT, numBedsTotal, numBedsInUse;
	double capacityFlag;
	//pop=P.N;

	if (t >= P.FuneralControlTimeStart)
	{
		//code to see if contact tracing capacity should be increased
		for (i = 0; i < P.NumAdunits; i++)
		{
			if (AdUnits[i].nextTimeToSDB < t)
			{
				capacityFlag = (double)(State.cumSDB_adunit[i]) / (double)(AdUnits[i].maxSDB);
				//if this has changed (more specifically, if it has got bigger because we've crossed the threshold for adding beds again)
				if (capacityFlag > P.CapacityToMoreSDB) //we also can't add more beds until we've added the last set
				{
					//set time for next lot of beds to be added
					AdUnits[i].nextTimeToSDB = t + P.DelayToSDB;
				}
			}
		}
		// check to see if new beds should be allocated
		for (i = 0; i < P.NumAdunits; i++)
		{
			//check to see if new beds should be added
			if ((int)t == (int)AdUnits[i].nextTimeToSDB)
			{
				AdUnits[i].maxSDB += P.incCapacitySDB;
				//State.NumBeds += AdUnits[i].nextETUBeds;
				//State.NumBeds_adunits[i] = AdUnits[i].totalETUBeds;
			}
		}
	}

	// something to get rid of capacity
	for (i = 0; i < P.NumAdunits; i++)
	{
		//check to see if new beds should be added
		if ((AdUnits[i].lastCaseDay > 0) && ((t - AdUnits[i].lastCaseDay) > P.DaysToRemoveCapacity))
		{
			AdUnits[i].maxSDB = 0;
		}
	}
}

/*
 * Function: UpdateVaccination
 *
 * Purpose: updates vaccination parameters such as number of rings to vaccinate, proportion of rings to vaccinate, at each time step
 * Parameters: none
 * Returns: none
 *
 * Author: ggilani, 29/05/2019
 */
void UpdateVaccination(double t, int n)
{

	int i, weeklyDC;//,nHospBeds,pop,newBeds;
	double avDailyDC; //average daily cases


	//first update current proportion of ring to vaccinate
	for (i = P.CurrIndPropRingVacc; i < P.NPropRingVacc; i++)
	{
		if (t >= P.ChangePointPropRingVacc[i])
		{
			P.PropRingVacc = P.ListPropRingVacc[i];
			P.CurrIndPropRingVacc++;
		}
	}

	if ((t >= P.TimeToIncVaccRing) && (P.NVaccRingsActive < 3))
	{
		P.NVaccRingsActive++;
	}

	//update vacc dose per day
	//first find total number of detected cases over the past 7 days
	if ((t >= P.VaccTimeStart) && (P.UpdateVaccDosePerDay))
	{
		weeklyDC = 0;
		for (i = max(n - 6, 0); i <= n; i++)
		{
			weeklyDC += TimeSeries[i].incDC;
		}
		avDailyDC = (double)weeklyDC / 7;
		if (avDailyDC > P.VaccDoseFlag)
		{
			P.VaccDoseFlag++;
			P.VaccDosePerDay = P.VaccDoseFlag * P.BaseVaccDosePerDay;
			P.VaccGeoDosePerDay = P.VaccDoseFlag * P.BaseVaccGeoDosePerDay;
		}
		if (P.VaccDosePerDay > P.MaxVaccDosePerDay)
		{
			P.VaccDosePerDay = P.MaxVaccDosePerDay;
			P.VaccGeoDosePerDay = P.MaxVaccGeoDosePerDay;
			P.UpdateVaccDosePerDay = 0;
		}
	}


	// something to get rid of capacity

	//check to see if new beds should be added
	if ((P.DayExtinct > 0) && ((t - P.DayExtinct) > P.DaysToRemoveCapacity))
	{
		P.VaccDosePerDay = 0;
		P.VaccGeoDosePerDay = 0;
	}


}


void UpdateCaseDetection(double t)
{
	int i, j;

	if (P.DoUpdateCaseDetectionByTime)
	{
		for (i = P.CurrIndUpdateCaseDetect; i < P.NUpdateCaseDetection; i++)
		{
			if (t >= P.TimeToUpdateCaseDetection[i])
			{
				if (P.DoAdUnits)
				{
					for (j = 0; j < P.NumAdunits; j++)
					{
						AdUnits[j].caseDetectRate = P.ListUpdateCaseDetection[i];
					}
				}
				else
				{
					P.CaseDetectionRate = P.ListUpdateCaseDetection[i];
				}
				P.CurrIndUpdateCaseDetect++;
			}
		}
	}
	else if ((P.DoUpdateCaseDetectionByCases) && (P.UpdateCaseDetectionByCasesFlag == 0))
	{
		if (State.cumDC > P.CaseThresholdUntilUpdateCaseDetection)
		{
			P.UpdateCaseDetectionByCasesFlag = 1;
			if (P.DoAdUnits)
			{
				for (j = 0; j < P.NumAdunits; j++)
				{
					AdUnits[j].caseDetectRate = P.CaseDetectionRateAfterThresholdReached;
				}
			}
			else
			{
				P.CaseDetectionRate = P.CaseDetectionRateAfterThresholdReached;
			}
		}
	}

	//if (t >= P.TimeToUpdateCaseDetection)
	//{
	//	if(P.DoAdUnits)
	//	{
	//		for (i = 0; i < P.NumAdunits; i++)
	//		{
	//			AdUnits[i].caseDetectRate = P.UpdatedCaseDetectionRate;
	//		}
	//	}
	//	else
	//	{
	//		P.CaseDetectionRate = P.UpdatedCaseDetectionRate;
	//	}
	//}
}
