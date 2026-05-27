/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "model_constants.h"
#include "model_structs.h"
#include "model_fns.h"
#include "model_globals_ext.h"
#include "binio.h"


void RecordSample(double t, int n)
{
	int i, j, k, S, L, I, R, D, cumC, cumTC, cumI, cumR, cumD, cumDC, cumDD, cumSDB, cumFC, cumFI; //added cumulative funeral infections: ggilani 24/10/14
	int cumETU; //add number of ETU cases, cumulative ETU cases: ggilani 28/10/14
	int cumH; //added cumulative number of hospitalisations - ggilani 01/07/24
	int cumCT; //added cumulative number of contact traced: ggilani 15/06/17
	int cumCC; //added cumulative number of cases who are contacts: ggilani 28/05/2019
	int cumHQ, cumAC, cumAH, cumAA, cumACS, cumAPC, cumAPA, cumAPCS, numPC, trigDC;
	int cumC_country[MAX_COUNTRIES]; //add cumulative cases per country
	int cumC_adunit[MAX_ADUNITS]; //added cumulative cases per adunit
	int cumDC_adunit[MAX_ADUNITS]; //added cumulative detected cases per adunit
	int cumDD_adunit[MAX_ADUNITS]; //added cumulative detected deaths per adunit <-this and cumulative detected recoveries will help calculate current active detected cases per adunit to use when determine bed capacity increases
	int cumDR_adunit[MAX_ADUNITS]; //added cumulative detected recoveries per adunit - ggilani 13/09/23
	int cumSDB_adunit[MAX_ADUNITS];
	int cumV_adunit[MAX_ADUNITS];
	int cumVG_adunit[MAX_ADUNITS];
	int cumD_adunit[MAX_ADUNITS];
	cell* ct;
	unsigned short int ts;
	int nCase; //added to calculate triggers for individual admin unit funeral controls - ggilani 07/03/2017
	float nsy, * temp2, * temp3, * temp7;

	nsy = DAYS_PER_YEAR / P.SampleStep;
	ts = (unsigned short int) (P.TimeStepsPerDay * t);
	S = L = I = R = D = cumI = cumC = cumDC = cumDD = cumSDB = cumTC = cumFC = cumFI = cumHQ = cumAC = cumAA = cumAH = cumACS = cumAPC = cumAPA = cumAPCS = cumD = cumETU = cumH = cumCT = cumCC = 0;
	for (i = 0; i < MAX_COUNTRIES; i++) cumC_country[i] = 0;
	for (i = 0; i < MAX_ADUNITS; i++) cumC_adunit[i] = cumDC_adunit[i] = cumD_adunit[i] = cumDD_adunit[i] = cumDR_adunit[i] = cumSDB_adunit[i] = cumV_adunit[i] = cumVG_adunit[i] = 0; //added cumulative detected cases, detected deaths, detected recoveries
#pragma omp parallel for private(i,ct) schedule(static,10000) reduction(+:S,L,I,R,D,cumTC) //added i to private
	for (i = 0; i < P.NCP; i++)
	{
		ct = CellLookup[i];
		S += (int)ct->S;
		L += (int)ct->L;
		I += (int)ct->I;
		R += (int)ct->R;
		D += (int)ct->D;
		cumTC += (int)ct->cumTC;
	}
	cumR = R;
	//cumD=D;
	State.sumRad2 = 0;
	for (j = 0; j < P.NumThreads; j++)
	{
		cumI += StateT[j].cumI;
		cumC += StateT[j].cumC;
		cumDC += StateT[j].cumDC;
		cumFC += StateT[j].cumFC;
		cumFI += StateT[j].cumFI; //added cumulative Funeral infections
		cumETU += StateT[j].cumETU; //added cumulative hospitalisation
		cumH += StateT[j].cumH;
		cumCT += StateT[j].cumCT; //added contact tracing
		cumCC += StateT[j].cumCC; //added cases who are contacts
		State.sumRad2 += StateT[j].sumRad2;
		cumHQ += StateT[j].cumHQ;
		cumAC += StateT[j].cumAC;
		cumAA += StateT[j].cumAA;
		cumAPC += StateT[j].cumAPC;
		cumAPA += StateT[j].cumAPA;
		cumAPCS += StateT[j].cumAPCS;
		cumAH += StateT[j].cumAH;
		cumACS += StateT[j].cumACS;
		cumD += StateT[j].cumD;
		cumDD += StateT[j].cumDD;
		cumSDB += StateT[j].cumSDB;
		//add up cumulative country counts: ggilani - 12/11/14
		for (i = 0; i < MAX_COUNTRIES; i++) cumC_country[i] += StateT[j].cumC_country[i];
		//add up cumulative adunit counts: ggilani - 12/11/14
		for (i = 0; i < P.NumAdunits; i++) cumC_adunit[i] += StateT[j].cumC_adunit[i];
		for (i = 0; i < P.NumAdunits; i++) cumDC_adunit[i] += StateT[j].cumDC_adunit[i];
		for (i = 0; i < P.NumAdunits; i++) cumD_adunit[i] += StateT[j].cumD_adunit[i];
		for (i = 0; i < P.NumAdunits; i++) cumSDB_adunit[i] += StateT[j].cumSDB_adunit[i];
		for (i = 0; i < P.NumAdunits; i++) cumV_adunit[i] += StateT[j].cumV_adunit[i];
		for (i = 0; i < P.NumAdunits; i++) cumVG_adunit[i] += StateT[j].cumVG_adunit[i]; //added geographically targeted vaccination by admin_unit
		for (i = 0; i < P.NumAdunits; i++) cumDD_adunit[i] += StateT[j].cumDD_adunit[i]; //detected deaths
		for (i = 0; i < P.NumAdunits; i++) cumDR_adunit[i] += StateT[j].cumDR_adunit[i]; //detected recoveries
		if (State.maxRad2 < StateT[j].maxRad2) State.maxRad2 = StateT[j].maxRad2;
	}
	for (j = 0; j < P.NumThreads; j++)
		StateT[j].maxRad2 = State.maxRad2;
	TimeSeries[n].t = t;
	TimeSeries[n].S = (double)S;
	TimeSeries[n].L = (double)L;
	TimeSeries[n].I = (double)I;
	TimeSeries[n].R = (double)R;
	TimeSeries[n].D = (double)D;
	TimeSeries[n].incI = (double)(cumI - State.cumI);
	TimeSeries[n].incC = (double)(cumC - State.cumC);
	TimeSeries[n].incFC = (double)(cumFC - State.cumFC);
	TimeSeries[n].incFI = (double)(cumFI - State.cumFI); //added funeral infections
	TimeSeries[n].incETU = (double)(cumETU - State.cumETU); //added incidence of hospitalisation
	TimeSeries[n].incH = (double)(cumH - State.cumH);
	TimeSeries[n].incCT = (double)(cumCT - State.cumCT); // added contact tracing
	TimeSeries[n].incCC = (double)(cumCC - State.cumCC); // added cases who are contacts
	TimeSeries[n].incDC = (double)(cumDC - State.cumDC); //added incidence of detected cases
	TimeSeries[n].incD = (double)(cumD - State.cumD); //added incidence of detected deaths
	TimeSeries[n].incDD = (double)(cumDD - State.cumDD); //added incidence of detected deaths
	TimeSeries[n].incSDB = (double)(cumSDB - State.cumSDB); //added incidence of safe burials
	TimeSeries[n].incTC = (double)(cumTC - State.cumTC);
	TimeSeries[n].incR = (double)(cumR - State.cumR);
	TimeSeries[n].incD = (double)(cumD - State.cumD);
	TimeSeries[n].incHQ = (double)(cumHQ - State.cumHQ);
	TimeSeries[n].incAC = (double)(cumAC - State.cumAC);
	TimeSeries[n].incAH = (double)(cumAH - State.cumAH);
	TimeSeries[n].incAA = (double)(cumAA - State.cumAA);
	TimeSeries[n].incACS = (double)(cumACS - State.cumACS);
	TimeSeries[n].incAPC = (double)(cumAPC - State.cumAPC);
	TimeSeries[n].incAPA = (double)(cumAPA - State.cumAPA);
	TimeSeries[n].incAPCS = (double)(cumAPCS - State.cumAPCS);
	TimeSeries[n].cumT = State.cumT;
	TimeSeries[n].cumUT = State.cumUT;
	TimeSeries[n].cumTP = State.cumTP;
	TimeSeries[n].cumV = State.cumV;
	if (t > P.VaccTimeStart)
	{
		TimeSeries[n].capV = P.VaccDosePerDay;
		TimeSeries[n].capVG = P.VaccGeoDosePerDay;
	}
	else
	{
		TimeSeries[n].capV = 0;
		TimeSeries[n].capVG = 0;
	}
	TimeSeries[n].cumVG = State.cumVG; //added VG;

	TimeSeries[n].cumDC = cumDC;
	TimeSeries[n].cumDD = cumDD;
	TimeSeries[n].cumSDB = cumSDB;
	//incidence per country
	for (i = 0; i < MAX_COUNTRIES; i++) TimeSeries[n].incC_country[i] = (double)(cumC_country[i] - State.cumC_country[i]);
	trigDC = cumDC;
	if (n >= P.TriggersSamplingInterval) trigDC -= TimeSeries[n - P.TriggersSamplingInterval].cumDC;
	if (trigDC > State.trigDC) State.trigDC = trigDC;
	State.S = S;
	State.L = L;
	State.I = I;
	State.R = R;
	State.D = D;
	State.cumI = cumI;
	State.cumDC = cumDC;
	State.cumDD = cumDD;
	State.cumSDB = cumSDB;
	State.cumTC = cumTC;
	State.cumFC = cumFC;
	State.cumFI = cumFI; //added cumulative funeral infections
	State.cumETU = cumETU; //added cumulative hospitalisation
	State.cumH = cumH;
	State.cumCT = cumCT; //added cumulative contact tracing
	State.cumCC = cumCC; //added cumulative cases who are contacts
	State.cumC = cumC;
	State.cumR = cumR;
	State.cumD = cumD;
	State.cumHQ = cumHQ;
	State.cumAC = cumAC;
	State.cumAH = cumAH;
	State.cumAA = cumAA;
	State.cumACS = cumACS;
	State.cumAPC = cumAPC;
	State.cumAPA = cumAPA;
	State.cumAPCS = cumAPCS;
	//update cumulative cases per country
	for (i = 0; i < MAX_COUNTRIES; i++) State.cumC_country[i] = cumC_country[i];
	//update overall state variable for cumulative cases per adunit
	for (i = 0; i < P.NumAdunits; i++) State.cumC_adunit[i] += cumC_adunit[i]; //here we increment rather than set equal to as the cumulative admin unit counts per thread are reset to zero after recording a sample
	for (i = 0; i < P.NumAdunits; i++) State.cumDC_adunit[i] += cumDC_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumD_adunit[i] += cumD_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumDD_adunit[i] += cumDD_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumDR_adunit[i] += cumDR_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumDD_adunit[i] += cumDD_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumSDB_adunit[i] += cumSDB_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumV_adunit[i] += cumV_adunit[i]; //added vaccination: ggilani 14/01/25
	for (i = 0; i < P.NumAdunits; i++) State.cumVG_adunit[i] += cumVG_adunit[i]; //added geographically targeted vaccination: ggilani 22/08/25
	TimeSeries[n].rmsRad = (State.cumI > 0) ? sqrt(State.sumRad2 / ((double)State.cumI)) : 0;
	TimeSeries[n].maxRad = sqrt(State.maxRad2);
	TimeSeries[n].extinct = ((((P.SmallEpidemicCases >= 0) && (State.R <= P.SmallEpidemicCases)) || (P.SmallEpidemicCases < 0)) && (State.I + State.L == 0)) ? 1 : 0;
	if ((TimeSeries[n].extinct == 1) && (P.DayExtinct == 0))
	{
		P.DayExtinct = (double)n;
		P.VaccNewCoursesEndTime = (double)n;
	}
	TimeSeries[n].detected = P.OutbreakDetected;
	for (i = 0; i < NUM_AGE_GROUPS; i++)
	{
		TimeSeries[n].incCa[i] = TimeSeries[n].incIa[i] = TimeSeries[n].incDa[i] = TimeSeries[n].incDCa[i] = TimeSeries[n].incETUa[i] = TimeSeries[n].incHa[i] = TimeSeries[n].incVa[i] = 0;
		for (j = 0; j < P.NumThreads; j++)
		{
			TimeSeries[n].incCa[i] += (double)StateT[j].cumCa[i];
			TimeSeries[n].incIa[i] += (double)StateT[j].cumIa[i];
			TimeSeries[n].incDa[i] += (double)StateT[j].cumDa[i];
			TimeSeries[n].incDCa[i] += (double)StateT[j].cumDCa[i];
			TimeSeries[n].incETUa[i] += (double)StateT[j].cumETUa[i];
			TimeSeries[n].incHa[i] += (double)StateT[j].cumHa[i];
			//TimeSeries[n].incVa[i] += (double)StateT[j].cumVa[i];
		}
		TimeSeries[n].incVa[i] += State.cumVa[i]; //changed vacc by age to here, because it is incremented outside of the threads: ggilani 22/02/22
	}
	for (i = 0; i < P.EvolResistNumTypes; i++)
	{
		TimeSeries[n].incC_resist[i] = TimeSeries[n].incI_resist[i] = TimeSeries[n].cumT_resist[i] = 0;
		for (j = 0; j < P.NumThreads; j++)
		{
			TimeSeries[n].incC_resist[i] += (double)StateT[j].cumC_resist[i];
			TimeSeries[n].incI_resist[i] += (double)StateT[j].cumI_resist[i];
			TimeSeries[n].cumT_resist[i] += (double)StateT[j].cumT_resist[i];
			StateT[j].cumC_resist[i] = StateT[j].cumI_resist[i] = 0;
		}
	}
	for (i = 0; i < 2; i++)
	{
		TimeSeries[n].incC_keyworker[i] = TimeSeries[n].incI_keyworker[i] = TimeSeries[n].cumT_keyworker[i] = TimeSeries[n].incD_keyworker[i] = 0;
		for (j = 0; j < P.NumThreads; j++)
		{
			TimeSeries[n].incC_keyworker[i] += (double)StateT[j].cumC_keyworker[i];
			TimeSeries[n].incI_keyworker[i] += (double)StateT[j].cumI_keyworker[i];
			TimeSeries[n].cumT_keyworker[i] += (double)StateT[j].cumT_keyworker[i];
			TimeSeries[n].incD_keyworker[i] += (double)StateT[j].cumD_keyworker[i];
			StateT[j].cumC_keyworker[i] = StateT[j].cumI_keyworker[i] = StateT[j].cumD_keyworker[i] = 0;
		}
	}

	for (i = 0; i < INFECT_TYPE_MASK; i++)
	{
		TimeSeries[n].incItype[i] = 0;
		for (j = 0; j < P.NumThreads; j++)
		{
			TimeSeries[n].incItype[i] += (double)StateT[j].cumItype[i];
			StateT[j].cumItype[i] = 0;
		}
	}
	if (P.DoAdUnits)
		for (i = 0; i <= P.NumAdunits; i++)
		{
			TimeSeries[n].incI_adunit[i] = TimeSeries[n].incC_adunit[i] = TimeSeries[n].incDC_adunit[i] = TimeSeries[n].incD_adunit[i] = TimeSeries[n].incDD_adunit[i] = TimeSeries[n].incDR_adunit[i] = TimeSeries[n].incSDB_adunit[i] = TimeSeries[n].incV_adunit[i] = TimeSeries[n].incVG_adunit[i] = TimeSeries[n].incETU_adunit[i] = TimeSeries[n].incH_adunit[i] = TimeSeries[n].cumT_adunit[i] = TimeSeries[n].incCT_adunit[i] = TimeSeries[n].incCC_adunit[i] = TimeSeries[n].capETU_adunit[i] = TimeSeries[n].ETU_adunit[i] = 0; //added detected cases: ggilani 03/02/15
			for (j = 0; j < P.NumThreads; j++)
			{
				TimeSeries[n].incI_adunit[i] += (double)StateT[j].cumI_adunit[i];
				TimeSeries[n].incC_adunit[i] += (double)StateT[j].cumC_adunit[i];
				TimeSeries[n].incDC_adunit[i] += (double)StateT[j].cumDC_adunit[i]; //added detected cases: ggilani 03/02/15
				TimeSeries[n].incD_adunit[i] += (double)StateT[j].cumD_adunit[i];
				TimeSeries[n].incDD_adunit[i] += (double)StateT[j].cumDD_adunit[i]; //added detected deaths: ggilani 03/02/15
				TimeSeries[n].incDR_adunit[i] += (double)StateT[j].cumDR_adunit[i]; //added detected recoveries: ggilani 03/02/15
				TimeSeries[n].incSDB_adunit[i] += (double)StateT[j].cumSDB_adunit[i]; //added safe burials: ggilani 05/10/23
				if ((t >= P.FuneralControlTimeStart) && (AdUnits[i].contactTraceThresholdCrossed))
				{
					TimeSeries[n].capSDB_adunit[i] = AdUnits[i].maxSDB;
				}
				else
				{
					TimeSeries[n].capSDB_adunit[i] = 0;
				}
				if ((t >= P.ContactTracingTimeStart) && (AdUnits[i].contactTraceThresholdCrossed))
				{
					TimeSeries[n].capCT_adunit[i] = AdUnits[i].contactTraceCapacity;
				}
				else
				{
					TimeSeries[n].capCT_adunit[i] = 0;
				}

				TimeSeries[n].incV_adunit[i] += (double)StateT[j].cumV_adunit[i]; //added vaccination: ggilani 05/10/23
				TimeSeries[n].incVG_adunit[i] += (double)StateT[j].cumVG_adunit[i]; //added vaccination: ggilani 05/10/23
				TimeSeries[n].incETU_adunit[i] += (double)StateT[j].cumETU_adunit[i]; //added hospitalisation
				TimeSeries[n].incH_adunit[i] += (double)StateT[j].cumH_adunit[i];
				TimeSeries[n].incCT_adunit[i] += (double)StateT[j].cumCT_adunit[i]; //added contact tracing: ggilani 15/06/17
				TimeSeries[n].incCC_adunit[i] += (double)StateT[j].cumCC_adunit[i]; //added cases who are contacts: ggilani 28/05/2019
				TimeSeries[n].cumT_adunit[i] += (double)StateT[j].cumT_adunit[i];
				StateT[j].cumI_adunit[i] = StateT[j].cumC_adunit[i] = StateT[j].cumH_adunit[i] = StateT[j].cumDC_adunit[i] = StateT[j].cumD_adunit[i] = StateT[j].cumDD_adunit[i] = StateT[j].cumDR_adunit[i] = StateT[j].cumSDB_adunit[i] = StateT[j].cumCT_adunit[i] = StateT[j].cumV_adunit[i] = StateT[j].cumVG_adunit[i] = StateT[j].cumETU_adunit[i] = StateT[j].cumCC_adunit[i] = 0; //added hospitalisation, detected cases, contact tracing: ggilani 03/02/15, 15/06/17
			}

			if ((P.DoHospitalisation) && (P.DoETUByAdUnit)) //added this to print out total number of beds in use at each time point: ggilani 31/10/14
			{
				for (j = 0; j < P.NumThreads; j++)
				{
					TimeSeries[n].ETU_adunit[i] += StateT[j].ETU_adunit[i];
					//TimeSeries[n].H_adunit[i] += StateT[j].H_adunit[i];

					if (StateT[j].capETU_adunit[i])
					{
						TimeSeries[n].capETU_adunit[i] = 1;
						StateT[j].capETU_adunit[i] = 0;
					}
				}

				TimeSeries[n].nBeds = State.NumBeds;
				TimeSeries[n].nBeds_adunit[i] = State.NumBeds_adunits[i];

			}
		}

	if (P.DoContactTracing) //added this to print out total number of beds in use at each time point: ggilani 31/10/14
	{
		for (i = 0; i < P.NumAdunits; i++)
		{
			TimeSeries[n].CT_adunit[i] = (double)AdUnits[i].nct;
		}
	}
	if ((P.DoPlaces) && (t >= P.PlaceCloseTimeStart))
	{
		for (i = 0; i < NUM_PLACE_TYPES; i++)
		{
			numPC = 0;
			for (j = 0; j < P.Nplace[i]; j++)
				if (PLACE_CLOSED(i, j)) numPC++;
			State.NumPlacesClosed[i] = numPC;
			TimeSeries[n].PropPlacesClosed[i] = ((double)numPC) / ((double)P.Nplace[i]);
		}
	}
	//if (State.cumC > P.NumUndetectedInfPreOutbreakAlert)
	//{
	//	P.ControlPropCasesId = P.PostAlertControlPropCasesId;
	//}
	if (State.cumDC >= P.PreControlClusterIdCaseThreshold)
	{
		if (P.OutbreakDetected == 0)
		{
			P.OutbreakDetected = 1; //mart outbreak as detected
			P.PropHospSeek *= P.RelChangeHospSeekPostOutbreak;
		}
		if (P.DoGlobalTriggers)
		{
			if (P.DoPerCapitaTriggers)
				D = (int)floor(((double)State.trigDC) * P.GlobalIncThreshPop / ((double)P.N));
			else
				D = State.trigDC;
			if (D >= P.TreatCellIncThresh)
			{
				if (P.TreatTimeStart >= 1e10) P.TreatTimeStart = t + P.TreatTimeStartBase;
				if (P.CaseIsolationTimeStart >= 1e10) P.CaseIsolationTimeStart = t + P.CaseIsolationTimeStartBase;
				if (P.HQuarantineTimeStart >= 1e10) P.HQuarantineTimeStart = t + P.HQuarantineTimeStartBase;
			}
			if (D >= P.VaccCellIncThresh)
			{
				if (P.VaccTimeStart >= 1e10)
				{
					P.VaccTimeStart = t + P.VaccTimeStartBase;
				}
				if (P.VaccNewCoursesStartTime >= 1e10)
				{
					P.VaccNewCoursesStartTime = t + P.VaccNewCoursesStartTimeBase;
					P.VaccNewCoursesBoostStartTime = t + P.VaccNewCoursesBoostStartTimeBase;
				}
			}
			if (D >= P.SocDistCellIncThresh)
			{
				if (P.SocDistTimeStart >= 1e10) P.SocDistTimeStart = t + P.SocDistTimeStartBase;
			}
			if (D >= P.PlaceCloseCellIncThresh)
			{
				if (P.PlaceCloseTimeStart >= 1e10) P.PlaceCloseTimeStart = t + P.PlaceCloseTimeStartBase;
				if ((P.PlaceCloseTimeStart2 >= 1e10) && (t >= P.PlaceCloseDuration + P.PlaceCloseTimeStart))
				{
					P.PlaceCloseTimeStart = t + P.PlaceCloseTimeStartBase2 - P.PlaceCloseTimeStartBase;
					P.PlaceCloseDuration = P.PlaceCloseDuration2;
				}
			}
			if (D >= P.MoveRestrCellIncThresh)
			{
				if (P.MoveRestrTimeStart >= 1e10) P.MoveRestrTimeStart = t + P.MoveRestrTimeStartBase;
			}
			if (D >= P.KeyWorkerProphCellIncThresh)
			{
				if (P.KeyWorkerProphTimeStart >= 1e10) P.KeyWorkerProphTimeStart = t + P.KeyWorkerProphTimeStartBase;
			}
			/*if (D >= P.GeoVaccCellIncThresh)
			{
				if (P.GeoVaccTimeStart >= 1e10) P.GeoVaccTimeStart = t + P.GeoVaccTimeStartBase;
			}*/
			//added this to enable funeral control start time when not doing funeral controls on an admin unit basis
			if ((D >= P.FuneralControlCellIncThresh))
			{
				if (P.FuneralControlTimeStart >= 1e10)
				{
					P.FuneralControlTimeStart = t + P.FuneralControlTimeStartBase;
				}
			}

			if (D >= P.ContactTracingCellIncThresh)
			{
				if (P.ContactTracingTimeStart >= 1e10)
				{
					P.ContactTracingTimeStart = t + P.ContactTracingTimeStartBase;
				}
			}
			/*if (D >= P.RingVaccCellIncThresh)
			{
				if (P.RingVaccTimeStart >= 1e10)
				{
					P.RingVaccTimeStart = t + P.RingVaccTimeStartBase;
				}
			}*/
			if (D >= P.ETUCellIncThresh)
			{
				if (P.ETUTimeStart >= 1e10)
				{
					P.ETUTimeStart = t + P.ETUTimeStartBase;
				}
			}
		}
		else
		{
			if (P.TreatTimeStart >= 1e10) P.TreatTimeStart = t + P.TreatTimeStartBase;
			if (P.CaseIsolationTimeStart >= 1e10) P.CaseIsolationTimeStart = t + P.CaseIsolationTimeStartBase;
			if (P.HQuarantineTimeStart >= 1e10) P.HQuarantineTimeStart = t + P.HQuarantineTimeStartBase;
			if (P.VaccTimeStart >= 1e10)
			{
				P.VaccTimeStart = t + P.VaccTimeStartBase;
				//fprintf(stderr, "t=%lg, P.VaccTimeStart=%lg, P.VaccTimeStartBase=%lg\n", t, P.VaccTimeStart, P.VaccTimeStartBase);
			}
			if (P.VaccNewCoursesStartTime >= 1e10)
			{
				P.VaccNewCoursesStartTime = t + P.VaccNewCoursesStartTimeBase;
				P.VaccNewCoursesBoostStartTime = t + P.VaccNewCoursesBoostStartTimeBase;
			}
			if (P.SocDistTimeStart >= 1e10) P.SocDistTimeStart = t + P.SocDistTimeStartBase;
			if (P.PlaceCloseTimeStart >= 1e10) P.PlaceCloseTimeStart = t + P.PlaceCloseTimeStartBase;
			if (P.MoveRestrTimeStart >= 1e10) P.MoveRestrTimeStart = t + P.MoveRestrTimeStartBase;
			if (P.KeyWorkerProphTimeStart >= 1e10) P.KeyWorkerProphTimeStart = t + P.KeyWorkerProphTimeStartBase;
			if (P.FuneralControlTimeStart >= 1e10) P.FuneralControlTimeStart = t + P.FuneralControlTimeStartBase;
			if (P.ContactTracingTimeStart >= 1e10) P.ContactTracingTimeStart = t + P.ContactTracingTimeStartBase;
			//if (P.GeoVaccTimeStart >= 1e10) P.GeoVaccTimeStart = t + P.GeoVaccTimeStartBase;
			//if (P.RingVaccTimeStart >= 1e10) P.RingVaccTimeStart = t + P.RingVaccTimeStartBase;
			if (P.ETUTimeStart >= 1e10) P.ETUTimeStart = t + P.ETUTimeStartBase;
			//added this to correctly set funeral control start times if everything is happening on the admin unit level - ggilani 07/03/17
			//if(P.DoFuneralByAdUnit)
			//{
			//	for(i=0;i<P.NumAdunits;i++)
			//	{
			//		//if we haven't passed the threshold yet and started funeral controls, check to see if we should
			//		if(AdUnits[i].startFuneralControl>=1e10)
			//		{
			//			//reset nCase to zero
			//		nCase=0;
			//			//add up cases in the admin unit so far...
			//			for(j=0;j<=n;j++)
			//			{
			//				nCase=nCase+TimeSeries[j].incDC_adunit[i];
			//			}
			//			//if number of cases exceeds the admin unit trigger, set the start of funeral control time
			//			if(nCase>AdUnits[i].caseDetPreFuneralControl)
			//			{
			//				AdUnits[i].startFuneralControl=t+AdUnits[i].timeToSafeFuneral;
			//			}
			//		}
			//	}
			//}
		}

		if (P.AirportCloseTimeStart >= 1e10) P.AirportCloseTimeStart = t + P.AirportCloseTimeStartBase;
	}
	if ((P.PlaceCloseIndepThresh > 0) && (((double)State.cumDC) >= P.PlaceCloseIndepThresh))
	{
		if (P.PlaceCloseTimeStart >= 1e10) P.PlaceCloseTimeStart = t + P.PlaceCloseTimeStartBase;
	}
	if (P.OutputBitmap >= 1)
	{
		TSMean = TSMeanNE; TSVar = TSVarNE;
		CaptureBitmap(n, 0);
		OutputBitmap(t, 0);
	}

}

void RecordInfTypes(void)
{
	int i, j, k, l, lc, lc2, b, c, n, nf, i2;
	double* res, * res_av, * res_var, t, s;

	for (n = 0; n < P.NumSamples; n++)
	{
		for (i = 0; i < INFECT_TYPE_MASK; i++) TimeSeries[n].Rtype[i] = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++) TimeSeries[n].Rage[i] = 0;
		TimeSeries[n].Rdenom = 0;
	}
	for (i = 0; i < INFECT_TYPE_MASK; i++) inftype[i] = 0;
	for (i = 0; i < MAX_COUNTRIES; i++) infcountry[i] = 0;
	for (i = 0; i < MAX_SEC_REC; i++)
		for (j = 0; j < MAX_GEN_REC; j++)
			indivR0[i][j] = 0;
	for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
		for (j = 0; j <= MAX_HOUSEHOLD_SIZE; j++)
			inf_household[i][j] = case_household[i][j] = 0;
	for (b = 0; b < P.NC; b++)
		if ((Cells[b].S != Cells[b].n) || (Cells[b].R > 0))
			for (c = 0; c < Cells[b].n; c++)
				Hosts[Cells[b].members[c]].listpos = 0;
	//	for(b=0;b<P.NC;b++)
	//		if((Cells[b].S!=Cells[b].n)||(Cells[b].R>0))
	{
		j = k = l = lc = lc2 = 0; t = 1e10;
		//			for(c=0;c<Cells[b].n;c++)
		for (i = 0; i < P.N; i++)
		{
			//				i=Cells[b].members[c];
			if (j == 0) j = k = Households[Hosts[i].hh].nh;
			if ((Hosts[i].inf != 0) && (Hosts[i].inf != 4))
			{
				if (Hosts[i].latent_time * P.TimeStep <= P.SampleTime)
					TimeSeries[(int)(Hosts[i].latent_time * P.TimeStep / P.SampleStep)].Rdenom++;
				infcountry[Mcells[Hosts[i].mcell].country]++;
				if (abs(Hosts[i].inf) < 3)
					l = -1;
				else if (l >= 0)
					l++;
				if ((l >= 0) && ((Hosts[i].inf == -3) || (Hosts[i].inf == -5)))
				{
					lc2++;
					if (Hosts[i].latent_time * P.TimeStep <= t) // This convoluted logic is to pick up households where the index is symptomatic
					{
						lc = 1; t = Hosts[i].latent_time * P.TimeStep;
					}
				}
				else if ((l > 0) && (Hosts[i].latent_time * P.TimeStep < t))
				{
					lc = 0; t = Hosts[i].latent_time * P.TimeStep;
				}
				i2 = Hosts[i].infector;
				if (i2 >= 0)
				{
					Hosts[i2].listpos++;
					if (Hosts[i2].latent_time * P.TimeStep <= P.SampleTime)
					{
						TimeSeries[(int)(Hosts[i2].latent_time * P.TimeStep / P.SampleStep)].Rtype[Hosts[i].infect_type % INFECT_TYPE_MASK]++;
						TimeSeries[(int)(Hosts[i2].latent_time * P.TimeStep / P.SampleStep)].Rage[HOST_AGE_GROUP(i)]++;
					}
				}
			}
			inftype[Hosts[i].infect_type % INFECT_TYPE_MASK]++;
			j--;
			if (j == 0)
			{
				if (l < 0) l = 0;
				inf_household[k][l]++;
				case_household[k][lc2]++; //now recording total symptomatic cases, rather than infections conditional on symptomatic index
				l = lc = lc2 = 0; t = 1e10;
			}
		}
	}
	for (b = 0; b < P.NC; b++)
		if ((Cells[b].S != Cells[b].n) || (Cells[b].R > 0))
			for (c = 0; c < Cells[b].n; c++)
			{
				i = Cells[b].members[c];
				if ((abs(Hosts[i].inf) == 3) || (abs(Hosts[i].inf) == 5))
				{
					l = Hosts[i].infect_type / INFECT_TYPE_MASK;
					if ((l < MAX_GEN_REC) && (Hosts[i].listpos < MAX_SEC_REC)) indivR0[Hosts[i].listpos][l]++;
				}
			}
	/* 	if(!TimeSeries[P.NumSamples-1].extinct) */
	{
		for (i = 0; i < INFECT_TYPE_MASK; i++) inftype_av[i] += inftype[i];
		for (i = 0; i < MAX_COUNTRIES; i++)
		{
			infcountry_av[i] += infcountry[i];
			if (infcountry[i] > 0) infcountry_num[i]++;
		}
		for (i = 0; i < MAX_SEC_REC; i++)
			for (j = 0; j < MAX_GEN_REC; j++)
				indivR0_av[i][j] += indivR0[i][j];
		for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
			for (j = 0; j <= MAX_HOUSEHOLD_SIZE; j++)
			{
				inf_household_av[i][j] += inf_household[i][j];
				case_household_av[i][j] += case_household[i][j];
			}
	}
	for (n = 0; n < P.NumSamples; n++)
	{
		s = 0;
		if (TimeSeries[n].Rdenom == 0) TimeSeries[n].Rdenom = 1e-10;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			TimeSeries[n].Rage[i] /= TimeSeries[n].Rdenom;
		for (i = 0; i < INFECT_TYPE_MASK; i++)
			s += (TimeSeries[n].Rtype[i] /= TimeSeries[n].Rdenom);
		TimeSeries[n].Rdenom = s;
	}
	nf = (sizeof(results) - 3 * sizeof(float*)) / sizeof(double);
	if (!P.DoAdUnits) nf -= MAX_ADUNITS;
	fprintf(stderr, "extinct=%lg (%i)\n", TimeSeries[P.NumSamples - 1].extinct, P.NumSamples - 1);
	if (TimeSeries[P.NumSamples - 1].extinct)
	{
		TSMean = TSMeanE; TSVar = TSVarE; P.NRactE++;
	}
	else
	{
		TSMean = TSMeanNE; TSVar = TSVarNE; P.NRactNE++;
	}
	s = 0;
	for (n = 0; n < P.NumSamples; n++)
	{
		if (s < TimeSeries[n].incC) { s = TimeSeries[n].incC; t = P.SampleStep * ((double)n); }
		res = (double*)&TimeSeries[n];
		res_av = (double*)&TSMean[n];
		res_var = (double*)&TSVar[n];
		for (i = 0; i < nf; i++)
		{
			res_av[i] += res[i];
			res_var[i] += res[i] * res[i];
		}
		if (TSMean[n].cumTmax < TimeSeries[n].cumT) TSMean[n].cumTmax = TimeSeries[n].cumT;
		if (TSMean[n].cumVmax < TimeSeries[n].cumV) TSMean[n].cumVmax = TimeSeries[n].cumV;
	}
	PeakHeightSum += s;
	PeakHeightSS += s * s;
	PeakTimeSum += t;
	PeakTimeSS += t * t;
}

/* Function: RecordEvent(t, ai)
 * Records an infection event in the event log
 *
 * Parameters:
 *	t: time of infection event
 *	ai: index of infectee
 *	nEventsPoint: pointer to number of events
 *
 * Returns: void
 *
 * Author: ggilani, Date: 10/10/2014
 */
void RecordEvent(double t, int ai, int run, int type, int tn) //added int as argument to RecordEvent to record run number: ggilani - 15/10/14
{
	//Declare int to store infector's index
	int bi;

	bi = Hosts[ai].infector;

	//Save information to event
#pragma omp critical (inf_event)
	{
		InfEventLog[*nEvents].run = run;
		InfEventLog[*nEvents].type = type;
		InfEventLog[*nEvents].t = t;
		InfEventLog[*nEvents].infectee_ind = ai;
		InfEventLog[*nEvents].infectee_adunit = Mcells[Hosts[ai].mcell].adunit;
		InfEventLog[*nEvents].infectee_x = Households[Hosts[ai].hh].loc_x + P.SpatialBoundingBox[0];
		InfEventLog[*nEvents].infectee_y = Households[Hosts[ai].hh].loc_y + P.SpatialBoundingBox[1];
		InfEventLog[*nEvents].listpos = Hosts[ai].listpos;
		InfEventLog[*nEvents].infectee_cell = Hosts[ai].pcell;
		InfEventLog[*nEvents].thread = tn;
		if (type == 0) //infection event - record time of onset of infector and infector
		{
			InfEventLog[*nEvents].infector_ind = bi;
			if (bi < 0)
			{
				InfEventLog[*nEvents].t_infector = -1;
				InfEventLog[*nEvents].infector_cell = -1;
			}
			else
			{
				InfEventLog[*nEvents].t_infector = (int)(Hosts[bi].infection_time / P.TimeStepsPerDay);
				InfEventLog[*nEvents].infector_cell = Hosts[bi].pcell;
			}
		}
		else if (type == 1) //onset event - record infectee's onset time
		{
			InfEventLog[*nEvents].t_infector = (int)(Hosts[ai].infection_time / P.TimeStepsPerDay);
		}
		else if ((type == 2) || (type == 3)) //recovery or death event - record infectee's onset time
		{
			InfEventLog[*nEvents].t_infector = (int)(Hosts[ai].latent_time / P.TimeStepsPerDay);
		}

		//increment the index of the infection event
		(*nEvents)++;
	}

}


/** function: CalcOriginDestMatrix_adunit()
 *
 * purpose: to output the origin destination matrix between admin units
 *
 * parameters: none
 *
 * returns: none
 *
 * author: ggilani, date: 28/01/15
 */
void CalcOriginDestMatrix_adunit()
{
	int tn, i, j, k, l, m, p;
	double total_flow, flow;
	int cl_from, cl_to, cl_from_mcl, cl_to_mcl, mcl_from, mcl_to;
	double pop_dens_from[MAX_ADUNITS], pop_dens_to[MAX_ADUNITS];

#pragma omp parallel for private(tn,i,j,k,l,m,p,total_flow,mcl_from,mcl_to,cl_from,cl_to,cl_from_mcl,cl_to_mcl,pop_dens_from,pop_dens_to,flow) schedule(static) //reduction(+:s,t2)
	for (tn = 0; tn < P.NumThreads; tn++)
	{
		for (i = tn; i < P.NCP; i += P.NumThreads)
		{
			//reset pop density matrix to zero
			for (k = 0; k < P.NumAdunits; k++)
			{
				pop_dens_from[k] = 0.0;
			}

			//find index of cell from which flow travels
			cl_from = CellLookup[i] - Cells;
			cl_from_mcl = (cl_from / P.nch) * P.NMCL * P.nmch + (cl_from % P.nch) * P.NMCL;

			//loop over microcells in these cells to find populations in each admin unit and so flows
			for (k = 0; k < P.NMCL; k++)
			{
				for (l = 0; l < P.NMCL; l++)
				{
					//get index of microcell
					mcl_from = cl_from_mcl + l + k * P.nmch;
					if (Mcells[mcl_from].n > 0)
					{
						//get proportion of each population of cell that exists in each admin unit
						pop_dens_from[Mcells[mcl_from].adunit] += (((double)Mcells[mcl_from].n) / ((double)Cells[cl_from].n));
					}
				}
			}

			for (j = i; j < P.NCP; j++)
			{
				//reset pop density matrix to zero
				for (m = 0; m < P.NumAdunits; m++)
				{
					pop_dens_to[m] = 0.0;
				}

				//find index of cell which flow travels to
				cl_to = CellLookup[j] - Cells;
				cl_to_mcl = (cl_to / P.nch) * P.NMCL * P.nmch + (cl_to % P.nch) * P.NMCL;
				//calculate distance and kernel between the cells
				//total_flow=Cells[cl_from].max_trans[j]*Cells[cl_from].n*Cells[cl_to].n;
				if (j == 0)
				{
					total_flow = Cells[cl_from].cum_trans[j] * Cells[cl_from].n;
				}
				else
				{
					total_flow = (Cells[cl_from].cum_trans[j] - Cells[cl_from].cum_trans[j - 1]) * Cells[cl_from].n;
				}

				//loop over microcells within destination cell
				for (m = 0; m < P.NMCL; m++)
				{
					for (p = 0; p < P.NMCL; p++)
					{
						//get index of microcell
						mcl_to = cl_to_mcl + p + m * P.nmch;
						if (Mcells[mcl_to].n > 0)
						{
							//get proportion of each population of cell that exists in each admin unit
							pop_dens_to[Mcells[mcl_to].adunit] += (((double)Mcells[mcl_to].n) / ((double)Cells[cl_to].n));
						}
					}
				}

				for (m = 0; m < P.NumAdunits; m++)
				{
					for (p = 0; p < P.NumAdunits; p++)
					{
						if (m != p)
						{
							if (AdUnits[m].id / P.CountryDivisor == AdUnits[p].id / P.CountryDivisor)
							{
								flow = total_flow * pop_dens_from[m] * pop_dens_to[p];
							}
							else
							{
								flow = total_flow * pop_dens_from[m] * pop_dens_to[p] * P.PropCrossBorderInf;
							}
							StateT[tn].origin_dest[m][p] += flow;
							StateT[tn].origin_dest[p][m] += flow;
						}
					}
				}
			}

			////loop over microcells within cell to find the proportion of the cell population in each admin unit
			//k=(cl_from/P.nch)*P.NMCL*P.nmch+(cl_from%P.nch)*P.NMCL;
			//for(l=0;l<P.NMCL;l++)
			//{
			//	for(m=0;m<P.NMCL;m++)
			//	{
			//		mcl_from=k+m+l*P.nmch;
			//		pop_cell_from[Mcells[mcl_from].adunit]+=Mcells[mcl_from].n;
			//	}
			//}
			////loop over cells
			//for(p=(i+1);p<P.NCP;p++)
			//{
			//	//reset population array
			//	for(j=0;j<P.NumAdunits;j++)
			//	{
			//		pop_cell_to[j]=0.0;
			//	}
			//	cl_to=CellLookup[p]-Cells;
			//	//loop over microcells within cell to find the proportion of the cell population in each admin unit
			//	q=(cl_to/P.nch)*P.NMCL*P.nmch+(cl_to%P.nch)*P.NMCL;
			//	for(l=0;l<P.NMCL;l++)
			//	{
			//		for(m=0;m<P.NMCL;m++)
			//		{
			//			mcl_to=q+m+l*P.nmch;
			//			pop_cell_to[Mcells[mcl_to].adunit]+=Mcells[mcl_to].n;
			//		}
			//	}

			//	//find distance and kernel function between cells
			//	dist=dist2_cc_min(Cells+cl_from,Cells+cl_to);
			//	dist_kernel=numKernel(dist);

			//	//add flow between adunits based on how population is distributed
			//	for(l=0;l<P.NumAdunits;l++)
			//	{
			//		for(m=(l+1);m<P.NumAdunits;m++)
			//		{
			//			AdUnits[l].origin_dest[m]+=pop_cell_from[l]*pop_cell_to[m]*dist_kernel;
			//			AdUnits[m].origin_dest[l]+=pop_cell_from[l]*pop_cell_to[m]*dist_kernel;
			//		}
			//	}

		}
	}

	//Sum up flow between adunits across threads
	for (i = 0; i < P.NumAdunits; i++)
	{
		for (j = 0; j < P.NumAdunits; j++)
		{
			for (k = 0; k < P.NumThreads; k++)
			{
				AdUnits[i].origin_dest[j] += StateT[k].origin_dest[i][j];
			}
		}
	}

}

/** function: CalcOriginDestMatrix_adunit()
 *
 * purpose: to output the origin destination matrix between admin units
 *
 * parameters: none
 *
 * returns: none
 *
 * author: ggilani, date: 28/01/15
 */
 //void CalcOriginDestMatrix_adunit()
 //{
 //	int tn,i,j;
 //	double flow;
 //	int mcl_from,mcl_to;
 //
 //#pragma omp parallel for private(i,j,flow,mcl_from,mcl_to) schedule(static,1000) //reduction(+:s,t2)
 //	for(i=0;i<P.NMCP;i++)
 //	{
 //		//loop over microcells
 //		for(j=(i+1);j<P.NMCP;j++)
 //		{
 //			mcl_from=McellLookup[i]-Mcells;
 //			mcl_to=McellLookup[j]-Mcells;
 //
 //			if(Mcells[mcl_from].adunit!=Mcells[mcl_to].adunit)
 //			{
 //				flow=numKernel(dist2_mm_min(Mcells+mcl_from,Mcells+mcl_to))*Mcells[mcl_from].n*Mcells[mcl_to].n;
 //#pragma omp atomic
 //				AdUnits[Mcells[mcl_from].adunit].origin_dest[Mcells[mcl_to].adunit]+=flow;
 //#pragma omp atomic
 //				AdUnits[Mcells[mcl_to].adunit].origin_dest[Mcells[mcl_from].adunit]+=flow;
 //			}
 //
 //
 //
 //		}
 //	}
 //
 //}



void SaveAgeDistrib(void)
{
	int i;
	FILE* dat;
	char outname[1024];

	sprintf(outname, "%s.agedist.xls", OutFile);
	if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
	if (P.DoDeath)
	{
		fprintf(dat, "age\tfreq\tlifeexpect\n");
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, "%i\t%lg\t%lg\n", i, AgeDist[i], AgeDist2[i]);
		fprintf(dat, "\np\tlife_expec\tage\n");
		for (i = 0; i <= 1000; i++)
			fprintf(dat, "%lg\t%lg\t%i\n", ((double)i) / 1000, P.InvLifeExpecDist[0][i], State.InvAgeDist[0][i]);
	}
	else
	{
		fprintf(dat, "age\tfreq\n");
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, "%i\t%lg\n", i, AgeDist[i]);
	}

	fclose(dat);
}

void SaveDistribs(void)
{
	int i, j, k;
	FILE* dat;
	char outname[1024];
	double s;

	if (P.DoPlaces)
	{
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != HOTEL_PLACE_TYPE)
			{
				for (i = 0; i < P.Nplace[j]; i++)
					Places[j][i].n = 0;
				for (i = 0; i < P.N; i++)
				{
					if (Hosts[i].PlaceLinks[j] >= P.Nplace[j])
						fprintf(stderr, "*%i %i: %i %i", i, j, Hosts[i].PlaceLinks[j], P.Nplace[j]);
					else if (Hosts[i].PlaceLinks[j] >= 0)
						Places[j][Hosts[i].PlaceLinks[j]].n++;
				}
			}
		for (j = 0; j < P.PlaceTypeNum; j++)
			for (i = 0; i < MAX_DIST; i++)
				PlaceDistDistrib[j][i] = 0;
		for (i = 0; i < P.N; i++)
			for (j = 0; j < P.PlaceTypeNum; j++)
				if ((j != HOTEL_PLACE_TYPE) && (Hosts[i].PlaceLinks[j] >= 0))
					if (Hosts[i].PlaceLinks[j] >= P.Nplace[j])
						fprintf(stderr, "*%i %i: %i ", i, j, Hosts[i].PlaceLinks[j]);
					else if ((!P.DoOutputPlaceDistForOneAdunit) ||
						((AdUnits[Mcells[Hosts[i].mcell].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor == (P.OutputPlaceDistAdunit % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor))
					{
						k = Hosts[i].PlaceLinks[j];
						s = sqrt(dist2_raw(Households[Hosts[i].hh].loc_x, Households[Hosts[i].hh].loc_y, Places[j][k].loc_x, Places[j][k].loc_y)) / OUTPUT_DIST_SCALE;
						k = (int)s;
						if (k < MAX_DIST) PlaceDistDistrib[j][k]++;
					}
		for (j = 0; j < P.PlaceTypeNum; j++)
			for (i = 0; i < MAX_PLACE_SIZE; i++)
				PlaceSizeDistrib[j][i] = 0;
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != HOTEL_PLACE_TYPE)
				for (i = 0; i < P.Nplace[j]; i++)
					if (Places[j][i].n < MAX_PLACE_SIZE)
						PlaceSizeDistrib[j][Places[j][i].n]++;
		sprintf(outname, "%s.placedist.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "dist");
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != HOTEL_PLACE_TYPE)
				fprintf(dat, "\tfreq_p%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < MAX_DIST; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 0; j < P.PlaceTypeNum; j++)
				if (j != HOTEL_PLACE_TYPE)
					fprintf(dat, "\t%i", PlaceDistDistrib[j][i]);
			fprintf(dat, "\n");
		}
		fclose(dat);
		sprintf(outname, "%s.placesize.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "size");
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != HOTEL_PLACE_TYPE)
				fprintf(dat, "\tfreq_p%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < MAX_PLACE_SIZE; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 0; j < P.PlaceTypeNum; j++)
				if (j != HOTEL_PLACE_TYPE)
					fprintf(dat, "\t%i", PlaceSizeDistrib[j][i]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	if ((P.DoAdunitBoundaries) && (P.PrivateTreatPropCases > 0) && (P.DoHouseholds))
	{
		sprintf(outname, "%s.privatestock.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "adunit_code\thouseholds\thouseholds_with_stock\n");
		for (j = 0; j < P.NumAdunits; j++)
			fprintf(dat, "%i\t%i\t%i\n", AdUnits[j].id, P.HouseholdsByAdunit[j], P.PrivateStockByAdunit[j]);
		fclose(dat);
	}
}

/** function: SaveOriginDestMatrix
 *
 * purpose: to save the calculated origin destination matrix to file
 * parameters: none
 * returns: none
 *
 * author: ggilani, 13/02/15
 */

void SaveOriginDestMatrix(void)
{
	int i, j;
	FILE* dat;
	char outname[1024];

	sprintf(outname, "%s.origdestmat.csv", OutFile);
	if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat, "0,");
	for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "%i,", (AdUnits[i].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor);
	fprintf(dat, "\n");
	for (i = 0; i < P.NumAdunits; i++)
	{
		fprintf(dat, "%i,", (AdUnits[i].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor);
		for (j = 0; j < P.NumAdunits; j++)
		{
			fprintf(dat, "%lg,", AdUnits[i].origin_dest[j]);
		}
		fprintf(dat, "\n");
	}
	fclose(dat);
}

void SaveResults(void)
{
	int i, j, k;
	double t;
	FILE* dat;
	char outname[1024];

	sprintf(outname, "%s.csv", OutFile);
	if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat, "t,S,L,I,R,D,incI,incR,incFC,incFI,incC,incDC,incD,incDD,incSDB,incTC,incETU,incH,incCT,incCC,cumT,cumTP,cumV,capV,cumVG,capVG,nBeds,Extinct,Detected,rmsRad,maxRad\n");//\t\t%lg\t%lg\t%lg\n",P.R0household,P.R0places,P.R0spatial);
	for (i = 0; i < P.NumSamples; i++)
	{
		fprintf(dat, "%lg,%lf,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
			TimeSeries[i].t, TimeSeries[i].S, TimeSeries[i].L, TimeSeries[i].I,
			TimeSeries[i].R, TimeSeries[i].D, TimeSeries[i].incI,
			TimeSeries[i].incR, TimeSeries[i].incFC, TimeSeries[i].incFI, TimeSeries[i].incC, TimeSeries[i].incDC, TimeSeries[i].incD, TimeSeries[i].incDD, TimeSeries[i].incSDB, TimeSeries[i].incTC, TimeSeries[i].incETU, TimeSeries[i].incH, TimeSeries[i].incCT, TimeSeries[i].incCC, //added incidence funeral transmissions and hospitalisation
			TimeSeries[i].cumT, TimeSeries[i].cumTP, TimeSeries[i].cumV, TimeSeries[i].capV, TimeSeries[i].cumVG, TimeSeries[i].capVG, TimeSeries[i].nBeds, TimeSeries[i].extinct, TimeSeries[i].detected, TimeSeries[i].rmsRad, TimeSeries[i].maxRad);
	}
	fclose(dat);
	if (P.DoControlOutput)
	{
		sprintf(outname, "%s.controls.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t\tS\tincC\tincTC\tincFC\tincFI\tincETU\tincH\tcumT\tcumUT\tcumTP\tcumV\tincHQ\tincAC\tincAH\tincAA\tincACS\tincAPC\tincAPA\tincAPCS");
		for (j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\tprClosed_%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg\t%lf\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg",
				TimeSeries[i].t, TimeSeries[i].S, TimeSeries[i].incC, TimeSeries[i].incTC, TimeSeries[i].incFC, TimeSeries[i].incFI, TimeSeries[i].incETU, TimeSeries[i].incH, //added incidence of funeral transmissions and hospitalisation
				TimeSeries[i].cumT, TimeSeries[i].cumUT, TimeSeries[i].cumTP, TimeSeries[i].cumV, TimeSeries[i].incHQ,
				TimeSeries[i].incAC, TimeSeries[i].incAH, TimeSeries[i].incAA, TimeSeries[i].incACS,
				TimeSeries[i].incAPC, TimeSeries[i].incAPA, TimeSeries[i].incAPCS);
			for (j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\t%lg", TimeSeries[i].PropPlacesClosed[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	if (P.DoAgeOutput)
	{
		sprintf(outname, "%s.age.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",I%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",C%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",D%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",DC%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",ETU%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",H%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",V%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg", TimeSeries[i].t);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", TimeSeries[i].incIa[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", TimeSeries[i].incCa[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", TimeSeries[i].incDa[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", TimeSeries[i].incDCa[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", TimeSeries[i].incETUa[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", TimeSeries[i].incHa[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", TimeSeries[i].incVa[j]);
			fprintf(dat, "\n");
		}
		fprintf(dat, "dist");
		for (j = 0; j < NUM_AGE_GROUPS; j++)
			fprintf(dat, ",%lg", AgeDist[j]);
		fprintf(dat, "\n");
		fclose(dat);
	}
	if ((P.DoAdUnits) && (P.DoAdunitOutput))
	{
		sprintf(outname, "%s.adunit.csv", OutFile); //modifying to csv file
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t,");
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "I_%s,", AdUnits[i].ad_name); //"\tI%i"
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "C_%s,", AdUnits[i].ad_name); //"\tC%i"
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "DC_%s,", AdUnits[i].ad_name); //added detected cases: ggilani 03/02/15
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "D_%s,", AdUnits[i].ad_name); //added deaths: ggilani 05/10/23
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "DD_%s,", AdUnits[i].ad_name); //added detected deaths: ggilani 05/10/23
		//for(i=0;i<P.NumAdunits;i++) fprintf(dat,"T_%s,",AdUnits[i].ad_name); //"\tT%i"

		if (P.DoFuneralTransmission)
		{
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "SDB_%s,", AdUnits[i].ad_name); //added safe burials: ggilani 05/10/23
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "capSDB_%s,", AdUnits[i].ad_name); //added safe burials: ggilani 05/10/23
		}

		if ((P.DoHospitalisation) & (P.DoETUByAdUnit))
		{
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "Beds_%s,", AdUnits[i].ad_name); //"\tT%i" //added number of beds
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "incETU_%s,", AdUnits[i].ad_name); //"\tT%i" //added incidence of hospitalisation
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "ETU_%s,", AdUnits[i].ad_name); //"\tT%i" //added hospitalisation
			if (P.DoOutputETUCapacity)
			{
				for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "capETU_%s,", AdUnits[i].ad_name); //"\tT%i" //added hospital capacity indicator: ggilani 25/04/22
			}
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "incH_%s,", AdUnits[i].ad_name); //"\tT%i" //added incidence of hospitalisation
			//for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "H_%s,", AdUnits[i].ad_name); //"\tT%i" //added hospitalisation
		}
		if (P.DoContactTracing)
		{
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "incCT_%s,", AdUnits[i].ad_name); //"\tT%i" //added contact tracing
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "CT_%s,", AdUnits[i].ad_name); //"\tT%i" //added contact tracing
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "capCT_%s,", AdUnits[i].ad_name); //"\tT%i" //added contact tracing
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "incCC_%s,", AdUnits[i].ad_name); //"\tT%i" //added incidence of cases who are contacts

			//for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "CTC_%s,", AdUnits[i].ad_name); //"\tT%i" //added contact tracing

		}
		if (P.DoRingVaccination)
		{
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "incV_%s,", AdUnits[i].ad_name); //"\tT%i" //added vaccination
		}
		if (P.DoGeoVaccination)
		{
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "incVG_%s,", AdUnits[i].ad_name); //"\tT%i" //added vaccination
		}
		//for(i=0;i<P.NumAdunits;i++) fprintf(dat,"%lg,",P.PopByAdunit[i][0]); //"\t%lg"
		//for(i=0;i<P.NumAdunits;i++) fprintf(dat,"%lg,",P.PopByAdunit[i][1]); //"\t%lg"
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg,", TimeSeries[i].t); //"%lg"
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "%lg,", TimeSeries[i].incI_adunit[j]); //"\t%lg"
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "%lg,", TimeSeries[i].incC_adunit[j]); //"\t%lg"
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "%lg,", TimeSeries[i].incDC_adunit[j]); //added detected cases: ggilani 03/02/15 
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "%lg,", TimeSeries[i].incD_adunit[j]); //"\t%lg"
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "%lg,", TimeSeries[i].incDD_adunit[j]); //"\t%lg"

			if (P.DoFuneralTransmission)
			{
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].incSDB_adunit[j]); //added safe burials: ggilani 05/10/23
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].capSDB_adunit[j]); //added safe burials: ggilani 05/10/23
			}

			if ((P.DoHospitalisation) & (P.DoETUByAdUnit))
			{
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].nBeds_adunit[j]); //"\t%lg" //added number of beds
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].incETU_adunit[j]); //"\t%lg" //added incidence hospitalisation
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].ETU_adunit[j]); //"\t%lg" //added hospitalisation
				if (P.DoOutputETUCapacity)
				{
					for (j = 0; j < P.NumAdunits; j++)
						fprintf(dat, "%lg,", TimeSeries[i].capETU_adunit[j]); //"\t%lg" //added output for hospital capacity: ggilani 25/04/22
				}
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].incH_adunit[j]); //"\t%lg" //added incidence hospitalisation
				//for (j = 0; j < P.NumAdunits; j++)
					//fprintf(dat, "%lg,", TimeSeries[i].H_adunit[j]); //"\t%lg" //added hospitalisation
			}
			if (P.DoContactTracing)
			{
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].incCT_adunit[j]); //"\t%lg" //added contact tracing
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].CT_adunit[j]); //"\t%lg" //added contact tracing
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].capCT_adunit[j]); //"\t%lg" //added contact tracing
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].incCC_adunit[j]); //"\t%lg" //added cases who are contacts
				//for (j = 0; j < P.NumAdunits; j++)
				//	fprintf(dat, "%lg,", TimeSeries[i].CC_adunit[j]); //"\t%lg" //added cases who are contacts
			}
			if (P.DoRingVaccination)
			{
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].incV_adunit[j]); //"\t%lg" //added vaccination
			}
			if (P.DoGeoVaccination)
			{
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "%lg,", TimeSeries[i].incVG_adunit[j]); //"\t%lg" //added geographic vaccination
			}
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	if (P.EvolResistNumTypes > 1)
	{
		sprintf(outname, "%s.resist.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for (i = 0; i < P.EvolResistNumTypes; i++) fprintf(dat, "\tI%i", i);
		for (i = 0; i < P.EvolResistNumTypes; i++) fprintf(dat, "\tC%i", i);
		for (i = 0; i < P.EvolResistNumTypes; i++) fprintf(dat, "\tT%i", i);
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg", TimeSeries[i].t);
			for (j = 0; j < P.EvolResistNumTypes; j++)
				fprintf(dat, "\t%lg", TimeSeries[i].incI_resist[j]);
			for (j = 0; j < P.EvolResistNumTypes; j++)
				fprintf(dat, "\t%lg", TimeSeries[i].incC_resist[j]);
			for (j = 0; j < P.EvolResistNumTypes; j++)
				fprintf(dat, "\t%lg", TimeSeries[i].cumT_resist[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	//if(P.KeyWorkerProphTimeStartBase<P.SampleTime)
	if (P.DoKeyworkerOutput)
	{
		sprintf(outname, "%s.keyworker.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for (i = 0; i < 2; i++) fprintf(dat, ",I%i", i);
		for (i = 0; i < 2; i++) fprintf(dat, ",C%i", i);
		for (i = 0; i < 2; i++) fprintf(dat, ",T%i", i);
		for (i = 0; i < 2; i++) fprintf(dat, ",D%i", i);
		fprintf(dat, ",%i,%i\n", P.KeyWorkerNum, P.KeyWorkerIncHouseNum);
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg", TimeSeries[i].t);
			for (j = 0; j < 2; j++)
				fprintf(dat, ",%lg", TimeSeries[i].incI_keyworker[j]);
			for (j = 0; j < 2; j++)
				fprintf(dat, ",%lg", TimeSeries[i].incC_keyworker[j]);
			for (j = 0; j < 2; j++)
				fprintf(dat, ",%lg", TimeSeries[i].cumT_keyworker[j]);
			for (j = 0; j < 2; j++)
				fprintf(dat, ",%lg", TimeSeries[i].incD_keyworker[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	if (P.DoInftypeOutput)
	{
		sprintf(outname, "%s.inftype.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t,R");
		for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, ",Rtype_%i", j);
		for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, ",incItype_%i", j);
		for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, ",Rage_%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg,%lg", TimeSeries[i].t, TimeSeries[i].Rdenom);
			for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, ",%lg", TimeSeries[i].Rtype[j]);
			for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, ",%lg", TimeSeries[i].incItype[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, ",%lg", TimeSeries[i].Rage[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	if (P.DoROutput)
	{
		sprintf(outname, "%s.R0.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		for (i = 0; i < MAX_SEC_REC; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 0; j < MAX_GEN_REC; j++)
				fprintf(dat, ",%lg", indivR0[i][j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
	{
		t = 0;
		for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
			t += inf_household[i][j];
		inf_household[i][0] = denom_household[i] - t;
	}
	for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
	{
		t = 0;
		for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
			t += case_household[i][j];
		case_household[i][0] = denom_household[i] - t;
	}
	if (P.DoHouseholdOutput)
	{
		sprintf(outname, "%s.household.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
			fprintf(dat, "\t%i", i);
		fprintf(dat, "\n");
		for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				fprintf(dat, "\t%lg", inf_household[j][i]);
			fprintf(dat, "\n");
		}
		fprintf(dat, "\n");
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
			fprintf(dat, "\t%i", i);
		fprintf(dat, "\n");
		for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				fprintf(dat, "\t%lg", case_household[j][i]);
			fprintf(dat, "\n");
		}
		fprintf(dat, "\n");
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
			fprintf(dat, "\t%lg", denom_household[i]);
		fprintf(dat, "\n");
		fclose(dat);
	}

	if ((P.DoVaccOutput) & (P.OutbreakDetected))
	{
		sprintf(outname, "%s.vacc.dose.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "num vaccine doses");
		for (i = 0; i < NUM_VACCDOSE_GROUPS; i++)
			fprintf(dat, ",%i-%i", VACCDOSE_WIDTH * i, VACCDOSE_WIDTH * (i + 1));
		fprintf(dat, "\n");
		fprintf(dat, "frequency");
		for (i = 0; i < NUM_VACCDOSE_GROUPS; i++)
		{
			fprintf(dat, ",%lg", vaccdose_dist[i]);
		}
		fprintf(dat, "\n");
		fclose(dat);

		sprintf(outname, "%s.vacc.dosecell.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "num vaccine doses");
		for (i = 0; i < NUM_VACCDOSE_GROUPS; i++)
			fprintf(dat, ",%i-%i", VACCDOSECELL_WIDTH * i, VACCDOSE_WIDTH * (i + 1));
		fprintf(dat, "\n");
		fprintf(dat, "frequency");
		for (i = 0; i < NUM_VACCDOSE_GROUPS; i++)
		{
			fprintf(dat, ",%lg", vaccdosecell_dist[i]);
		}
		fprintf(dat, "\n");
		fclose(dat);

		sprintf(outname, "%s.vacc.dosering.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "num vaccine doses");
		for (i = 0; i < NUM_VACCDOSE_GROUPS; i++)
			fprintf(dat, ",%i-%i", VACCDOSE_WIDTH * i, VACCDOSE_WIDTH * (i + 1));
		fprintf(dat, "\n");
		fprintf(dat, "frequency");
		for (i = 0; i < NUM_VACCDOSE_GROUPS; i++)
		{
			fprintf(dat, ",%lg", vaccdosering_dist[i]);
		}
		fprintf(dat, "\n");
		fclose(dat);

		sprintf(outname, "%s.vacc.dist.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "max distance from trigger case");
		for (i = 0; i < NUM_VACCDIST_GROUPS; i++)
			fprintf(dat, ",%i-%i", VACCDIST_WIDTH * i, VACCDIST_WIDTH * (i + 1));
		fprintf(dat, "\n");
		fprintf(dat, "frequency");
		for (i = 0; i < NUM_VACCDIST_GROUPS; i++)
		{
			fprintf(dat, ",%lg", vaccdistance_dist[i]);
		}
		fprintf(dat, "\n");
		fclose(dat);

		sprintf(outname, "%s.vacc.pop.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "log population in trigger cell");
		for (i = 0; i < NUM_VACCDIST_GROUPS; i++)
			fprintf(dat, ",%i-", POPDIST_WIDTH * i);
		//fprintf(dat, ",%i - %i", POPDIST_WIDTH * i, POPDIST_WIDTH * (i + 1));
		fprintf(dat, "\n");
		fprintf(dat, "frequency");
		for (i = 0; i < NUM_POP_GROUPS; i++)
		{
			fprintf(dat, ",%lg", vaccpop_dist[i]);
		}
		fprintf(dat, "\n");
		fclose(dat);

		sprintf(outname, "%s.vacc.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "Mcell ID,Population,Num triggers,Total doses,Population vacced,Min dist,Doses at min dist,Time at min dist,Max dist,Doses at max dist,Time at max dist");
		for (i = 0; i < P.NMCP; i++)
		{
			j = (int)(McellLookup[i] - Mcells);
			if (Mcells[j].ntriggervacc > 0)
			{
				fprintf(dat, "\n");
				fprintf(dat, "%i,%i,%i,%i,%i,%lg,%i,%lg,%lg,%i,%lg", j, Mcells[j].n, Mcells[j].ntriggervacc, Mcells[j].totalvacc, Mcells[j].popvacc, Mcells[j].minvaccdist, Mcells[j].minvaccdist_dose, Mcells[j].minvaccdist_t, Mcells[j].maxvaccdist, Mcells[j].maxvaccdist_dose, Mcells[j].maxvaccdist_t);
			}
		}
		fprintf(dat, "\n");
		fclose(dat);
	}

	if (P.DoCountryOutput)
	{
		sprintf(outname, "%s.country.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		for (i = 0; i < MAX_COUNTRIES; i++)
			fprintf(dat, "%i\t%lg\n", i, infcountry[i]);
		fclose(dat);
	}
	//added code to output events, if required: ggilani - 10/10/2014
	if (P.DoRecordInfEvents)
	{
		sprintf(outname, "%s.infevents.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		//fprintf(dat,"t\tind_infectee\tx_infectee\ty_infectee\tt_infector\tind_infector\tx_infector\tiy_infector\n");
		for (i = 0; i < *nEvents; i++)
		{
			fprintf(dat, "%lg,%i,%lg,%lg,%lg,%i,%lg,%lg\n",
				InfEventLog[i].t, InfEventLog[i].infectee_ind, InfEventLog[i].infectee_x, InfEventLog[i].infectee_y, InfEventLog[i].t_infector, InfEventLog[i].infector_ind, InfEventLog[i].infector_x, InfEventLog[i].infector_y);
		}
		fclose(dat);
	}
	if (P.DoInfectionTree)
	{
		sprintf(outname, "%s.tree.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "id,onset,infector,infector_onset,district,detected,outcome\n");
		for (i = 0; i < P.N; i++)
			if (Hosts[i].infect_type % INFECT_TYPE_MASK > 0)
				if (Hosts[i].infector > 0)
					fprintf(dat, "%i,%i,%i,%i,%i,%i,%i\n", i, Hosts[i].latent_time, Hosts[i].infector, Hosts[Hosts[i].infector].latent_time, AdUnits[Mcells[Hosts[i].mcell].adunit].id, Hosts[i].detected, Hosts[i].to_die);
				else
					fprintf(dat, "%i,%i,%i,%i,%i,%i,%i\n", i, Hosts[i].latent_time, Hosts[i].infector, -1, AdUnits[Mcells[Hosts[i].mcell].adunit].id, Hosts[i].detected, Hosts[i].to_die);
		//fprintf(dat,"%i,%i,%i,%i\n",i,Hosts[i].infector,Hosts[i].infect_type%INFECT_TYPE_MASK,(int) HOST_AGE_YEAR(i));
		fclose(dat);
	}
#if defined(WIN32_BM) || defined(IMAGE_MAGICK)
	static int dm[13] = { 0,31,28,31,30,31,30,31,31,30,31,30,31 };
	int d, m, y, dml, f;
#ifdef WIN32_BM
	if (P.OutputBitmap == 1) CloseAvi(avi);
	if ((TimeSeries[P.NumSamples - 1].extinct) && (P.OutputNonExtinct))
	{
		sprintf(outname, "%s.avi", OutFile);
		DeleteFile(outname);
	}
#endif
	if (P.OutputBitmap >= 1)
	{
		// Generate Google Earth .kml file
#ifdef WIN32_BM
		sprintf(outname, "%s.ge.kml", OutFile); //sprintf(outname,"%s.ge\\%s.kml",OutFileBase,OutFile);
#else	
		sprintf(outname, "%s.ge.kml", OutFile);
#endif
		if (!(dat = fopen(outname, "w")))
		{
			ERR_CRITICAL("Unable to open output kml file\n");
		}
		fprintf(dat, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://earth.google.com/kml/2.2\">\n<Document>\n");
		fprintf(dat, "<name>%s</name>\n", OutFile);
		y = 2009;
		m = 1;
		d = 1;
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "<GroundOverlay>\n<name>Snapshot %i</name>\n", i + 1);
			fprintf(dat, "<TimeSpan>\n<begin>%i-%02i-%02iT00:00:00Z</begin>\n", y, m, d);
			d += (int)P.SampleStep; // SampleStep has to be an integer here.
			do
			{
				f = 1;
				dml = dm[m];
				if ((m == 2) && (y % 4 == 0)) dml = 29;
				if (d > dml)
				{
					m++;
					if (m > 12)
					{
						m -= 12;
						y++;
					}
					d -= dml;
					f = 0;
				}
			} while (!f);
			fprintf(dat, "<end>%i-%02i-%02iT00:00:00Z</end>\n</TimeSpan>\n", y, m, d);
			sprintf(outname, "%s.%i.png", OutFile, i + 1);
			fprintf(dat, "<Icon>\n<href>%s</href>\n</Icon>\n", outname);
			fprintf(dat, "<LatLonBox>\n<north>%lg</north>\n<south>%lg</south>\n<east>%lg</east>\n<west>%lg</west>\n</LatLonBox>\n",
				P.SpatialBoundingBox[3], P.SpatialBoundingBox[1], P.SpatialBoundingBox[2], P.SpatialBoundingBox[0]);
			fprintf(dat, "</GroundOverlay>\n");
		}
		fprintf(dat, "</Document>\n</kml>\n");
		fclose(dat);
	}
#endif
}


void SaveSummaryResults(void)
{
	int i, j;
	double c, t;
	FILE* dat;
	char outname[1024];

	c = 1 / ((double)(P.NRactE + P.NRactNE));
	sprintf(outname, "%s.csv", OutFile);
	if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat, "t,S,L,I,R,D,incI,incR,incFC,incFI,incC,incDC,incTC,incETU,cumT,cumTmax,cumTP,cumV,cumVmax,Extinct,Detected,rmsRad,maxRad,vS,vI,vR,vD,vincI,vincR,vincFC,vincFI,tvincC,vincDC,vincTC,vincH,vrmsRad,vmaxRad,%i,%i,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
		P.NRactNE, P.NRactE, P.R0household, P.R0places, P.R0spatial, c * PeakHeightSum, c * PeakHeightSS - c * c * PeakHeightSum * PeakHeightSum, c * PeakTimeSum, c * PeakTimeSS - c * c * PeakTimeSum * PeakTimeSum);
	c = 1 / ((double)P.NRactual);
	//added this for sake of test border control
	//c=1;
	for (i = 0; i < P.NumSamples; i++)
	{
		fprintf(dat, "%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,",
			c * TSMean[i].t, c * TSMean[i].S, c * TSMean[i].L, c * TSMean[i].I, c * TSMean[i].R,
			c * TSMean[i].D, c * TSMean[i].incI, c * TSMean[i].incR, c * TSMean[i].incFC, c * TSMean[i].incFI, c * TSMean[i].incC, c * TSMean[i].incDC, c * TSMean[i].incTC, c * TSMean[i].incETU, //added incidence of funeral transmission and hospitalisation
			c * TSMean[i].cumT, TSMean[i].cumTmax, c * TSMean[i].cumTP, c * TSMean[i].cumV, TSMean[i].cumVmax, c * TSMean[i].extinct, c * TSMean[i].detected, c * TSMean[i].rmsRad, c * TSMean[i].maxRad);
		fprintf(dat, "%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
			c * TSVar[i].S - c * c * TSMean[i].S * TSMean[i].S,
			c * TSVar[i].I - c * c * TSMean[i].I * TSMean[i].I,
			c * TSVar[i].R - c * c * TSMean[i].R * TSMean[i].R,
			c * TSVar[i].D - c * c * TSMean[i].D * TSMean[i].D,
			c * TSVar[i].incI - c * c * TSMean[i].incI * TSMean[i].incI,
			c * TSVar[i].incR - c * c * TSMean[i].incR * TSMean[i].incR,
			c * TSVar[i].incD - c * c * TSMean[i].incFC * TSMean[i].incFC,
			c * TSVar[i].incFI - c * c * TSMean[i].incFI * TSMean[i].incFI, //added funeral transmissions
			c * TSVar[i].incC - c * c * TSMean[i].incC * TSMean[i].incC,
			c * TSVar[i].incDC - c * c * TSMean[i].incDC * TSMean[i].incDC, //added detected cases
			c * TSVar[i].incTC - c * c * TSMean[i].incTC * TSMean[i].incTC,
			c * TSVar[i].incETU - c * c * TSMean[i].incETU * TSMean[i].incETU, //added hospitalisation
			c * TSVar[i].rmsRad - c * c * TSMean[i].rmsRad * TSMean[i].rmsRad,
			c * TSVar[i].maxRad - c * c * TSMean[i].maxRad * TSMean[i].maxRad);
	}
	fclose(dat);
	if (P.DoControlOutput)
	{
		sprintf(outname, "%s.controls.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t\tS\tincC\tincTC\tincFC\tincFI\tincH\tcumT\tcumUT\tcumTP\tcumV\tincHQ\tincAC\tincAH\tincAA\tincACS\tincAPC\tincAPA\tincAPCS");
		for (j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\tprClosed_%i", j);
		fprintf(dat, "t\tvS\tvincC\tvincTC\tvincFC\tvincFI\tvincETU\tvcumT\tvcumUT\tvcumTP\tvcumV");
		for (j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\tvprClosed_%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg\t%lf\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg",
				c * TSMean[i].t, c * TSMean[i].S, c * TSMean[i].incC, c * TSMean[i].incTC, c * TSMean[i].incFC, c * TSMean[i].incFI, c * TSMean[i].incETU,
				c * TSMean[i].cumT, c * TSMean[i].cumUT, c * TSMean[i].cumTP, c * TSMean[i].cumV, c * TSMean[i].incHQ,
				c * TSMean[i].incAC, c * TSMean[i].incAH, c * TSMean[i].incAA, c * TSMean[i].incACS,
				c * TSMean[i].incAPC, c * TSMean[i].incAPA, c * TSMean[i].incAPCS);
			for (j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\t%lg", c * TSMean[i].PropPlacesClosed[j]);
			fprintf(dat, "\t%lf\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg",
				c * TSVar[i].S - c * c * TSMean[i].S * TSMean[i].S,
				c * TSVar[i].incC - c * c * TSMean[i].incC * TSMean[i].incC,
				c * TSVar[i].incTC - c * c * TSMean[i].incTC * TSMean[i].incTC,
				c * TSVar[i].incFC - c * c * TSMean[i].incFC * TSMean[i].incFC,
				c * TSVar[i].incFI - c * c * TSMean[i].incFI * TSMean[i].incFI,
				c * TSVar[i].incETU - c * c * TSMean[i].incETU * TSMean[i].incETU,
				c * TSVar[i].cumT - c * c * TSMean[i].cumT * TSMean[i].cumT,
				c * TSVar[i].cumUT - c * c * TSMean[i].cumUT * TSMean[i].cumUT,
				c * TSVar[i].cumTP - c * c * TSMean[i].cumTP * TSMean[i].cumTP,
				c * TSVar[i].cumV - c * c * TSMean[i].cumV * TSMean[i].cumV);
			for (j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\t%lg", TSVar[i].PropPlacesClosed[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	if (P.DoAgeOutput)
	{
		sprintf(outname, "%s.age.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",I%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",C%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",D%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",DC%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",ETU%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, ",V%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg", c * TSMean[i].t);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", c * TSMean[i].incIa[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", c * TSMean[i].incCa[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", c * TSMean[i].incDa[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", c * TSMean[i].incDCa[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", c * TSMean[i].incETUa[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, ",%lg", c * TSMean[i].incVa[j]);
			fprintf(dat, "\n");
		}
		fprintf(dat, "dist");
		for (j = 0; j < NUM_AGE_GROUPS; j++)
			fprintf(dat, ",%lg", AgeDist[j]);
		fprintf(dat, "\n");
		fclose(dat);
	}
	if ((P.DoAdUnits) && (P.DoAdunitOutput))
	{
		sprintf(outname, "%s.adunit.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",I_%s", AdUnits[i].ad_name);
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",C_%s", AdUnits[i].ad_name);
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",DC_%s", AdUnits[i].ad_name); //added detected cases: ggilani 03/02/15
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",T_%s", AdUnits[i].ad_name);
		if ((P.DoHospitalisation) & (P.DoETUByAdUnit))
		{
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",ETU_%s", AdUnits[i].ad_name); //added hospitalisation
		}
		if (P.DoContactTracing)
		{
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",CT_%s", AdUnits[i].ad_name); //added contact tracing
		}
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",%lg", P.PopByAdunit[i][0]);
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",%lg", P.PopByAdunit[i][1]);
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg", c * TSMean[i].t);
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, ",%lg", c * TSMean[i].incI_adunit[j]);
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, ",%lg", c * TSMean[i].incC_adunit[j]);
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, ",%lg", c * TSMean[i].incDC_adunit[j]); //added detected cases: ggilani 03/02/15
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, ",%lg", c * TSMean[i].cumT_adunit[j]);
			if ((P.DoHospitalisation) & (P.DoETUByAdUnit))
			{
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, ",%lg", c * TSMean[i].incETU_adunit[j]); //added hospitalisation
			}
			if (P.DoContactTracing)
			{
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, ",%lg", c * TSMean[i].incCT_adunit[j]); //added contact tracing
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, ",%lg", c * TSMean[i].incCC_adunit[j]); //added cases who are contacts
			}
			fprintf(dat, "\n");
		}
		fclose(dat);
		sprintf(outname, "%s.adunitVar.csv", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",I_%s", AdUnits[i].ad_name);
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",C_%s", AdUnits[i].ad_name);
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",DC_%s", AdUnits[i].ad_name); //added detected cases: ggilani 03/02/15
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",T_%s", AdUnits[i].ad_name);
		if ((P.DoHospitalisation) & (P.DoETUByAdUnit))
		{
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",ETU_%s", AdUnits[i].ad_name); //added hospitalisation
		}
		if (P.DoContactTracing)
		{
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, ",CT_%s", AdUnits[i].ad_name); //added contact tracing
		}
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg", c * TSMean[i].t);
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, ",%lg", c * TSVar[i].incI_adunit[j] - c * c * TSMean[i].incI_adunit[j] * TSMean[i].incI_adunit[j]);
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, ",%lg", c * TSVar[i].incC_adunit[j] - c * c * TSMean[i].incC_adunit[j] * TSMean[i].incC_adunit[j]);
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, ",%lg", c * TSVar[i].incDC_adunit[j] - c * c * TSMean[i].incDC_adunit[j] * TSMean[i].incDC_adunit[j]); //added detected cases: ggilani 03/02/15
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, ",%lg", c * TSVar[i].cumT_adunit[j] - c * c * TSMean[i].cumT_adunit[j] * TSMean[i].cumT_adunit[j]);
			if ((P.DoHospitalisation) & (P.DoETUByAdUnit))
			{
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, ",%lg", c * TSVar[i].incETU_adunit[j] - c * c * TSMean[i].incETU_adunit[j] * TSMean[i].incETU_adunit[j]); //added hospitalisation
			}
			if (P.DoContactTracing)
			{
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, ",%lg", c * TSVar[i].incCT_adunit[j] - c * c * TSMean[i].incCT_adunit[j] * TSMean[i].incCT_adunit[j]); //added contact tracing
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, ",%lg", c * TSVar[i].incCC_adunit[j] - c * c * TSMean[i].incCC_adunit[j] * TSMean[i].incCC_adunit[j]); //added cases who are contacts
			}
			fprintf(dat, "\n");
		}
		fclose(dat);
	}

	if (P.EvolResistNumTypes > 1)
	{
		sprintf(outname, "%s.resist.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for (i = 0; i < P.EvolResistNumTypes; i++) fprintf(dat, "\tI%i", i);
		for (i = 0; i < P.EvolResistNumTypes; i++) fprintf(dat, "\tC%i", i);
		for (i = 0; i < P.EvolResistNumTypes; i++) fprintf(dat, "\tT%i", i);
		for (i = 0; i < P.EvolResistNumTypes; i++) fprintf(dat, "\tvI%i", i);
		for (i = 0; i < P.EvolResistNumTypes; i++) fprintf(dat, "\tvC%i", i);
		for (i = 0; i < P.EvolResistNumTypes; i++) fprintf(dat, "\tvT%i", i);
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg", c * TSMean[i].t);
			for (j = 0; j < P.EvolResistNumTypes; j++)
				fprintf(dat, "\t%lg", c * TSMean[i].incI_resist[j]);
			for (j = 0; j < P.EvolResistNumTypes; j++)
				fprintf(dat, "\t%lg", c * TSMean[i].incC_resist[j]);
			for (j = 0; j < P.EvolResistNumTypes; j++)
				fprintf(dat, "\t%lg", c * TSMean[i].cumT_resist[j]);
			for (j = 0; j < P.EvolResistNumTypes; j++)
				fprintf(dat, "\t%lg", c * TSVar[i].incI_resist[j] - c * c * TSMean[i].incI_resist[j] * TSMean[i].incI_resist[j]);
			for (j = 0; j < P.EvolResistNumTypes; j++)
				fprintf(dat, "\t%lg", c * TSVar[i].incC_resist[j] - c * c * TSMean[i].incC_resist[j] * TSMean[i].incC_resist[j]);
			for (j = 0; j < P.EvolResistNumTypes; j++)
				fprintf(dat, "\t%lg", c * TSVar[i].cumT_resist[j] - c * c * TSMean[i].cumT_resist[j] * TSMean[i].cumT_resist[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	if (P.KeyWorkerProphTimeStartBase < P.SampleTime)
	{
		sprintf(outname, "%s.keyworker.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for (i = 0; i < 2; i++) fprintf(dat, "\tI%i", i);
		for (i = 0; i < 2; i++) fprintf(dat, "\tC%i", i);
		for (i = 0; i < 2; i++) fprintf(dat, "\tT%i", i);
		for (i = 0; i < 2; i++) fprintf(dat, "\tvI%i", i);
		for (i = 0; i < 2; i++) fprintf(dat, "\tvC%i", i);
		for (i = 0; i < 2; i++) fprintf(dat, "\tvT%i", i);
		fprintf(dat, "\t%i\t%i\n", P.KeyWorkerNum, P.KeyWorkerIncHouseNum);
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg", c * TSMean[i].t);
			for (j = 0; j < 2; j++)
				fprintf(dat, "\t%lg", c * TSMean[i].incI_keyworker[j]);
			for (j = 0; j < 2; j++)
				fprintf(dat, "\t%lg", c * TSMean[i].incC_keyworker[j]);
			for (j = 0; j < 2; j++)
				fprintf(dat, "\t%lg", c * TSMean[i].cumT_keyworker[j]);
			for (j = 0; j < 2; j++)
				fprintf(dat, "\t%lg", c * TSVar[i].incI_keyworker[j] - c * c * TSMean[i].incI_keyworker[j] * TSMean[i].incI_keyworker[j]);
			for (j = 0; j < 2; j++)
				fprintf(dat, "\t%lg", c * TSVar[i].incC_keyworker[j] - c * c * TSMean[i].incC_keyworker[j] * TSMean[i].incC_keyworker[j]);
			for (j = 0; j < 2; j++)
				fprintf(dat, "\t%lg", c * TSVar[i].cumT_keyworker[j] - c * c * TSMean[i].cumT_keyworker[j] * TSMean[i].cumT_keyworker[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	if (P.DoInftypeOutput)
	{
		sprintf(outname, "%s.inftype.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t,R");
		for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\tRtype_%i", j);
		for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\tincItype_%i", j);
		for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\tRage_%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg,%lg", c * TSMean[i].t, c * TSMean[i].Rdenom);
			for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\t%lg", c * TSMean[i].Rtype[j]);
			for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\t%lg", c * TSMean[i].incItype[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%lg", c * TSMean[i].Rage[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	if (P.DoROutput)
	{
		sprintf(outname, "%s.R0.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		for (i = 0; i < MAX_SEC_REC; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 0; j < MAX_GEN_REC; j++)
				fprintf(dat, "\t%lg", c * indivR0_av[i][j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	if (P.DoHouseholdOutput)
	{
		sprintf(outname, "%s.household.xls", OutFile);
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			t = 0;
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				t += inf_household_av[i][j];
			inf_household_av[i][0] = denom_household[i] / c - t;
		}
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			t = 0;
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				t += case_household_av[i][j];
			case_household_av[i][0] = denom_household[i] / c - t;
		}
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
			fprintf(dat, "\t%i", i);
		fprintf(dat, "\n");
		for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				fprintf(dat, "\t%lg", inf_household_av[j][i] * c);
			fprintf(dat, "\n");
		}
		fprintf(dat, "\n");
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
			fprintf(dat, "\t%i", i);
		fprintf(dat, "\n");
		for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				fprintf(dat, "\t%lg", case_household_av[j][i] * c);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
	if (P.DoCountryOutput)
	{
		sprintf(outname, "%s.country.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		for (i = 0; i < MAX_COUNTRIES; i++)
			fprintf(dat, "%i\t%lg\t%lg\n", i, infcountry_av[i] * c, infcountry_num[i] * c);
		fclose(dat);
	}

}

/* function: SaveRandomSeeds(void)
 *
 * Purpose: outputs the random seeds used for each run to a file
 * Parameter: none
 * Returns: none
 *
 * Author: ggilani, 09/03/17
 */
void SaveRandomSeeds(void)
{
	FILE* dat;
	char outname[1024];

	sprintf(outname, "%s.seeds.csv", OutFile);
	if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat, "%i,%i\n", P.newseed1, P.newseed2);
	fclose(dat);
}


/* function: SaveEvents(void)
 *
 * Purpose: outputs event log to a csv file if required
 * Parameters: none
 * Returns: none
 *
 * Author: ggilani, 15/10/2014
 */
void SaveEvents(void)
{
	int i;
	FILE* dat;
	char outname[1024];

	sprintf(outname, "%s.infevents.csv", OutFile);
	if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat, "type,t,thread,ind_infectee,cell_infectee,listpos_infectee,adunit_infectee,x_infectee,y_infectee,t_infector,ind_infector,cell_infector\n");
	for (i = 0; i < *nEvents; i++)
	{
		fprintf(dat, "%i,%lg,%i,%i,%i,%i,%i,%lg,%lg,%lg,%i,%i\n",
			InfEventLog[i].type, InfEventLog[i].t, InfEventLog[i].thread, InfEventLog[i].infectee_ind, InfEventLog[i].infectee_cell, InfEventLog[i].listpos, InfEventLog[i].infectee_adunit, InfEventLog[i].infectee_x, InfEventLog[i].infectee_y, InfEventLog[i].t_infector, InfEventLog[i].infector_ind, InfEventLog[i].infector_cell);
	}
	fclose(dat);
}
