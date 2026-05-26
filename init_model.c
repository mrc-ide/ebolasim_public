/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "model_constants.h"
#include "model_structs.h"
#include "model_fns.h"
#include "model_globals_ext.h"
#include "binio.h"


void InitModel(int run) //passing run number so we can save run number in the infection event log: ggilani - 15/10/2014
{
	int i, j, k, l, m, i2, j2, tn, nim, stt, stp, b, nhcws, nflws;
	double t, s;
	char buf[200];

	if (P.OutputBitmap)
	{
#ifdef WIN32_BM
		if (P.OutputBitmap == 1)
		{
			sprintf(buf, "%s.avi", OutFile);
			avi = CreateAvi(buf, P.BitmapMovieFrame, NULL);
		}
#endif
		for (i = 0; i < bmh->imagesize; i++)
		{
			bmi2[i] = bmi3[i] = bmi4[i] = 0;
		}
	}

	ns = 0;
	State.S = P.N;
	State.L = State.I = State.R = 0;
	State.cumI = State.cumR = State.cumC = State.cumFC = State.cumFI = State.cumETU = State.cumH = State.cumCT = State.cumCC = State.cumTC = State.cumD = State.cumDC = State.trigDC = State.cumSDB = State.cumDD
		= State.cumInf_h = State.cumInf_n = State.cumInf_s = State.cumHQ
		= State.cumAC = State.cumAH = State.cumAA = State.cumACS
		= State.cumAPC = State.cumAPA = State.cumAPCS = 0;
	State.cumT = State.cumUT = State.cumTP = State.cumV = State.sumRad2 = State.maxRad2 = State.cumV_daily = State.cumVG = 0; //added State.cumVG
	State.mvacc_cum = 0;
	State.NumBeds = 0;
	State.nringvacc_queue = State.ngeovacc_queue = State.nvacc_queue = State.ringvacc_cum = State.geovacc_cum = State.ringvacc_ind = State.geovacc_ind = State.vacc_ind = State.vacc_cum = 0; //reset ring and geo vaccination variables, 21/08/19

	for (i = 0; i < NUM_AGE_GROUPS; i++) State.cumCa[i] = State.cumIa[i] = State.cumDa[i] = State.cumDCa[i] = State.cumETUa[i] = State.cumHa[i] = State.cumVa[i] = 0; //adding det case, hosp, vacc by age: ggilani 22/02/22
	for (i = 0; i < P.EvolResistNumTypes; i++) State.cumC_resist[i] = State.cumI_resist[i] = State.cumT_resist[i] = 0;
	for (i = 0; i < 2; i++) State.cumC_keyworker[i] = State.cumI_keyworker[i] = State.cumT_keyworker[i] = 0;
	for (i = 0; i < NUM_PLACE_TYPES; i++) State.NumPlacesClosed[i] = 0;
	for (i = 0; i < INFECT_TYPE_MASK; i++) State.cumItype[i] = 0;
	//initialise cumulative case counts per country to zero: ggilani 12/11/14
	for (i = 0; i < MAX_COUNTRIES; i++) State.cumC_country[i] = 0;
	if (P.DoAdUnits)
		for (i = 0; i <= P.NumAdunits; i++)
		{
			State.cumI_adunit[i] = State.cumC_adunit[i] = State.cumT_adunit[i] = State.cumETU_adunit[i] = State.cumH_adunit[i] = State.cumDC_adunit[i] = State.cumDD_adunit[i] = State.cumSDB_adunit[i] = State.cumDR_adunit[i] = State.cumCT_adunit[i] = State.cumV_adunit[i] = State.cumVG_adunit[i] = State.cumC_adunit[i] = State.cumCC_adunit[i] = 0; //added hospitalisation, added detected cases, contact tracing per adunit, cases who are contacts: ggilani 03/02/15, 15/06/17
			State.ETU_adunit[i] = State.H_adunit[i] = State.NumBeds_adunits[i] = 0;
			AdUnits[i].place_close_trig = 0;
			AdUnits[i].revacc = 0;
			AdUnits[i].currentETUBeds = 0; //reset occupied beds to zero;
			AdUnits[i].ETUbedsActive = 0;
			AdUnits[i].totalETUBeds = 0;
			AdUnits[i].lastCaseDay = 0;
			//added this for reactive hospital beds
			if ((P.DoReactETUBeds))//&&(P.StartTimeReactiveBeds<1))
			{
				//AdUnits[i].totalBeds = 0; // AdUnits[i].initialBeds;
				//AdUnits[i].timeBedsAvailable = 1e10;
				AdUnits[i].flagETUBeds = 0;
				//AdUnits[i].bedsActive = 0;
				AdUnits[i].nextETUBeds = 0;
				AdUnits[i].nextTimeToETUBeds = 1e10;
			}
			//if (P.DoUpdateCaseDetection)
			//{
			//	AdUnits[i].caseDetectRate = AdUnits[i].caseDetectInit;
			//	P.CurrIndUpdateCaseDetect = 0;
			//	P.UpdateCaseDetectionByCasesFlag = 0;
			//	P.CaseDetectionRate = P.RR1;
			//}
			AdUnits[i].contactTraceThresholdCrossed = 0;
			AdUnits[i].nct = 0; //no-one being contact traced at beginning of run
			AdUnits[i].nct_queue = 0; //no-one in contact tracing queue at beginning of run
			AdUnits[i].nh_queue = 0; //no-one in hospital queue at beginning of run
			AdUnits[i].contactTraceStartDay = 1e6;
			AdUnits[i].contactTraceCapacity = P.AdunitCTCapacity;
			AdUnits[i].nextTimeToSDB = 0;
			AdUnits[i].maxSDB = P.AdunitSDBCapacity;
		}
	for (j = 0; j < MAX_NUM_THREADS; j++)
	{
		StateT[j].L = StateT[j].I = StateT[j].R = 0;
		StateT[j].cumI = StateT[j].cumR = StateT[j].cumC = StateT[j].cumFC = StateT[j].cumFI = StateT[j].cumETU = StateT[j].cumH = StateT[j].cumCT = StateT[j].cumCC = StateT[j].cumTC = StateT[j].cumD = StateT[j].cumDC = StateT[j].cumDD = StateT[j].cumSDB = 0; //added setting funeral infections cumFI to zero: ggilani 24/10/14
		StateT[j].cumInf_h = StateT[j].cumInf_n = StateT[j].cumInf_s = StateT[j].cumHQ = StateT[j].cumAC = StateT[j].cumACS = StateT[j].cumAH = StateT[j].cumAA = StateT[j].cumAPC = StateT[j].cumAPA = StateT[j].cumAPCS = 0;
		StateT[j].cumT = StateT[j].cumUT = StateT[j].cumTP = StateT[j].cumV = StateT[j].sumRad2 = StateT[j].maxRad2 = StateT[j].cumV_daily = 0;
		StateT[j].nringvacc_queue = StateT[j].nvacc_queue = StateT[j].vacc_cum = StateT[j].ngeovacc_queue = StateT[j].ringvacc_cum = StateT[j].geovacc_cum = 0; //reset ring and geo vaccination variables, 21/08/19
		for (i = 0; i < NUM_AGE_GROUPS; i++) StateT[j].cumCa[i] = StateT[j].cumIa[i] = StateT[j].cumDa[i] = StateT[j].cumDCa[i] = StateT[j].cumETUa[i] = StateT[j].cumHa[i] = StateT[j].cumVa[i] = 0; //adding det cases, hosp, vacc by age: ggilani 22/02/22
		for (i = 0; i < P.EvolResistNumTypes; i++) StateT[j].cumC_resist[i] = StateT[j].cumI_resist[i] = StateT[j].cumT_resist[i] = 0;
		for (i = 0; i < 2; i++) StateT[j].cumC_keyworker[i] = StateT[j].cumI_keyworker[i] = StateT[j].cumT_keyworker[i] = 0;
		for (i = 0; i < NUM_PLACE_TYPES; i++) StateT[j].NumPlacesClosed[i] = 0;
		for (i = 0; i < INFECT_TYPE_MASK; i++) StateT[j].cumItype[i] = 0;
		//initialise cumulative case counts per country per thread to zero: ggilani 12/11/14
		for (i = 0; i < MAX_COUNTRIES; i++) StateT[j].cumC_country[i] = 0;
		if (P.DoAdUnits)
			for (i = 0; i <= P.NumAdunits; i++)
				StateT[j].cumI_adunit[i] = StateT[j].cumC_adunit[i] = StateT[j].cumT_adunit[i] = StateT[j].cumETU_adunit[i] = StateT[j].cumH_adunit[i] = StateT[j].ETU_adunit[i] = StateT[j].H_adunit[i] = StateT[j].cumDC_adunit[i] = StateT[j].cumD_adunit[i] = StateT[j].cumDD_adunit[i] = StateT[j].cumSDB_adunit[i] = StateT[j].cumDR_adunit[i] = StateT[j].cumCT_adunit[i] = StateT[j].cumV_adunit[i] = StateT[j].cumV_adunit[i] = StateT[j].cumC_adunit[i] = StateT[j].cumCC_adunit[i] = StateT[j].nct_queue[i] = 0; //added hospitalisation, detected cases, contact tracing per adunit, cases who are contacts: ggilani 03/02/15, 15/06/17
	}
	nim = 0;

	//added this to reorder the susceptible array to be in the same order as the member array - ggilani 12/03/17
	//put this back in to try and fix a bug - it seems to be working but I don't think I should need this - ggilani  22/08/19
	//for (i = 0; i < P.NC; i++)
	//{
	//	for (j = 0; j < Cells[i].n; j++) Cells[i].susceptible[j] = Cells[i].members[j];
	//}


#pragma omp parallel for private(i,b,j,k,l,m,tn,stt,stp) reduction(+:nim) schedule(static,1)
	for (tn = 0; tn < P.NumThreads; tn++)
	{
		stp = P.NC / P.NumThreads + 1;
		stt = (tn + 1) * stp;
		if (stt > P.NC) stt = P.NC;
		for (i = tn * stp; i < stt; i++)
		{
			if ((Cells[i].tot_treat != 0) || (Cells[i].tot_vacc != 0) || (Cells[i].S != Cells[i].n) || (Cells[i].D > 0) || (Cells[i].R > 0))
			{
				Cells[i].S = Cells[i].n;
				Cells[i].L = Cells[i].I = Cells[i].R = Cells[i].cumTC = Cells[i].D = 0;
				Cells[i].infected = Cells[i].latent = Cells[i].susceptible + Cells[i].S;
				Cells[i].tot_treat = Cells[i].tot_vacc = 0;

				for (l = 0; l < MAX_INTERVENTION_TYPES; l++) Cells[i].CurInterv[l] = -1;
				for (j = 0; j < Cells[i].n; j++)
				{
					k = Cells[i].members[j];
					Cells[i].susceptible[j] = k; //added this in here instead
					if (Households[Hosts[k].hh].stockpile != 0) Households[Hosts[k].hh].stockpile = 1;
					if (P.DoAirports) Hosts[k].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
					Hosts[k].quar_start_time =
						Hosts[k].isolation_start_time = Hosts[k].absent_start_time = USHRT_MAX - 1;
					Hosts[k].absent_stop_time = 0;
					Hosts[k].quar_comply = 2;
					if ((AdUnits[Mcells[Hosts[k].mcell].adunit].id / P.CountryDivisor) == P.TargetCountry)
					{
						Hosts[k].susc = P.RelativeSusceptibilityGuinea;
					}
					else if ((AdUnits[Mcells[Hosts[k].mcell].adunit].id / P.CountryDivisor) == P.TargetCountry2)
					{
						Hosts[k].susc = P.RelativeSusceptibilityLiberia;
					}
					else
					{
						Hosts[k].susc = 1.0;
					}
					Hosts[k].to_die = 0;
					Hosts[k].Travelling = 0;
					Hosts[k].detected = 0; //set detected to zero initially: ggilani - 19/02/15
					Hosts[k].dayDetected = USHRT_MAX - 1; //set day on which each case is detected to zero initially: ggilani - 23/06/15
					Hosts[k].inf = 0;
					Hosts[k].listpos = j;
					Hosts[k].treat_stop_time = Hosts[k].num_treats = Hosts[k].contactTraced_end_time = 0;
					Hosts[k].vacc_start_time = Hosts[k].treat_start_time = Hosts[k].contactTraced_start_time = USHRT_MAX - 1;
					Hosts[k].revacc = 0;
					Hosts[k].latent_time = Hosts[k].recovery_time = Hosts[k].hospital_time = Hosts[k].detect_time = 0; //also set hospitalisation time to zero: ggilani 28/10/2014
					Hosts[k].hospitalised = Hosts[k].etu = 0; //set hospitalised flag to zero: ggilani 28/10/14
					Hosts[k].contactTraced = 0;
					Hosts[k].vaccRing = 0; //set flag for vacc ring to zero: ggilani 29/05/19
					Hosts[k].ringCase = -1; //set flag for ring case to -1;
					Hosts[k].resist = 0;
					Hosts[k].infector = -1;
					Hosts[k].infect_type = 0;
					Hosts[k].infectiousMult = 1; //reset to 1 - this is changed when funeral transmission temporarily increases infectiousness

				}
				// Next loop needs to count down for DoImmune host list reordering to work
				for (j = Cells[i].n - 1; j >= 0; j--)
				{
					k = Cells[i].members[j];
					if (P.DoWholeHouseholdImmunity)
					{
						if (P.InitialImmunity[0] != 0)
						{
							if (Households[Hosts[k].hh].FirstPerson == k)
							{
								for (m = 0; m < Households[Hosts[k].hh].nh; m++)
									Hosts[k + m].inf = 0;
								if ((P.InitialImmunity[0] == 1) || (ranf_mt(tn) < P.InitialImmunity[0]))
									for (m = Households[Hosts[k].hh].nh - 1; m >= 0; m--)
										DoImmune(k + m);
							}
						}
					}
					Hosts[k].inf = 0;
					if (P.DoInitEquilib)
					{
						if (P.DoSIS)
						{
							if (P.SuscReductionFactorPerInfection > 0)
								Hosts[k].susc = exp(-P.MeanAnnualDeathRate * ((double)HOST_AGE_YEAR(k)));
							else
								Hosts[k].susc = (ranf_mt(tn) > exp(-P.MeanAnnualDeathRate * ((double)HOST_AGE_YEAR(k)))) ? 0 : 1;
						}
						else if (ranf_mt(tn) > exp(-P.MeanAnnualDeathRate * ((double)HOST_AGE_YEAR(k))))
							DoImmune(k);
					}
					else
					{
						m = HOST_AGE_GROUP(k);
						if (P.DoSIS)
						{
							if (P.SuscReductionFactorPerInfection > 0)
								Hosts[k].susc = 1 - P.InitialImmunity[m];
							else
								nim += (Hosts[k].susc = (ranf_mt(tn) < P.InitialImmunity[m]) ? 0 : 1);
						}
						else
						{
							if ((P.InitialImmunity[m] == 1) || ((P.InitialImmunity[m] > 0) && (ranf_mt(tn) < P.InitialImmunity[m]))) { DoImmune(k); nim += 1; }
						}
					}
				}
			}
		}
	}

	fprintf(stderr, "Finished cell init - %i people assigned as immune.\n", nim);


#pragma omp parallel for private(i,j,k,i2,j2,l) schedule(static,500)
	for (l = 0; l < P.NMCP; l++)
	{
		i = (int)(McellLookup[l] - Mcells);
		Mcells[i].vacc_start_time = Mcells[i].treat_start_time = USHRT_MAX - 1;
		//added these things, just to do some bookkeeping when testing out geo vaccination
		Mcells[i].ntriggervacc = 0;
		Mcells[i].popvacc = 0;
		Mcells[i].totalvacc = 0;
		Mcells[i].minvaccdist = 1e6;
		Mcells[i].maxvaccdist = -1;
		Mcells[i].minvaccdist_t = 1e6;
		Mcells[i].maxvaccdist_t = -1;
		Mcells[i].minvaccdist_dose = 1e6;
		Mcells[i].maxvaccdist_dose = -1;
		Mcells[i].treat_end_time = 0;
		Mcells[i].treat_trig = Mcells[i].vacc_trig = Mcells[i].vacc = Mcells[i].treat = 0;
		Mcells[i].place_trig = Mcells[i].move_trig = Mcells[i].socdist_trig = Mcells[i].keyworkerproph_trig =
			Mcells[i].placeclose = Mcells[i].moverest = Mcells[i].socdist = Mcells[i].keyworkerproph = 0;
		Mcells[i].move_start_time = USHRT_MAX - 1;
		Mcells[i].place_end_time = Mcells[i].move_end_time =
			Mcells[i].socdist_end_time = Mcells[i].keyworkerproph_end_time = 0;
		if (P.DoPlaces)
			for (j = 0; j < P.PlaceTypeNum; j++)
				for (k = 0; k < Mcells[i].np[j]; k++)
				{
					j2 = Mcells[i].places[j][k];
					Places[j][j2].treat = Places[j][j2].control_trig = 0;
					Places[j][j2].treat_end_time = Places[j][j2].close_end_time = 0;
					Places[j][j2].close_start_time = USHRT_MAX - 1;
#ifdef ABSENTEEISM_PLACE_CLOSURE
					Places[j][j2].AbsentLastUpdateTime = 0;
					for (i2 = 0; i2 < MAX_ABSENT_TIME; i2++) Places[j][j2].Absent[i2] = 0;
#endif
				}
	}


	for (i = 0; i < MAX_NUM_THREADS; i++)
	{
		for (j = 0; j < MAX_NUM_THREADS; j++)
			StateT[i].n_queue[j] = 0;
		for (j = 0; j < P.PlaceTypeNum; j++)
			StateT[i].np_queue[j] = 0;
	}
	if (DoInitUpdateProbs)
	{
		UpdateProbs(0);
		DoInitUpdateProbs = 0;
	}
	//initialise event log to zero at the beginning of every run: ggilani - 10/10/2014. UPDATE: 15/10/14 - we are now going to store all events from all realisations in one file
	if ((P.DoRecordInfEvents) && (P.RecordInfEventsPerRun))
	{
		*nEvents = 0;
		for (i = 0; i < P.MaxInfEvents; i++)
		{
			InfEventLog[i].t = InfEventLog[i].infectee_x = InfEventLog[i].infectee_y = InfEventLog[i].t_infector = 0.0;
			InfEventLog[i].infectee_ind = InfEventLog[i].infector_ind = 0;
			InfEventLog[i].infectee_adunit = InfEventLog[i].listpos = InfEventLog[i].infectee_cell = InfEventLog[i].infector_cell = InfEventLog[i].thread = 0;
		}
	}

	SeedInfection(0, P.NumInitialInfections, 0, run);
	P.ControlPropCasesId = P.PostAlertControlPropCasesId;
	//P.ControlPropCasesId=P.PreAlertControlPropCasesId;
	P.TreatTimeStart = 1e10;

	//  Mass Vacc now starts after outbreak alert trigger
	//	if(P.DoMassVacc)
	//		P.VaccTimeStart=P.VaccTimeStartBase;
	//	else
	P.DayExtinct = 0;
	P.VaccTimeStart = 1e10;
	P.VaccNewCoursesStartTime = 1e10;
	P.VaccNewCoursesEndTime = 1e10;
	P.VaccNewCoursesBoostStartTime = 1e10;
	P.MoveRestrTimeStart = 1e10;
	P.PlaceCloseTimeStart = 1e10;
	P.PlaceCloseTimeStart2 = 1e10;
	P.SocDistTimeStart = 1e10;
	P.AirportCloseTimeStart = 1e10;
	P.CaseIsolationTimeStart = 1e10;
	P.HQuarantineTimeStart = 1e10;
	P.KeyWorkerProphTimeStart = 1e10;
	//P.GeoVaccTimeStart = 1e10; //reset additional start times - ggilani 13/09/23
	//P.RingVaccTimeStart = 1e10;
	P.FuneralControlTimeStart = 1e10;
	P.ContactTracingTimeStart = 1e10;
	P.ETUTimeStart = 1e10;
	P.TreatMaxCourses = P.TreatMaxCoursesBase;
	P.VaccMaxCourses = P.VaccMaxCoursesBase;
	P.PlaceCloseDuration = P.PlaceCloseDurationBase;
	//if do hospitalisation, reset a couple of hospitalisation parameters
	P.CurrIndETUBeds = 0;
	P.CurrIndMeanTimeToHosp = 0;
	P.CurrIndMeanTimeToHospCT = 0;
	P.ResetSeedsFlag = 0; //added this to allow resetting seeds part way through run: ggilani 27/11/2019
	P.NVaccRingsActive = P.NVaccRings;
	P.CurrIndPropRingVacc = 0;
	P.OutbreakDetected = 0;
	P.PropHospSeek = P.PropHospSeekPreOutbreak;
	P.VaccDosePerDay = P.BaseVaccDosePerDay;
	P.VaccGeoDosePerDay = P.BaseVaccGeoDosePerDay;
	P.VaccDoseFlag = 1;
	P.UpdateVaccDosePerDay = 1;

	//vaccinate HCWs and FLWs
	if (P.IncludeHospitalPlaceType)
	{
		for (i = 0; i < P.Nplace[P.HospPlaceTypeNum]; i++)
		{
			for (j = 0; j < Places[P.HospPlaceTypeNum][i].nhcws; j++)
			{
				if (Hosts[Places[P.HospPlaceTypeNum][i].members[j]].vacc_accept < P.PropHCWFLWVacc)
				{
					DoVaccNoDelay(Places[P.HospPlaceTypeNum][i].members[j], -1 * ((int)P.TimeStepsPerDay * P.DayHCWFLWVacc));
				}
			}
			if (P.IncludeFLWs)
			{
				for (j = 0; j < Places[P.HospPlaceTypeNum][i].nflws; j++)
				{
					if (Hosts[Places[P.HospPlaceTypeNum][i].members[j]].vacc_accept < P.PropHCWFLWVacc)
					{
						DoVaccNoDelay(Places[P.HospPlaceTypeNum][i].flwmembers[j], -1 * ((int)P.TimeStepsPerDay * P.DayHCWFLWVacc));
					}
				}
			}
		}
	}

	for (i = 0; i < NUM_VACCDOSE_GROUPS; i++)
	{
		vaccdose_dist[i] = 0;
	}
	for (i = 0; i < NUM_VACCDOSECELL_GROUPS; i++)
	{
		vaccdose_dist[i] = 0;
	}
	for (i = 0; i < NUM_VACCDIST_GROUPS; i++)
	{
		vaccdistance_dist[i] = 0;
	}
	for (i = 0; i < NUM_VACCDIST_GROUPS; i++)
	{
		vaccdosering_dist[i] = 0;
	}
	for (i = 0; i < NUM_POP_GROUPS; i++)
	{
		vaccpop_dist[i] = 0;
	}

	fprintf(stderr, "Finished InitModel.\n");
}

void SeedInfection(double t, int* nsi, int rf, int run) //adding run number to pass it to event log
{
	int i, j, k, l, m, f, f2, n;

	f2 = ((t >= 0) && (P.DoImportsViaAirports));
	n = ((rf == 0) ? P.NumSeedLocations : 1);
	for (i = 0; i < n; i++)
	{
		if ((!P.DoRandomInitialInfectionLoc) || ((P.DoAllInitialInfectioninSameLoc) && (rf)))
		{
			//if we are importing cases to a specific location, and this is not initial seeding
			if ((rf == 1) && (P.DoImportToSpecLocation))
			{
				k = (int)(P.ImportLocation[0] / P.mcwidth);
				l = (int)(P.ImportLocation[1] / P.mcheight);
			}
			else
			{
				k = (int)(P.LocationInitialInfection[i][0] / P.mcwidth);
				l = (int)(P.LocationInitialInfection[i][1] / P.mcheight);
			}
			j = k * P.nmch + l;
			m = 0;
			for (k = 0; (k < nsi[i]) && (m < 10000); k++)
			{
				l = Mcells[j].members[(int)(ranf() * ((double)Mcells[j].n))];
				if (Hosts[l].inf == 0)
				{
					if (CalcPersonSusc(l, 0, 0, 0) > 0)
					{
						//only reset the initial location if rf==0, i.e. when initial seeds are being set, not when imported cases are being set
						if (rf == 0)
						{
							P.LocationInitialInfection[i][0] = Households[Hosts[l].hh].loc_x;
							P.LocationInitialInfection[i][1] = Households[Hosts[l].hh].loc_y;
						}
						Hosts[l].infector = -2; Hosts[l].infect_type = INFECT_TYPE_MASK - 1;
						Hosts[l].base_inf_level = 1.0;
						Hosts[l].resist = CalcSeedResist();
						DoInfect(l, t, 0, run);
						m = 0;
					}
				}
				else
				{
					k--; m++;
				}
			}
		}
		else if (P.DoAllInitialInfectioninSameLoc)
		{
			f = 0;
			do
			{
				m = 0;
				do
				{
					l = (int)(ranf() * ((double)P.N));
					j = Hosts[l].mcell;
					//fprintf(stderr,"%i ",AdUnits[Mcells[j].adunit].id);
				} while ((Mcells[j].n < nsi[i]) || (Mcells[j].n > P.MaxPopDensForInitialInfection)
					|| (Mcells[j].n < P.MinPopDensForInitialInfection) || (((int)(AdUnits[Mcells[j].adunit].id / P.CountryDivisor) != P.TargetCountry) && (P.TargetCountry >= 0)) // from Mcells[j].country to Mcells[j].adunit/P.CountryDivisor
					|| ((P.InitialInfectionsAdminUnit[i] > 0) && ((AdUnits[Mcells[j].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor != (P.InitialInfectionsAdminUnit[i] % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor)));
				for (k = 0; (k < nsi[i]) && (m < 10000); k++)
				{
					l = Mcells[j].members[(int)(ranf() * ((double)Mcells[j].n))];
					if (Hosts[l].inf == 0)
					{
						if (CalcPersonSusc(l, 0, 0, 0) > 0)
						{
							P.LocationInitialInfection[i][0] = Households[Hosts[l].hh].loc_x;
							P.LocationInitialInfection[i][1] = Households[Hosts[l].hh].loc_y;
							Hosts[l].infector = -2; Hosts[l].infect_type = INFECT_TYPE_MASK - 1;
							Hosts[l].base_inf_level = 1.0;
							Hosts[l].resist = CalcSeedResist();
							DoInfect(l, t, 0, run);
							m = 0;
						}
					}
					else
					{
						k--; m++;
					}
				}
				if (m)
					f++;
				else
					f = 0;
			} while ((f > 0) && (f < 1000));
		}
		else
		{
			m = 0;
			for (k = 0; (k < nsi[i]) && (m < 10000); k++)
			{
				do
				{
					l = (int)(ranf() * ((double)P.N));
					j = Hosts[l].mcell;
					//fprintf(stderr,"%i ",AdUnits[Mcells[j].adunit].id);
				} while ((Mcells[j].n == 0) || (Mcells[j].n > P.MaxPopDensForInitialInfection)
					|| (Mcells[j].n < P.MinPopDensForInitialInfection) || (((int)(AdUnits[Mcells[j].adunit].id / P.CountryDivisor) != P.TargetCountry) && (P.TargetCountry >= 0)) // from Mcells[j].country to Mcells[j].adunit/P.CountryDivisor
					|| ((P.InitialInfectionsAdminUnit[i] > 0) && ((AdUnits[Mcells[j].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor != (P.InitialInfectionsAdminUnit[i] % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor)));
				l = Mcells[j].members[(int)(ranf() * ((double)Mcells[j].n))];
				if (Hosts[l].inf == 0)
				{
					if (CalcPersonSusc(l, 0, 0, 0) > 0)
					{
						P.LocationInitialInfection[i][0] = Households[Hosts[l].hh].loc_x;
						P.LocationInitialInfection[i][1] = Households[Hosts[l].hh].loc_y;
						Hosts[l].infector = -2; Hosts[l].infect_type = INFECT_TYPE_MASK - 1;
						Hosts[l].base_inf_level = 1.0;
						Hosts[l].resist = CalcSeedResist();
						DoInfect(l, t, 0, run);
						m = 0;
					}
					else
					{
						k--; m++;
					}
				}
				else
				{
					k--; m++;
				}
			}
		}
	}
	if (m > 0) fprintf(stderr, "### Seeding error ###\n");
}

int CalcSeedResist(void)
{
	int i, j;
	double t;

	if (MAX_NUM_RESIST_TYPES == 1)
		i = 0;
	else
	{
		t = ranf();
		for (i = 0; (i < MAX_NUM_RESIST_TYPES) && (t >= 0); i++)
		{
			t -= P.EvolResistSeedProp[i];
		}
		i--;
	}
	return i;
}

void LoadSnapshot(void)
{
	FILE* dat;
	int i, j, tsi, * CellMemberArray, * CellSuscMemberArray;
	long l;
	long long CM_offset, CSM_offset;
	double t, st;
	int** Array_InvCDF;
	float* Array_tot_prob, ** Array_cum_trans, ** Array_max_trans;

	if (!(dat = fopen(SnapshotLoadFile, "rb"))) ERR_CRITICAL("Unable to open snapshot file\n");
	fprintf(stderr, "Loading snapshot.");
	if (!(Array_InvCDF = (int**)malloc(P.NCP * sizeof(int*)))) ERR_CRITICAL("Unable to allocate temp cell storage\n");
	if (!(Array_max_trans = (float**)malloc(P.NCP * sizeof(float*)))) ERR_CRITICAL("Unable to temp allocate cell storage\n");
	if (!(Array_cum_trans = (float**)malloc(P.NCP * sizeof(float*)))) ERR_CRITICAL("Unable to temp allocate cell storage\n");
	if (!(Array_tot_prob = (float*)malloc(P.NCP * sizeof(float)))) ERR_CRITICAL("Unable to temp allocate cell storage\n");
	for (i = 0; i < P.NCP; i++)
	{
		Array_InvCDF[i] = Cells[i].InvCDF;
		Array_max_trans[i] = Cells[i].max_trans;
		Array_cum_trans[i] = Cells[i].cum_trans;
		Array_tot_prob[i] = Cells[i].tot_prob;
	}

	fread_big((void*)&i, sizeof(int), 1, dat); if (i != P.N) { fprintf(stderr, "Incorrect N (%i %i) in snapshot file.\n", P.N, i); exit(1); }
	fread_big((void*)&i, sizeof(int), 1, dat); if (i != P.NH) ERR_CRITICAL("Incorrect NH in snapshot file.\n");
	fread_big((void*)&i, sizeof(int), 1, dat); if (i != P.NC) { fprintf(stderr, "## %i neq %i\n", i, P.NC); ERR_CRITICAL("Incorrect NC in snapshot file.\n"); }
	fread_big((void*)&i, sizeof(int), 1, dat); if (i != P.NCP) ERR_CRITICAL("Incorrect NCP in snapshot file.\n");
	fread_big((void*)&i, sizeof(int), 1, dat); if (i != P.ncw) ERR_CRITICAL("Incorrect ncw in snapshot file.\n");
	fread_big((void*)&i, sizeof(int), 1, dat); if (i != P.nch) ERR_CRITICAL("Incorrect nch in snapshot file.\n");
	fread_big((void*)&l, sizeof(long), 1, dat); if (l != P.seed1) ERR_CRITICAL("Incorrect seed1 in snapshot file.\n");
	fread_big((void*)&l, sizeof(long), 1, dat); if (l != P.seed2) ERR_CRITICAL("Incorrect seed2 in snapshot file.\n");
	fread_big((void*)&t, sizeof(double), 1, dat); if (t != P.TimeStep) ERR_CRITICAL("Incorrect TimeStep in snapshot file.\n");
	fread_big((void*)&(P.SnapshotLoadTime), sizeof(double), 1, dat);
	P.NumSamples = 1 + (int)ceil((P.SampleTime - P.SnapshotLoadTime) / P.SampleStep);
	fprintf(stderr, ".");
	fread_big((void*)&CellMemberArray, sizeof(int*), 1, dat);
	fprintf(stderr, ".");
	fread_big((void*)&CellSuscMemberArray, sizeof(int*), 1, dat);
	fprintf(stderr, ".");
	CM_offset = State.CellMemberArray - CellMemberArray;
	CSM_offset = State.CellSuscMemberArray - CellSuscMemberArray;

	zfread_big((void*)Hosts, sizeof(person), (size_t)P.N, dat);
	fprintf(stderr, ".");
	zfread_big((void*)Households, sizeof(household), (size_t)P.NH, dat);
	fprintf(stderr, ".");
	zfread_big((void*)Cells, sizeof(cell), (size_t)P.NC, dat);
	fprintf(stderr, ".");
	zfread_big((void*)Mcells, sizeof(microcell), (size_t)P.NMC, dat);
	fprintf(stderr, ".");
	zfread_big((void*)State.CellMemberArray, sizeof(int), (size_t)P.N, dat);
	fprintf(stderr, ".");
	zfread_big((void*)State.CellSuscMemberArray, sizeof(int), (size_t)P.N, dat);
	fprintf(stderr, ".");
	for (i = 0; i < P.NC; i++)
	{
		if (Cells[i].n > 0)
		{
			Cells[i].members += CM_offset;
			Cells[i].susceptible += CSM_offset;
			Cells[i].latent += CSM_offset;
			Cells[i].infected += CSM_offset;
		}
		for (j = 0; j < MAX_INTERVENTION_TYPES; j++) Cells[i].CurInterv[j] = -1; // turn interventions off in loaded image
	}
	for (i = 0; i < P.NMC; i++)
		if (Mcells[i].n > 0)
			Mcells[i].members += CM_offset;

	for (i = 0; i < P.NCP; i++)
	{
		Cells[i].InvCDF = Array_InvCDF[i];
		Cells[i].max_trans = Array_max_trans[i];
		Cells[i].cum_trans = Array_cum_trans[i];
		Cells[i].tot_prob = Array_tot_prob[i];
	}
	free(Array_tot_prob);
	free(Array_cum_trans);
	free(Array_max_trans);
	free(Array_InvCDF);
	fprintf(stderr, "\n");
	fclose(dat);
}

void SaveSnapshot(void)
{
	FILE* dat;
	int i = 1, j, n, st;
	double sz;

	if (!(dat = fopen(SnapshotSaveFile, "wb"))) ERR_CRITICAL("Unable to open snapshot file\n");

	fwrite_big((void*)&(P.N), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)&(P.NH), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)&(P.NC), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)&(P.NCP), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)&(P.ncw), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)&(P.nch), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)&(P.seed1), sizeof(long), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)&(P.seed2), sizeof(long), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)&(P.TimeStep), sizeof(double), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)&(P.SnapshotSaveTime), sizeof(double), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)&(State.CellMemberArray), sizeof(int*), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)&(State.CellSuscMemberArray), sizeof(int*), 1, dat);
	fprintf(stderr, "## %i\n", i++);

	zfwrite_big((void*)Hosts, sizeof(person), (size_t)P.N, dat);

	fprintf(stderr, "## %i\n", i++);
	zfwrite_big((void*)Households, sizeof(household), (size_t)P.NH, dat);
	fprintf(stderr, "## %i\n", i++);
	zfwrite_big((void*)Cells, sizeof(cell), (size_t)P.NC, dat);
	fprintf(stderr, "## %i\n", i++);
	zfwrite_big((void*)Mcells, sizeof(microcell), (size_t)P.NMC, dat);
	fprintf(stderr, "## %i\n", i++);

	zfwrite_big((void*)State.CellMemberArray, sizeof(int), (size_t)P.N, dat);
	fprintf(stderr, "## %i\n", i++);
	zfwrite_big((void*)State.CellSuscMemberArray, sizeof(int), (size_t)P.N, dat);
	fprintf(stderr, "## %i\n", i++);

	fclose(dat);
}
