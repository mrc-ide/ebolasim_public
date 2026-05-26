/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "model_constants.h"
#include "model_structs.h"
#include "model_fns.h"
#include "model_globals_ext.h"
#include "binio.h"



void DoInfect(int ai, double t, int tn, int run) //added int as argument to DoInfect to record run number: ggilani - 15/10/14
{

	int i, j;
	unsigned short int ts;
	double q, x, y;
	person* a;

	a = Hosts + ai;

	if (a->inf == 0)
	{
		ts = (unsigned short int) (P.TimeStepsPerDay * t);
		a->inf = 1;
		a->infection_time = (unsigned short int) ts;
		if (a->infector >= 0)
		{
			a->resist = Hosts[a->infector].resist;
			a->base_inf_level = Hosts[a->infector].base_inf_level;
			if (abs(Hosts[a->infector].inf) == 6) //if host has an infector (i.e. they are not a seed case), check to see if infection status flag of host is 6, i.e. that it is a funeral contact, ggilani: 24/10/14
			{
				StateT[tn].cumFI++;
			}
		}
		StateT[tn].cumI++;
		StateT[tn].cumItype[a->infect_type % INFECT_TYPE_MASK]++;
		StateT[tn].cumIa[HOST_AGE_GROUP(ai)]++;
		x = (Households[a->hh].loc_x - P.LocationInitialInfection[0][0]);
		y = (Households[a->hh].loc_y - P.LocationInitialInfection[0][1]);
		q = x * x + y * y;
		StateT[tn].sumRad2 += q;
		if (q > StateT[tn].maxRad2) StateT[tn].maxRad2 = q;
		{
			Cells[a->pcell].S--;
			if (Cells[a->pcell].S > 0)
			{
				Cells[a->pcell].susceptible[a->listpos] = Cells[a->pcell].susceptible[Cells[a->pcell].S];
				Hosts[Cells[a->pcell].susceptible[a->listpos]].listpos = a->listpos;
			}
			Cells[a->pcell].latent--;
			Cells[a->pcell].L++;
			a->listpos = Cells[a->pcell].S;
			Cells[a->pcell].latent[0] = ai;
		}
		if ((HOST_TREATED(ai)) && (a->resist < (MAX_NUM_RESIST_TYPES - 1)) && (ranf_mt(tn) < P.EvolResistProphMutationRate)) a->resist++;
		StateT[tn].cumI_resist[a->resist]++;
		StateT[tn].cumI_keyworker[a->keyworker]++;
		if (P.DoLatent)
		{
			i = (int)floor((q = ranf_mt(tn) * CDF_RES));
			q -= ((double)i);
			a->latent_time = (unsigned short int) floor(0.5 + (t - P.LatentPeriod * log(q * P.latent_icdf[i + 1] + (1.0 - q) * P.latent_icdf[i])) * P.TimeStepsPerDay);
		}
		else
			a->latent_time = (unsigned short int) (t * P.TimeStepsPerDay);
		if (P.DoAdUnits)
		{
			StateT[tn].cumI_adunit[Mcells[a->mcell].adunit]++;
		}
		////Add some case detection code here which determines whether someone will be detected as soon has they become infected - however, we still have some
		////case detection code in DoCase, to ensure that they are added to the detected cases at the right time.
		////This is to allow us to output bitmaps with only detected cases - ggilani 04/08/2015
		//if(P.DoCaseDetection)
		//{
		//	if((P.DoCaseDetectionAdunit)&&(P.DoAdUnits))
		//	{
		//		if((Hosts[ai].rep_rate<AdUnits[Mcells[a->mcell].adunit].caseDetectRate)||(Hosts[ai].contactTraced>0)) //add something to ensure that hosts who are being currently contact traced are always detected, changed ranf_mt(tn) to Hosts[ai].rep_rate 13/06/19
		//		{
		//			Hosts[ai].detected=1;
		//		}
		//	}
		//	else
		//	{
		//		if(Hosts[ai].rep_rate<P.CaseDetectionRate||(Hosts[ai].contactTraced>0))
		//		{
		//			Hosts[ai].detected=1;
		//		}
		//	}
		//	
		//}
		if (P.OutputBitmap)
		{
			if ((P.OutputBitmapDetected == 0) || ((P.OutputBitmapDetected == 1) && (Hosts[ai].detected == 1)))
			{
				x = ((int)(Households[a->hh].loc_x * P.scalex)) - P.bminx;
				y = ((int)(Households[a->hh].loc_y * P.scaley)) - P.bminy;
				if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
				{
					j = y * bmh->width + x;
					if ((j < bmh->imagesize) && (j >= 0))
					{
#pragma omp atomic
						bmi2[j]++;
					}
				}
			}
		}
		//added this to record event if flag is set to 1 : ggilani - 10/10/2014
		if (P.DoRecordInfEvents)
		{
			if (*nEvents < P.MaxInfEvents)
			{
				RecordEvent(t, ai, run, 0, tn); //added int as argument to RecordEvent to record run number: ggilani - 15/10/14
			}
		}
		if ((t > 0) && (P.DoOneGen))
		{
			DoIncub(ai, ts, tn, run);
			DoCase(ai, t, ts, tn);
			DoRecover(ai, run, tn);
		}
	}
}

void DoIncub(int ai, unsigned short int ts, int tn, int run)
{
	int i;
	person* a;
	double q, ti, cfr;
	int age, day;

	age = HOST_AGE_GROUP(ai);
	if (age >= NUM_AGE_GROUPS) age = NUM_AGE_GROUPS - 1;

	if (!P.DoEventMortality)
	{
		if (P.DoAgeMortality)
		{
			cfr = P.AgeMortality[age];
		}
		else
		{
			cfr = P.DiseaseMortalityVacc;
		}

		if (HOST_TO_BE_VACCED(ai) || HOST_VACCED(ai))
		{
			cfr = cfr * P.DiseaseMortalityVacc;
		}
	}

	a = Hosts + ai;
	if (a->inf == 1)
	{
		if (P.InfectiousnessSD == 0)
			a->infectiousness = P.AgeInfectiousness[age];
		else
			//a->infectiousness=P.AgeInfectiousness[age]*gen_beta_mt(P.InfectiousnessBetaA,P.InfectiousnessBetaB,tn);
			a->infectiousness = P.AgeInfectiousness[age] * gen_gamma_mt(P.InfectiousnessGamA, P.InfectiousnessGamR, tn);
		q = P.ProportionSymptomatic[age]
			* (HOST_TREATED(ai) ? P.EvolResistRelTreatSympDrop[Hosts[ai].resist] : 1)
			* (HOST_VACCED(ai) ? (1 - P.VaccSympDrop) : 1);
		if (ranf_mt(tn) < q)
		{
			a->inf = -1;
			a->infectiousness = -P.SymptInfectiousness * a->infectiousness;
		}
		else
		{
			a->inf = 2;
		}
		// commenting this out for now to redo cfr/mortality taking into account the impact of vaccination

		if (P.DoInfectiousnessProfile)
			a->recovery_time = a->latent_time + (unsigned short int) (P.InfectiousPeriod * P.TimeStepsPerDay);
		else
		{
			i = (int)floor(q = ranf_mt(tn) * CDF_RES);
			q -= ((double)i);
			ti = -P.InfectiousPeriod * log(q * P.infectious_icdf[i + 1] + (1.0 - q) * P.infectious_icdf[i]);
			a->recovery_time = a->latent_time + (unsigned short int) floor(0.5 + (ti * P.TimeStepsPerDay));
		}

		if (P.DoMortality) //added DoMortality to allow different approaches to assigning mortality: ggilani - 22/10/2014
		{
			if (P.DoEventMortality)
			{
				day = (int)ceil(P.TimeStep * (a->recovery_time - a->latent_time));
				if (ranf_mt(tn) <= P.RecoveryProb[day])
				{
					a->to_die = 0;
				}
				else
				{
					a->to_die = 1;
				}
			}
			else
			{

				if (ranf_mt(tn) < cfr)
				{
					a->recovery_time = (a->latent_time + (unsigned short int) (P.LethalInfectiousPeriod * ((double)(a->recovery_time - a->latent_time)))); //put LethalInfectiousPeriod as proportion of infectious period - otherwise time until death is too long
					a->to_die = 1;
				}
				else
				{
					a->to_die = 0;
				}
			}
		}


		////added time to hospitalisation based on onset to hospitalisation distribution, if we are considering hospitalisation: ggilani 28/10/14
		//if((P.DoHospitalisation)&(a->detected))
		//{
		//	//i=(int) floor(q=ranf_mt(tn)*CDF_RES);
		//	//q-=((double) i);
		//	//ti=-P.HospitalisationTime*log(q*P.hospital_icdf[i+1]+(1.0-q)*P.hospital_icdf[i]);
		//	//ti=P.HospitalisationTime*sexpo_mt(tn);
		//	if(a->contactTraced==0)
		//	{
		//		a->hospital_time=a->latent_time+(unsigned short int) floor(0.5+(P.HospitalisationTime*P.TimeStepsPerDay));
		//	}
		//	else
		//	{
		//		a->hospital_time=a->latent_time+(unsigned short int) floor(0.5+(P.HospitalisationTime_contactTrace*P.TimeStepsPerDay)); //different hospitalisation time for contact traced case: ggilani 05/07/2017
		//	}
		//}
		Cells[a->pcell].L--;
		if (Cells[a->pcell].L > 0)
		{
			Cells[a->pcell].susceptible[a->listpos] = Cells[a->pcell].latent[Cells[a->pcell].L];
			Hosts[Cells[a->pcell].susceptible[a->listpos]].listpos = a->listpos;
		}
		Cells[a->pcell].infected--;
		a->listpos = Cells[a->pcell].S + Cells[a->pcell].L;
		Cells[a->pcell].I++;
		Cells[a->pcell].infected[0] = ai;
		////added this to record event if flag is set to 1 and if host isn't initial seed, i.e. if Hosts[ai].infector>=0: ggilani - 10/10/2014
		//if(P.DoRecordInfEvents)
		//	{
		//		if(*nEvents<P.MaxInfEvents)
		//		{
		//			RecordEvent(((double)a->latent_time)/P.TimeStepsPerDay,ai,run,1); //added int as argument to RecordEvent to record run number: ggilani - 15/10/14
		//		}
		//	}
	}
}

void DoImmune(int ai)
// This transfers a person straight from susceptible to immune. Used to start a run with a partially immune population.
{
	person* a;
	int j, c;
	double x, y;

	a = Hosts + ai;
	if (a->inf == 0)
	{
		c = a->pcell;
		a->inf = 4;
		Cells[c].S--;
		if (a->listpos < Cells[c].S)
		{
			Cells[c].susceptible[a->listpos] = Cells[c].susceptible[Cells[c].S];
			Hosts[Cells[c].susceptible[a->listpos]].listpos = a->listpos;
		}
		if (Cells[c].L > 0)
		{
			Cells[c].susceptible[Cells[c].S] = Cells[c].susceptible[Cells[c].S + Cells[c].L];
			Hosts[Cells[c].susceptible[Cells[c].S]].listpos = Cells[c].S;
		}
		if (Cells[c].I > 0)
		{
			Cells[c].susceptible[Cells[c].S + Cells[c].L] = Cells[c].susceptible[Cells[c].S + Cells[c].L + Cells[c].I];
			Hosts[Cells[c].susceptible[Cells[c].S + Cells[c].L]].listpos = Cells[c].S + Cells[c].L;
		}
		if (a->listpos < Cells[c].S + Cells[c].L + Cells[c].I)
		{
			Cells[c].susceptible[Cells[c].S + Cells[c].L + Cells[c].I] = ai;
			a->listpos = Cells[c].S + Cells[c].L + Cells[c].I;
		}
		Cells[c].latent--;
		Cells[c].infected--;
		Cells[c].R++;
		if (P.OutputBitmap)
		{
			x = ((int)(Households[a->hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[a->hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				j = y * bmh->width + x;
				if ((j < bmh->imagesize) && (j >= 0))
				{
#pragma omp atomic
					bmi3[j]++;
				}
			}
		}
	}
}


void DoDetectedCase(int ai, double t, unsigned short int ts, int tn)
{

	int j, k, f, j1, j2, age;
	int i, ii, i1, i2, i3, k2, k3, l, vaccflag;
	int currentRing, cnt, nActiveCases;
	person* a;

	a = Hosts + ai;
	age = HOST_AGE_GROUP(ai);
	if ((!P.RestrictTreatToTarget) || (Mcells[a->mcell].country == P.TargetCountry))
	{
		//only increase triggers if outbreak has been declared
		//if (P.OutbreakDetected)
		//{
		if (Mcells[a->mcell].treat_trig < USHRT_MAX - 1) Mcells[a->mcell].treat_trig++;
		if ((!P.OnlyDoGeoVaccWhenNoRing) || (Hosts[ai].vacc_accept > P.ProbEstablishRing))
		{
			if (Mcells[a->mcell].vacc_trig < USHRT_MAX - 1) Mcells[a->mcell].vacc_trig++;
		}
		if (Mcells[a->mcell].move_trig < USHRT_MAX - 1) Mcells[a->mcell].move_trig++;
		if (Mcells[a->mcell].socdist_trig < USHRT_MAX - 1) Mcells[a->mcell].socdist_trig++;
		if (Mcells[a->mcell].keyworkerproph_trig < USHRT_MAX - 1) Mcells[a->mcell].keyworkerproph_trig++;
		//}
#ifndef ABSENTEEISM_PLACE_CLOSURE
#ifdef PLACE_CLOSE_ROUND_HOUSEHOLD
		if (Mcells[a->mcell].place_trig < USHRT_MAX - 1) Mcells[a->mcell].place_trig++;
#endif
		if (t >= P.PlaceCloseTimeStart)
			for (j = 0; j < P.PlaceTypeNum; j++)
				if ((j != HOTEL_PLACE_TYPE) && (a->PlaceLinks[j] >= 0))
				{
					if ((!P.RestrictTreatToTarget) || (Places[j][a->PlaceLinks[j]].country == P.TargetCountry))
					{
						DoPlaceClose(j, a->PlaceLinks[j], ts, tn, 0);
#ifndef PLACE_CLOSE_ROUND_HOUSEHOLD
						if (Mcells[Places[j][a->PlaceLinks[j]].mcell].place_trig < USHRT_MAX - 1)
						{
#pragma omp critical (place_trig)
							Mcells[Places[j][a->PlaceLinks[j]].mcell].place_trig++;
						}
#endif
					}
				}
#endif
		if (t >= P.TreatTimeStart)
		{
			if ((P.DoHouseholds) && (Households[Hosts[ai].hh].stockpile == 1))
			{
				if ((P.PrivateTreatPropCases == 1) || (ranf_mt(tn) < P.PrivateTreatPropCases))
				{
					DoPrivateTreatCase(ai, ts, tn);
					if ((t < P.TreatTimeStart + P.TreatHouseholdsDuration) && ((P.PrivateTreatPropCaseHouseholds == 1) || (ranf_mt(tn) < P.PrivateTreatPropCaseHouseholds)))
					{
						j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
						for (j = j1; j < j2; j++)
							if (!HOST_TO_BE_TREATED(j)) DoPrivateProph(j, ts, tn);
					}
				}
			}
			else if ((P.TreatPropCases == 1) || (ranf_mt(tn) < P.TreatPropCases))
			{
				DoTreatCase(ai, ts, tn);
				if (P.DoHouseholds)
				{
					if ((t < P.TreatTimeStart + P.TreatHouseholdsDuration) && ((P.TreatPropCaseHouseholds == 1) || (ranf_mt(tn) < P.TreatPropCaseHouseholds)))
					{
						j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
						for (j = j1; j < j2; j++)
							if (!HOST_TO_BE_TREATED(j)) DoProph(j, ts, tn);
					}
				}
				if (P.DoPlaces)
				{
					if (t < P.TreatTimeStart + P.TreatPlaceGeogDuration)
						for (j = 0; j < P.PlaceTypeNum; j++)
							if (a->PlaceLinks[j] >= 0)
							{
								if ((!P.RestrictTreatToTarget) || (Places[j][a->PlaceLinks[j]].country == P.TargetCountry))
								{
									if (P.DoPlaceGroupTreat)
									{
										if ((P.TreatPlaceProbCaseId[j] == 1) || (ranf_mt(tn) < P.TreatPlaceProbCaseId[j]))
										{
											StateT[tn].p_queue[j][StateT[tn].np_queue[j]] = a->PlaceLinks[j];
											StateT[tn].pg_queue[j][StateT[tn].np_queue[j]++] = a->PlaceGroupLinks[j];
										}
									}
									else
									{
										f = 0;
#pragma omp critical (starttreat)
										if (!Places[j][a->PlaceLinks[j]].treat) f = Places[j][a->PlaceLinks[j]].treat = 1;
										if (f)
										{
											if ((P.TreatPlaceProbCaseId[j] == 1) || (ranf_mt(tn) < P.TreatPlaceProbCaseId[j]))
												StateT[tn].p_queue[j][StateT[tn].np_queue[j]++] = a->PlaceLinks[j];
											else
												Places[j][a->PlaceLinks[j]].treat = 0;
										}
									}
								}
							}
				}
			}
		}
		if (P.DoHouseholds)
		{
			if ((!P.DoMassVacc) && (t >= P.VaccTimeStart) && (State.cumV < P.VaccMaxCourses))
			{
				if ((t < P.VaccTimeStart + P.VaccHouseholdsDuration) && ((P.VaccPropCaseHouseholds == 1) || (Hosts[ai].vacc_accept < P.VaccPropCaseHouseholds)))
				{
					j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
					for (j = j1; j < j2; j++)
					{
						StateT[tn].vacc_queue[StateT[tn].nvacc_queue] = j;
						StateT[tn].nvacc_queue++;
						//DoVacc(j, ts);
					}
				}

			}
			if ((t >= P.HQuarantineTimeStart) && (t < P.HQuarantineTimeStart + P.HQuarantinePolicyDuration))
			{
				j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
				if (!HOST_TO_BE_QUARANTINED(j1))
				{
					Hosts[j1].quar_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.HQuarantineHouseDelay));
					k = (ranf_mt(tn) < P.HQuarantinePropHouseCompliant) ? 1 : 0;
					if (k) StateT[tn].cumHQ++;

					Hosts[j1].quar_comply = ((k == 0) ? 0 : ((ranf_mt(tn) < P.HQuarantinePropIndivCompliant) ? 1 : 0));
					if ((Hosts[j1].quar_comply) && (!HOST_ABSENT(j1)))
					{
						if (HOST_AGE_YEAR(j1) >= P.CaseAbsentChildAgeCutoff)
						{
							if (Hosts[j1].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] >= 0) StateT[tn].cumAH++;
						}
						else
							StateT[tn].cumACS++;
					}
					for (j = j1 + 1; j < j2; j++)
					{
						Hosts[j].quar_start_time = Hosts[j1].quar_start_time;
						Hosts[j].quar_comply = ((k == 0) ? 0 : ((ranf_mt(tn) < P.HQuarantinePropIndivCompliant) ? 1 : 0));
						if ((Hosts[j].quar_comply) && (!HOST_ABSENT(j)))
						{
							if (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff)
							{
								if (Hosts[j].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] >= 0) StateT[tn].cumAH++;
							}
							else
								StateT[tn].cumACS++;
						}
					}
				}
			}
		}
		if ((t >= P.CaseIsolationTimeStart) && (t < P.CaseIsolationTimeStart + P.CaseIsolationPolicyDuration))
		{
			if ((P.CaseIsolationProp == 1) || (ranf_mt(tn) < P.CaseIsolationProp))
			{
				Hosts[ai].isolation_start_time = ts;
				if (HOST_ABSENT(ai))
				{
					if (a->absent_stop_time < ts + P.usCaseAbsenteeismDelay + P.usCaseIsolationDuration)
						a->absent_stop_time = ts + P.usCaseAbsenteeismDelay + P.usCaseIsolationDuration;
				}
				else if (P.DoRealSymptWithdrawal)
					/* This calculates adult absenteeism from work due to care of isolated children.  */
				{
					Hosts[ai].absent_start_time = ts + P.usCaseIsolationDelay;
					Hosts[ai].absent_stop_time = ts + P.usCaseIsolationDelay + P.usCaseIsolationDuration;
					if (P.DoPlaces)
					{
						if ((!HOST_QUARANTINED(ai)) && (Hosts[ai].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] >= 0) && (HOST_AGE_YEAR(ai) >= P.CaseAbsentChildAgeCutoff))
							StateT[tn].cumAC++;
					}
					if ((P.DoHouseholds) && (P.DoPlaces) && (HOST_AGE_YEAR(ai) < P.CaseAbsentChildAgeCutoff))
					{
						if (!HOST_QUARANTINED(ai)) StateT[tn].cumACS++;
						if ((P.CaseAbsentChildPropAdultCarers == 1) || (ranf_mt(tn) < P.CaseAbsentChildPropAdultCarers))
						{
							j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
							f = 0;
							for (j = j1; (j < j2) && (!f); j++)
								f = ((abs(Hosts[j].inf) != 5) && (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && (Hosts[j].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] < 0));
							if (!f)
							{
								for (j = j1; (j < j2) && (!f); j++)
									f = ((abs(Hosts[j].inf) != 5) && (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && ((HOST_ABSENT(j)) || (HOST_QUARANTINED(j))));
								if (!f)
								{
									f = 0; k = -1;
									for (j = j1; (j < j2) & (!f); j++)
										if ((HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && (abs(Hosts[j].inf) != 5)) { k = j; f = 1; }
									if (f)
									{
										Hosts[k].absent_start_time = ts + P.usCaseIsolationDelay;
										Hosts[k].absent_stop_time = ts + P.usCaseIsolationDelay + P.usCaseIsolationDuration;
										StateT[tn].cumAA++;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	StateT[tn].cumDC++;
	if (P.DoAdUnits)
	{
		StateT[tn].cumDC_adunit[Mcells[a->mcell].adunit]++;
		AdUnits[Mcells[a->mcell].adunit].lastCaseDay = t;
	}
	StateT[tn].cumDCa[age]++; //added detected case by age: ggilani 22/02/22
	Hosts[ai].dayDetected = (int)t;
	if (Hosts[ai].contactTraced)
	{
		StateT[tn].cumCC++;
		if (P.DoAdUnits)
		{
			StateT[tn].cumCC_adunit[Mcells[a->mcell].adunit]++;
		}
		Hosts[ai].contactTraced_end_time = ts;
	}



	//if ((P.DoCaseDetectionAdunit) && (P.DoAdUnits))
	//{
	//		//we already know whether the case will be detected based on some code in DoInfect
	//	if (Hosts[ai].detected == 1)
	//	{
	//		if (((!P.RestrictTreatToTarget) || (Mcells[a->mcell].country == P.TargetCountry)) && (t > P.VaccTimeStartGeo))
	//		{
	//			if (Mcells[Hosts[ai].mcell].vacc_trig < USHRT_MAX - 1) Mcells[Hosts[ai].mcell].vacc_trig++;
	//		}
	//		StateT[tn].cumDC++;
	//		StateT[tn].cumDC_adunit[Mcells[a->mcell].adunit]++;
	//		StateT[tn].cumDCa[age]++; //added detected case by age: ggilani 22/02/22
	//		//store day on which host is detected
	//		Hosts[ai].dayDetected = (int)t;
	//		if (Hosts[ai].contactTraced == 1)
	//		{
	//			StateT[tn].cumCC_adunit[Mcells[a->mcell].adunit]++;
	//			StateT[tn].cumCC++;
	//		}
	//	}
	//}
	//else
	//{
	//	if (Hosts[ai].detected == 1)
	//	{
	//		StateT[tn].cumDC++;
	//		StateT[tn].cumDCa[age]++;
	//		if (P.DoAdUnits)
	//		{
	//			StateT[tn].cumDC_adunit[Mcells[a->mcell].adunit]++;
	//		}
	//		if (Hosts[ai].contactTraced == 1)
	//		{
	//			StateT[tn].cumCC++;
	//			if (P.DoAdUnits)
	//			{
	//				StateT[tn].cumCC_adunit[Mcells[a->mcell].adunit]++;
	//			}
	//		}
	//	}
	//}

		//Add some code to take care of isolation of cases - ggilani 09/06/2017
		//If someone has been contact traced and they become a case (and are detected), then they will be isolated
		//if(Hosts[ai].contactTraced>0)
		//{
		//	Hosts[ai].isolation_start_time=ts;
		//}
		//Else detected cases are isolated with probability CaseIsolationProp
		//else if((P.CaseIsolationProp==1)||(ranf_mt(tn)<P.CaseIsolationProp))
		//{
		//	Hosts[ai].isolation_start_time=ts;
		//}


	//Add my ring vaccination code here similar to DoDetectedCase code - ggilani 15/02/17
	if ((Hosts[ai].detected) && (Hosts[ai].vacc_accept < P.ProbEstablishRing) && ((P.DoRingVaccination && (t >= P.VaccTimeStart)) || (P.DoContactTracing && (t >= P.ContactTracingTimeStart))))
	{
		nActiveCases = State.cumDC_adunit[Mcells[Hosts[ai].mcell].adunit];// -(State.cumDD_adunit[Mcells[Hosts[ai].mcell].adunit] + State.cumDR_adunit[Mcells[Hosts[ai].mcell].adunit]);
		//first check to see if we've reached the threshold to start contact tracing in the admin unit of the initial case
		if ((t >= P.ContactTracingTimeStart) && (nActiveCases >= AdUnits[Mcells[Hosts[ai].mcell].adunit].contactTraceCaseThreshold) && (AdUnits[Mcells[Hosts[ai].mcell].adunit].contactTraceThresholdCrossed == 0)) //changed to detected cases
		{
			//mark that we've now reached the threshold and mark the day on which it occurs
			AdUnits[Mcells[Hosts[ai].mcell].adunit].contactTraceThresholdCrossed = 1;
			AdUnits[Mcells[Hosts[ai].mcell].adunit].contactTraceStartDay = (int)t;
		}

		if (((t >= P.VaccTimeStart) && (State.cumV < P.VaccMaxCourses)) || AdUnits[Mcells[Hosts[ai].mcell].adunit].contactTraceThresholdCrossed == 1) //modified this to include criteria for having to be past contact tracing threshold
		{

			if (P.DoContactTracing && (AdUnits[Mcells[Hosts[ai].mcell].adunit].contactTraceThresholdCrossed == 1))
			{
				//set current ring to 1
				currentRing = 1;

				//first set counting variable to zero, then add household members
				cnt = 0;
				i1 = Households[Hosts[ai].hh].FirstPerson;
				i2 = i1 + Households[Hosts[ai].hh].nh;
				for (i = i1; i < i2; i++)
				{
					if ((i != ai))//&&(Hosts[i].contactTraced==0))//&&(Hosts[i].vaccRing==0))
					{
						StateT[tn].ringvacclist[cnt] = i;
						//Hosts[i].vaccRing = currentRing; //These people belong to ring 1
						//Hosts[i].ringCase = ai; //associated with the ring around case ai
						//StateT[tn].ringlist[cnt] = 1; //These people belong to ring 1
						cnt++;
					}
				}
				//add place members - find placetype that case belongs to
				for (j = 0; j < P.PlaceTypeNum; j++)
				{
					if ((Hosts[ai].PlaceLinks[j] >= 0) && (j != P.HospPlaceTypeNum))
					{
						for (k = 0; k < Places[j][Hosts[ai].PlaceLinks[j]].n; k++)
						{
							if ((Places[j][Hosts[ai].PlaceLinks[j]].members[k] != ai))//&&(Hosts[Places[j][Hosts[ai].PlaceLinks[j]].members[k]].contactTraced==0))
							{
								StateT[tn].ringvacclist[cnt] = Places[j][Hosts[ai].PlaceLinks[j]].members[k];
								//Hosts[Places[j][Hosts[ai].PlaceLinks[j]].members[k]].vaccRing = currentRing; //These people belong to ring 1
								//Hosts[Places[j][Hosts[ai].PlaceLinks[j]].members[k]].ringCase = ai; //associated with the ring around case ai
								cnt++;
							}
						}
					}
				}

				for (i = 0; i < cnt; i++)
				{
					//if (Hosts[StateT[tn].ringvacclist[i]].vaccRing == 1) //only contact tracing immediate contacts in the first ring
					//{
					if ((P.propContactTraced == 0) || (ranf_mt(tn) < P.propContactTraced))
					{
						StateT[tn].ct_queue[Mcells[Hosts[StateT[tn].ringvacclist[i]].mcell].adunit][StateT[tn].nct_queue[Mcells[Hosts[StateT[tn].ringvacclist[i]].mcell].adunit]++] = StateT[tn].ringvacclist[i];
					}
					//}
				}
			}

			if ((State.cumV < P.VaccMaxCourses) && (t >= P.VaccTimeStart))
			{
				//set current ring to 1
				currentRing = 1;
				//first set counting variable to zero, then add household members
				cnt = 0;

				//add some code here to revaccinate HCWs/FLWs in the admin unit if appropriate
				if ((P.RevaccHCWs) && (!AdUnits[Mcells[Hosts[ai].mcell].adunit].revacc))
				{
					//set this flag to 1 so that we don't check the same health workers again
					AdUnits[Mcells[Hosts[ai].mcell].adunit].revacc = 1;
					for (i = 0; i < P.Nplace[P.HospPlaceTypeNum]; i++)
					{
						if (Mcells[Places[P.HospPlaceTypeNum][i].mcell].adunit == Mcells[Hosts[ai].mcell].adunit)
						{
							//revacc this place
							//start with hcws
							for (j = 0; j < Places[P.HospPlaceTypeNum][i].nhcws; j++)
							{
								StateT[tn].ringvacclist[cnt] = Places[P.HospPlaceTypeNum][i].members[j];
								cnt++;
							}
							//now for front line workers
							for (j = 0; j < Places[P.HospPlaceTypeNum][i].nflws; j++)
							{
								StateT[tn].ringvacclist[cnt] = Places[P.HospPlaceTypeNum][i].flwmembers[j];
								cnt++;
							}
						}
					}
				}

				i1 = Households[Hosts[ai].hh].FirstPerson;
				i2 = i1 + Households[Hosts[ai].hh].nh;
				for (i = i1; i < i2; i++)
				{
					vaccflag = 1;
					if ((i != ai) && !HOST_TO_BE_VACCED(i))//&&(Hosts[i].vaccRing==0))
					{
						ii = 0;
						while ((ii < cnt) && (vaccflag == 1))
						{
							if (StateT[tn].ringvacclist[ii] == i)
							{
								vaccflag = 0;
							}
							ii++;
						}
						if (vaccflag)
						{
							StateT[tn].ringvacclist[cnt] = i;
							Hosts[i].vaccRing = currentRing; //These people belong to ring 1
							Hosts[i].ringCase = ai; //associated with the ring around case ai
							//StateT[tn].ringlist[cnt] = 1; //These people belong to ring 1
							cnt++;
						}

					}
				}
				//add place members - find placetype that case belongs to
				for (j = 0; j < P.PlaceTypeNum; j++)
				{
					if ((Hosts[ai].PlaceLinks[j] >= 0) && (j != P.HospPlaceTypeNum))
					{
						for (k = 0; k < Places[j][Hosts[ai].PlaceLinks[j]].n; k++)
						{
							vaccflag = 1;
							if ((Places[j][Hosts[ai].PlaceLinks[j]].members[k] != ai) && (!HOST_TO_BE_VACCED(Places[j][Hosts[ai].PlaceLinks[j]].members[k])))//&&(Hosts[Places[j][Hosts[ai].PlaceLinks[j]].members[k]].vaccRing==0))
							{
								ii = 0;
								while ((ii < cnt) && (vaccflag == 1))
								{
									if (StateT[tn].ringvacclist[ii] == Places[j][Hosts[ai].PlaceLinks[j]].members[k])
									{
										vaccflag = 0;
									}
									ii++;
								}
								if (vaccflag)
								{
									StateT[tn].ringvacclist[cnt] = Places[j][Hosts[ai].PlaceLinks[j]].members[k];
									Hosts[Places[j][Hosts[ai].PlaceLinks[j]].members[k]].vaccRing = currentRing; //These people belong to ring 1
									Hosts[Places[j][Hosts[ai].PlaceLinks[j]].members[k]].ringCase = ai; //associated with the ring around case ai
									cnt++;
								}
							}
						}
					}
				}

				//if number of rings is greater than 1, add extra members to the list
				if (P.NVaccRingsActive > 1)
				{
					//increment currentRing
					currentRing = 2;
					//First add group members of household contacts
					i1 = Households[Hosts[ai].hh].FirstPerson;
					i2 = i1 + Households[Hosts[ai].hh].nh;
					for (i = i1; i < i2; i++)
					{
						for (j = 0; j < P.PlaceTypeNum; j++)
						{
							if ((Hosts[i].PlaceLinks[j] >= 0) && (j != P.HospPlaceTypeNum))
							{
								for (k = 0; k < Places[j][Hosts[i].PlaceLinks[j]].n; k++)
								{
									vaccflag = 1;
									if ((Places[j][Hosts[i].PlaceLinks[j]].members[k] != ai) && (!HOST_TO_BE_VACCED(Places[j][Hosts[i].PlaceLinks[j]].members[k])))//&&(Hosts[Places[j][Hosts[i].PlaceLinks[j]].members[k]].vaccRing==0))
									{
										ii = 0;
										while ((ii < cnt) && (vaccflag == 1))
										{
											if (StateT[tn].ringvacclist[ii] == Places[j][Hosts[i].PlaceLinks[j]].members[k])
											{
												vaccflag = 0;
											}
											ii++;
										}
										if (vaccflag)
										{
											StateT[tn].ringvacclist[cnt] = Places[j][Hosts[i].PlaceLinks[j]].members[k];
											//StateT[tn].ringlist[cnt] = 2; //These people belong to ring 2
											if ((Hosts[Places[j][Hosts[i].PlaceLinks[j]].members[k]].vaccRing == 0) || (Hosts[Places[j][Hosts[i].PlaceLinks[j]].members[k]].vaccRing > currentRing))
											{
												Hosts[Places[j][Hosts[i].PlaceLinks[j]].members[k]].vaccRing = currentRing; //These people belong to ring 2
												Hosts[Places[j][Hosts[i].PlaceLinks[j]].members[k]].ringCase = ai; //associated with the ring around case ai
											}
											cnt++;
										}
									}
								}
							}
						}
					}
					//Now add household members of group contacts
					for (j = 0; j < P.PlaceTypeNum; j++)
					{
						if ((Hosts[ai].PlaceLinks[j] > 0) && (j != P.HospPlaceTypeNum))
						{
							for (k = 0; k < Places[j][Hosts[ai].PlaceLinks[j]].n; k++)
							{
								i1 = Households[Hosts[Places[j][Hosts[ai].PlaceLinks[j]].members[k]].hh].FirstPerson;
								i2 = i1 + Households[Hosts[Places[j][Hosts[ai].PlaceLinks[j]].members[k]].hh].nh;
								for (i = i1; i < i2; i++)
								{
									vaccflag = 1;
									if ((i != ai) && (!HOST_TO_BE_VACCED(i)))
									{
										ii = 0;
										while ((ii < cnt) && (vaccflag == 1))
										{
											if (StateT[tn].ringvacclist[ii] == i)
											{
												vaccflag = 0;
											}
											ii++;
										}
										if (vaccflag)
										{
											StateT[tn].ringvacclist[cnt] = i;
											if ((Hosts[i].vaccRing == 0) || (Hosts[i].vaccRing > currentRing))
											{
												Hosts[i].vaccRing = currentRing;
												Hosts[i].ringCase = ai;
											}
											//StateT[tn].ringlist[cnt] = currentRing; //These people belong to ring 2
											cnt++;
										}
									}
								}
							}
						}
					}
				}
				//if we want 3 rings, we also have to add household members of these group members
				if (P.NVaccRingsActive > 2)
				{
					currentRing = 3;
					//First add household members of group members of household contacts
					i1 = Households[Hosts[ai].hh].FirstPerson;
					i2 = i1 + Households[Hosts[ai].hh].nh;
					for (i = i1; i < i2; i++)
					{
						for (j = 0; j < P.PlaceTypeNum; j++)
						{
							if ((Hosts[ai].PlaceLinks[j] >= 0) && (j != P.HospPlaceTypeNum))
							{
								for (k = 0; k < Places[j][Hosts[i].PlaceLinks[j]].n; k++)
								{
									// we know the group member - now find people in their household
									k2 = Households[Hosts[Places[j][Hosts[i].PlaceLinks[j]].members[k]].hh].FirstPerson;
									k3 = k2 + Households[Hosts[Places[j][Hosts[i].PlaceLinks[j]].members[k]].hh].nh;
									for (l = k2; l < k3; l++)
									{
										vaccflag = 1;
										if ((l != ai) && (!HOST_TO_BE_VACCED(l)))// && (Hosts[l].vaccRing == 0))
										{
											ii = 0;
											while ((ii < cnt) && (vaccflag == 1))
											{
												if (StateT[tn].ringvacclist[ii] == l)
												{
													vaccflag = 0;
												}
												ii++;
											}
											if (vaccflag)
											{
												StateT[tn].ringvacclist[cnt] = l;
												//StateT[tn].ringlist[cnt] = 3; //These people belong to ring 3
												if ((Hosts[l].vaccRing == 0) || (Hosts[l].vaccRing > currentRing))
												{
													Hosts[l].vaccRing = currentRing; //These people belong to ring 3
													Hosts[l].ringCase = ai; //Associated with the ring around case ai
												}
												cnt++;
											}
										}
									}
								}
							}
						}
					}
					//Now add group members of household contacts of group members
					for (j = 0; j < P.PlaceTypeNum; j++)
					{
						if ((Hosts[ai].PlaceLinks[j] >= 0) && (j != P.HospPlaceTypeNum))
						{
							for (k = 0; k < Places[j][Hosts[ai].PlaceLinks[j]].n; k++)
							{
								i1 = Households[Hosts[Places[j][Hosts[ai].PlaceLinks[j]].members[k]].hh].FirstPerson;
								i2 = i1 + Households[Hosts[Places[j][Hosts[ai].PlaceLinks[j]].members[k]].hh].nh;
								for (i = i1; i < i2; i++)
								{
									for (k2 = 0; k2 < P.PlaceTypeNum; k2++)
									{
										if (Hosts[i].PlaceLinks[k2] >= 0)
										{
											for (k3 = 0; k3 < Places[k2][Hosts[i].PlaceLinks[k2]].n; k3++)
											{
												vaccflag = 1;
												if ((Places[k2][Hosts[i].PlaceLinks[k2]].members[k3] != ai) && (!HOST_TO_BE_VACCED(Places[k2][Hosts[i].PlaceLinks[k2]].members[k3])))// && (Hosts[Places[k2][Hosts[i].PlaceLinks[k2]].members[k3]].vaccRing == 0))
												{
													ii = 0;
													while ((ii < cnt) && (vaccflag == 1))
													{
														if (StateT[tn].ringvacclist[ii] == Places[k2][Hosts[i].PlaceLinks[k2]].members[k3])
														{
															vaccflag = 0;
														}
														ii++;
													}
													if (vaccflag)
													{
														StateT[tn].ringvacclist[cnt] = Places[k2][Hosts[i].PlaceLinks[k2]].members[k3];
														//StateT[tn].ringlist[cnt] = 3; //These people belong to ring 3
														if ((Hosts[Places[k2][Hosts[i].PlaceLinks[k2]].members[k3]].vaccRing == 0) || (Hosts[Places[k2][Hosts[i].PlaceLinks[k2]].members[k3]].vaccRing > currentRing))
														{
															Hosts[Places[k2][Hosts[i].PlaceLinks[k2]].members[k3]].vaccRing = currentRing;
															Hosts[Places[k2][Hosts[i].PlaceLinks[k2]].members[k3]].ringCase = ai;
														}
														cnt++;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}


				StateT[tn].vaccdosering_dist[VACCDOSE_GROUP((int)cnt)]++;
				//Finally do the vaccination step
				for (i = 0; i < cnt; i++)
				{
					//
					if (HOST_AGE_YEAR(StateT[tn].ringvacclist[i]) >= P.MinVaccAge)
					{
						//fprintf(stderr, "Max vaccine dose per day: %i\n", P.VaccDosePerDay);
						if ((Hosts[i].vacc_accept < P.PropRingVacc) && !((Hosts[StateT[tn].ringvacclist[i]].inf >= 5))) //changed from P.PropRingVacc to adjPropToVacc
						{

							StateT[tn].vacc_queue[StateT[tn].nvacc_queue] = StateT[tn].ringvacclist[i];
							StateT[tn].ring_queue[StateT[tn].nvacc_queue] = 1; //set this to one if added to the queue via ring vaccination - it will be zero for geo
							StateT[tn].nvacc_queue++;
						}
					}
				}
				//fprintf(stderr, "cnt=%i\t, State.cumV=%lg\t, Prop to vacc:%lg\n", cnt, State.cumV, P.PropRingVacc);
			}
		}
	}

}

void DoCase(int ai, double t, unsigned short int ts, int tn)
{
	int j, k, f, j1, j2, m, q;
	person* a;
	int age;
	//int *RingVaccList,*RingVaccNonHouseholdList; //added this for my ring vaccination code - ggilani 15/02/17
	//int *RingVaccRingList; //added this to keep track of which ring each contact belongs too - ggilani 29/05/2019
	int currentRing; //to keep track of the current ring
	int nVacc, casePlaceType, i, i1, i2, cnt, k2, k3, l, nAlreadyVacc, ii, vaccflag;
	double propVacc, adjPropToVacc;

	casePlaceType = 0; //added this to initialise casePlaceType memory - ggilani 21/10/19
	currentRing = 0; //to initialise - ggilani 29/10/19

	age = HOST_AGE_GROUP(ai);
	if (age >= NUM_AGE_GROUPS) age = NUM_AGE_GROUPS - 1;
	a = Hosts + ai;
	if (a->inf == -1)
	{
		a->inf = -2;
		if (HOST_ABSENT(ai))
		{
			if (a->absent_stop_time < ts + P.usCaseAbsenteeismDelay + P.usCaseAbsenteeismDuration)
				a->absent_stop_time = ts + P.usCaseAbsenteeismDelay + P.usCaseAbsenteeismDuration;
		}
		else
		{
			a->absent_start_time = USHRT_MAX - 1;
			if (P.DoPlaces)
			{
				for (j = 0; j < P.PlaceTypeNum; j++)
				{
					if ((a->PlaceLinks[j] >= 0) && (j != HOTEL_PLACE_TYPE) && (!HOST_ABSENT(ai)) && (P.SymptPlaceTypeWithdrawalProp[j] > 0))
					{
						if ((P.SymptPlaceTypeWithdrawalProp[j] == 1) || (ranf_mt(tn) < P.SymptPlaceTypeWithdrawalProp[j]))
						{
							a->absent_start_time = ts + P.usCaseAbsenteeismDelay;
							a->absent_stop_time = ts + P.usCaseAbsenteeismDelay + P.usCaseAbsenteeismDuration;
#ifdef ABSENTEEISM_PLACE_CLOSURE
							if (t >= P.PlaceCloseTimeStart)
								for (j = 0; j < P.PlaceTypeNum; j++)
									if ((j != HOTEL_PLACE_TYPE) && (a->PlaceLinks[j] >= 0))
									{
										if ((!P.RestrictTreatToTarget) || (Places[j][a->PlaceLinks[j]].country == P.TargetCountry))
											DoPlaceClose(j, a->PlaceLinks[j], ts, tn, 0);
									}
#endif
							if ((!HOST_QUARANTINED(ai)) && (Hosts[ai].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] >= 0) && (HOST_AGE_YEAR(ai) >= P.CaseAbsentChildAgeCutoff))
								StateT[tn].cumAC++;
							/* This calculates adult absenteeism from work due to care of sick children. Note, children not at school not counted (really this should
							be fixed in population setup by having adult at home all the time for such kids. */
							if ((P.DoHouseholds) && (HOST_AGE_YEAR(ai) < P.CaseAbsentChildAgeCutoff))
							{
								if (!HOST_QUARANTINED(ai)) StateT[tn].cumACS++;
								if ((P.CaseAbsentChildPropAdultCarers == 1) || (ranf_mt(tn) < P.CaseAbsentChildPropAdultCarers))
								{
									j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
									f = 0;
									for (j = j1; (j < j2) && (!f); j++)
										f = ((abs(Hosts[j].inf) != 5) && (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && (Hosts[j].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] < 0));
									if (!f)
									{
										for (j = j1; (j < j2) && (!f); j++)
											f = ((abs(Hosts[j].inf) != 5) && (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && ((HOST_ABSENT(j)) || (HOST_QUARANTINED(j))));
										if (!f)
										{
											f = 0; k = -1;
											for (j = j1; (j < j2) & (!f); j++)
												if ((HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && (abs(Hosts[j].inf) != 5)) { k = j; f = 1; }
											if (f)
											{
												Hosts[k].absent_start_time = ts + P.usCaseIsolationDelay;
												Hosts[k].absent_stop_time = ts + P.usCaseIsolationDelay + P.usCaseIsolationDuration;
												StateT[tn].cumAA++;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		//if((P.ControlPropCasesId==1)||(ranf_mt(tn)<P.ControlPropCasesId))
		//	{
		//	StateT[tn].cumDC++;
		//	DoDetectedCase(ai,t,ts,tn);
		//	}
		// 
		// 

		if (HOST_TREATED(ai)) Cells[Hosts[ai].pcell].cumTC++;
		StateT[tn].cumC++;
		StateT[tn].cumCa[age]++;
		StateT[tn].cumC_country[Mcells[Hosts[ai].mcell].country]++; //add to cumulative count of cases in that country: ggilani - 12/11/14
		StateT[tn].cumC_resist[a->resist]++;
		StateT[tn].cumC_keyworker[a->keyworker]++;

		//added some case detection code here: ggilani - 03/02/15, updated and moved: 11/05/22
		//Add some case detection code here which determines whether someone will be detected once they have become a case, and determines the time point at which they will be detected
		//based on an input distribution function. 
		if (P.OutbreakDetected)
		{
			if (Hosts[ai].rep_rate < P.ControlPropCasesId || (Hosts[ai].contactTraced) || (Hosts[ai].hcw) || (Hosts[ai].flw)) //changed this to proportion seeking care
			{
				Hosts[ai].detected = 1;
			}
		}
		else
		{
			if (Hosts[ai].rep_rate < (P.ControlPropCasesId * P.PropHospSeek))
			{
				Hosts[ai].detected = 1;
			}
		}

		//if (P.DoCaseDetection)
		//{
		//	if ((P.DoCaseDetectionAdunit) && (P.DoAdUnits))
		//	{
		//		if ((Hosts[ai].rep_rate < AdUnits[Mcells[a->mcell].adunit].caseDetectRate) || (Hosts[ai].contactTraced > 0)) //add something to ensure that hosts who are being currently contact traced are always detected, changed ranf_mt(tn) to Hosts[ai].rep_rate 13/06/19
		//		{
		//			Hosts[ai].detected = 1;
		//		}
		//	}
		//	else
		//	{
		//		if (Hosts[ai].rep_rate < P.CaseDetectionRate || (Hosts[ai].contactTraced > 0))
		//		{
		//			Hosts[ai].detected = 1;
		//		}
		//	}

		//assign a time to detection, similar to code in DoInfect, when latent period is determined, but time to detection must be before time to recovery/death and then they
		//leave the system and there's nothing we can do to train them and update the number of detected cases.
		if (Hosts[ai].detected == 1)
		{
			if (P.DoDetectDelay)
			{
				//--> this bit is for distribution
				//
				//do
				//{
				//i = (int)floor((q = ranf_mt(tn) * CDF_RES));
				//q -= ((double)i);
				//Hosts[ai].detect_time = Hosts[ai].latent_time + (unsigned short int) floor(0.5 + (t - P.DetectTime * log(q * P.detect_icdf[i + 1] + (1.0 - q) * P.detect_icdf[i])) * P.TimeStepsPerDay);
				//while(Hosts[ai].detect_time<=Hosts[ai].recoverytime);
				if (Hosts[ai].contactTraced == 0)
				{
					if ((unsigned short int) (Hosts[ai].latent_time + (P.DetectTime * P.TimeStepsPerDay)) < Hosts[ai].recovery_time)
					{
						Hosts[ai].detect_time = (unsigned short int) (Hosts[ai].latent_time + (P.DetectTime * P.TimeStepsPerDay) + (P.LatentToSymptDelay / P.TimeStep)); //currently using a fixed delay
					}
					else
					{
						Hosts[ai].detect_time = Hosts[ai].recovery_time; //if the delay to 
					}
				}
				else
				{

					Hosts[ai].detect_time = Hosts[ai].latent_time + ((int)(P.LatentToSymptDelay / P.TimeStep)); //if contact traced, detected immediately, set detect_time to be the same time
				}
			}
			else
			{
				Hosts[ai].detect_time = Hosts[ai].latent_time + ((int)(P.LatentToSymptDelay / P.TimeStep)); //if detected immediately, set detect_time to be the same time
			}

			if ((P.DoHospitalisation))//&&(t>=P.ETUTimeStart))
			{
				if ((P.OutbreakDetected) && (ranf_mt(tn) < P.PropHospSeek))
				{
					if (Hosts[ai].contactTraced == 0)
					{
						a->hospital_time = a->latent_time + (unsigned short int) floor(0.5 + (P.HospitalisationTime * P.TimeStepsPerDay));
					}
					else
					{
						a->hospital_time = a->latent_time + (unsigned short int) floor(0.5 + (P.HospitalisationTime_contactTrace * P.TimeStepsPerDay)); //different hospitalisation time for contact traced case: ggilani 05/07/2017
					}
				}
				else
				{
					a->hospital_time = a->latent_time + (unsigned short int) floor(0.5 + (P.HospitalisationTime * P.TimeStepsPerDay));
				}
			}
		}
		//}

		if (P.DoAdUnits) StateT[tn].cumC_adunit[Mcells[a->mcell].adunit]++;
	}
}

void DoFalseCase(int ai, double t, unsigned short int ts, int tn)
{
	int j, k, f, j1, j2;
	person* a;

	/* Arguably adult absenteeism to take care of sick kids could be included here, but then output absenteeism would not be 'excess' absenteeism */
	a = Hosts + ai;
	if ((P.ControlPropCasesId == 1) || (ranf_mt(tn) < P.ControlPropCasesId))
	{
		if ((!P.DoEarlyCaseDiagnosis) || (State.cumDC >= P.PreControlClusterIdCaseThreshold)) StateT[tn].cumDC++;
		DoDetectedCase(ai, t, ts, tn);
	}
	StateT[tn].cumFC++;
}

void DoRecover(int ai, int run, int tn)
{
	int i, j, x, y;
	person* a;


	a = Hosts + ai;
	if (abs(a->inf) == 2)
	{
		i = a->listpos;
		Cells[a->pcell].I--;
		if (P.DoSIS)
		{
			if (Cells[a->pcell].I > 0)
			{
				Cells[a->pcell].susceptible[i] = Cells[a->pcell].infected[0];
				Hosts[Cells[a->pcell].susceptible[i]].listpos = i;
				if (Cells[a->pcell].L > 0)
				{
					Cells[a->pcell].latent[Cells[a->pcell].L] = Cells[a->pcell].latent[0];
					Hosts[Cells[a->pcell].latent[Cells[a->pcell].L]].listpos = Cells[a->pcell].S + Cells[a->pcell].L;
				}
			}
			else if (Cells[a->pcell].L > 0)
			{
				Cells[a->pcell].susceptible[i] = Cells[a->pcell].latent[0];
				Hosts[Cells[a->pcell].susceptible[i]].listpos = i;
			}
			Cells[a->pcell].susceptible[Cells[a->pcell].S] = ai;
			a->listpos = Cells[a->pcell].S;
			Cells[a->pcell].S++;
			Cells[a->pcell].latent++;
			Cells[a->pcell].infected++;
			a->susc *= P.SuscReductionFactorPerInfection;
			a->inf = 0;
			a->infector = -1;
			if (P.OutputBitmap)
			{
				if ((P.OutputBitmapDetected == 0) || ((P.OutputBitmapDetected == 1) && (Hosts[ai].detected == 1)))
				{

					x = ((int)(Households[a->hh].loc_x * P.scalex)) - P.bminx;
					y = ((int)(Households[a->hh].loc_y * P.scaley)) - P.bminy;
					if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
					{
						j = y * bmh->width + x;
						if ((j < bmh->imagesize) && (j >= 0))
						{
#pragma omp atomic
							bmi2[j]--;
						}
					}
				}
			}

		}
		else
		{
			Cells[a->pcell].R++;
			if ((P.DoAdUnits) && (a->detected))
			{
				StateT[tn].cumDR_adunit[Mcells[a->mcell].adunit]++;
			}
			if (i < Cells[a->pcell].S + Cells[a->pcell].L + Cells[a->pcell].I)
			{
				Cells[a->pcell].susceptible[i] = Cells[a->pcell].infected[Cells[a->pcell].I];
				Hosts[Cells[a->pcell].susceptible[i]].listpos = i;
			}
			a->inf = 3 * a->inf / abs(a->inf);
			a->listpos = Cells[a->pcell].S + Cells[a->pcell].L + Cells[a->pcell].I;
			Cells[a->pcell].susceptible[a->listpos] = ai;
			if (P.OutputBitmap)
			{
				if ((P.OutputBitmapDetected == 0) || ((P.OutputBitmapDetected == 1) && (Hosts[ai].detected == 1)))
				{
					x = ((int)(Households[a->hh].loc_x * P.scalex)) - P.bminx;
					y = ((int)(Households[a->hh].loc_y * P.scaley)) - P.bminy;
					if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
					{
						j = y * bmh->width + x;
						if ((j < bmh->imagesize) && (j >= 0))
						{
#pragma omp atomic
							bmi3[j]++;
#pragma omp atomic
							bmi2[j]--;
						}
					}
				}
			}
		}
	}
	//added this to record event if flag is set to 1 and if host isn't initial seed, i.e. if Hosts[ai].infector>=0: ggilani - 10/10/2014
	//if(P.DoRecordInfEvents)
	//{
	//	if(*nEvents<P.MaxInfEvents)
	//	{
	//		RecordEvent(((double)a->recovery_time)*P.TimeStep,ai,run,3); //added int as argument to RecordEvent to record run number: ggilani - 15/10/14
	//	}
	//}
}

void DoDeath(int ai, int tn, int run)
{
	int j, x, y;
	person* a;

	a = Hosts + ai;
	if ((abs(a->inf) == 2) || (abs(a->inf) == 6)) //added infection status of 6 as well, as this also needs to deal with 'dead' infectious individuals: ggilani 25/10/14
	{
		a->inf = 5 * a->inf / abs(a->inf);
		Cells[a->pcell].D++;
		Cells[a->pcell].I--;
		if (Cells[a->pcell].I > 0)
		{
			Cells[a->pcell].susceptible[a->listpos] = Cells[a->pcell].infected[Cells[a->pcell].I];
			Hosts[Cells[a->pcell].susceptible[a->listpos]].listpos = a->listpos;
		}
		a->listpos = Cells[a->pcell].S + Cells[a->pcell].L + Cells[a->pcell].I;
		Cells[a->pcell].susceptible[a->listpos] = ai;
		/*		a->listpos=-1; */
		if (abs(a->inf) != 6) //if if is someone who is infectious after death, we'll do this bit of accounting in IncubRecoverySweep to make sure their death is recorded on the right day, even though they will stay in the infectious list
		{
			StateT[tn].cumD++;
			StateT[tn].cumDa[HOST_AGE_GROUP(ai)]++;
			if (P.DoAdUnits)
			{
				StateT[tn].cumD_adunit[Mcells[a->mcell].adunit]++;
				if (a->detected) StateT[tn].cumDD_adunit[Mcells[a->mcell].adunit]++;
			}
			StateT[tn].cumD_keyworker[a->keyworker]++;
		}
		if (P.OutputBitmap)
		{
			if ((P.OutputBitmapDetected == 0) || ((P.OutputBitmapDetected == 1) && (Hosts[ai].detected == 1)))
			{
				x = ((int)(Households[a->hh].loc_x * P.scalex)) - P.bminx;
				y = ((int)(Households[a->hh].loc_y * P.scaley)) - P.bminy;
				if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
				{
					j = y * bmh->width + x;
					if ((j < bmh->imagesize) && (j >= 0))
					{
#pragma omp atomic
						bmi3[j]++;
#pragma omp atomic
						bmi2[j]--;
					}
				}
			}
		}
	}
	////added this to record event if flag is set to 1 and if host isn't initial seed, i.e. if Hosts[ai].infector>=0: ggilani - 10/10/2014
	//if(P.DoRecordInfEvents)
	//{
	//	if(*nEvents<P.MaxInfEvents)
	//	{
	//		RecordEvent(((double)a->recovery_time)*P.TimeStep,ai,run,2); //added int as argument to RecordEvent to record run number: ggilani - 15/10/14
	//	}
	//}
}

void DoTreatCase(int ai, unsigned short int ts, int tn)
{
	int j;
	double x, y;

	if (State.cumT < P.TreatMaxCourses)
	{
		if ((!HOST_TREATED(ai)) && (Hosts[ai].resist < (MAX_NUM_RESIST_TYPES - 1)) && (ranf_mt(tn) < P.EvolResistTreatMutationRate)) Hosts[ai].resist++;
#ifdef NO_TREAT_PROPH_CASES
		if (!HOST_TO_BE_TREATED(ai))
#endif
		{
			Hosts[ai].treat_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.TreatDelayMean));
			Hosts[ai].treat_stop_time = ts + ((unsigned short int) (P.TimeStepsPerDay * (P.TreatDelayMean + P.TreatCaseCourseLength)));
			StateT[tn].cumT++;
			if ((abs(Hosts[ai].inf) > 0) && (Hosts[ai].inf != 5)) Cells[Hosts[ai].pcell].cumTC++;
			StateT[tn].cumT_resist[Hosts[ai].resist]++;
			StateT[tn].cumT_keyworker[Hosts[ai].keyworker]++;
			if ((++Hosts[ai].num_treats) < 2) StateT[tn].cumUT++;
			Cells[Hosts[ai].pcell].tot_treat++;
			if (P.DoAdUnits) StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit]++;
			if (P.OutputBitmap)
			{
				x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
				y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
				if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
				{
					j = y * bmh->width + x;
					if ((j < bmh->imagesize) && (j >= 0))
					{
#pragma omp atomic
						bmi4[j]++;
					}
				}
			}
		}
	}
}

void DoProph(int ai, unsigned short int ts, int tn)
{
	int j;
	double x, y;

	if (State.cumT < P.TreatMaxCourses)
	{
		if ((!HOST_TREATED(ai)) && (Hosts[ai].resist < (MAX_NUM_RESIST_TYPES - 1)) && ((abs(Hosts[ai].inf) == 1) || (abs(Hosts[ai].inf) == 2)) && (ranf_mt(tn) < P.EvolResistProphMutationRate)) Hosts[ai].resist++;
		Hosts[ai].treat_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.TreatDelayMean));
		Hosts[ai].treat_stop_time = ts + ((unsigned short int) (P.TimeStepsPerDay * (P.TreatDelayMean + P.TreatProphCourseLength)));
		StateT[tn].cumT++;
		//		StateT[tn].cumT_resist[Hosts[ai].resist]++;
		StateT[tn].cumT_keyworker[Hosts[ai].keyworker]++;
		if ((++Hosts[ai].num_treats) < 2) StateT[tn].cumUT++;
		if (P.DoAdUnits)	StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit]++;
#pragma omp critical (tot_treat)
		Cells[Hosts[ai].pcell].tot_treat++;
		if (P.OutputBitmap)
		{
			x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				j = y * bmh->width + x;
				if ((j < bmh->imagesize) && (j >= 0))
				{
#pragma omp atomic
					bmi4[j]++;
				}
			}
		}
	}
}

void DoPrivateTreatCase(int ai, unsigned short int ts, int tn)
{
	int j;
	double x, y;

	if ((!HOST_TREATED(ai)) && (Hosts[ai].resist < (MAX_NUM_RESIST_TYPES - 1)) && (ranf_mt(tn) < P.EvolResistTreatMutationRate)) Hosts[ai].resist++;
#ifdef NO_TREAT_PROPH_CASES
	if (!HOST_TO_BE_TREATED(ai))
#endif
	{
		Households[Hosts[ai].hh].stockpile = -1;
		Hosts[ai].treat_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.PrivateTreatDelayMean));
		Hosts[ai].treat_stop_time = ts + ((unsigned short int) (P.TimeStepsPerDay * (P.PrivateTreatDelayMean + P.TreatCaseCourseLength)));
		StateT[tn].cumTP++;
		if ((abs(Hosts[ai].inf) > 0) && (Hosts[ai].inf != 5))Cells[Hosts[ai].pcell].cumTC++;
		StateT[tn].cumT_resist[Hosts[ai].resist]++;
		StateT[tn].cumT_keyworker[Hosts[ai].keyworker]++;
		if ((++Hosts[ai].num_treats) < 2) StateT[tn].cumUT++;
		Cells[Hosts[ai].pcell].tot_treat++;
		if (P.DoAdUnits) StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit]++;
		if (P.OutputBitmap)
		{
			x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				j = y * bmh->width + x;
				if ((j < bmh->imagesize) && (j >= 0))
				{
#pragma omp atomic
					bmi4[j]++;
				}
			}
		}
	}
}

void DoPrivateProph(int ai, unsigned short int ts, int tn)
{
	int j;
	double x, y;

	if ((!HOST_TREATED(ai)) && (Hosts[ai].resist < (MAX_NUM_RESIST_TYPES - 1)) && ((abs(Hosts[ai].inf) == 1) || (abs(Hosts[ai].inf) == 2)) && (ranf_mt(tn) < P.EvolResistProphMutationRate)) Hosts[ai].resist++;
	Hosts[ai].treat_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.PrivateTreatDelayMean));
	Hosts[ai].treat_stop_time = ts + ((unsigned short int) (P.TimeStepsPerDay * (P.PrivateTreatDelayMean + P.TreatProphCourseLength)));
	StateT[tn].cumTP++;
	//StateT[tn].cumT_resist[Hosts[ai].resist]++;
	StateT[tn].cumT_keyworker[Hosts[ai].keyworker]++;
	if ((++Hosts[ai].num_treats) < 2) StateT[tn].cumUT++;
	if (P.DoAdUnits)	StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit]++;
#pragma omp critical (tot_treat)
	Cells[Hosts[ai].pcell].tot_treat++;
	if (P.OutputBitmap)
	{
		x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
		y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
		if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
		{
			j = y * bmh->width + x;
			if ((j < bmh->imagesize) && (j >= 0))
			{
#pragma omp atomic
				bmi4[j]++;
			}
		}
	}
}

void DoProphNoDelay(int ai, unsigned short int ts, int tn, int nc)
{
	int j;
	double x, y;

	if (State.cumT < P.TreatMaxCourses)
	{
		if ((!HOST_TREATED(ai)) && (Hosts[ai].resist < (MAX_NUM_RESIST_TYPES - 1)) && ((abs(Hosts[ai].inf) == 1) || (abs(Hosts[ai].inf) == 2)) && (ranf_mt(tn) < P.EvolResistProphMutationRate)) Hosts[ai].resist++;
		Hosts[ai].treat_start_time = ts;
		Hosts[ai].treat_stop_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.TreatProphCourseLength * nc));
		StateT[tn].cumT += nc;
		//StateT[tn].cumT_resist[Hosts[ai].resist]+=nc;
		StateT[tn].cumT_keyworker[Hosts[ai].keyworker] += nc;
		if ((++Hosts[ai].num_treats) < 2) StateT[tn].cumUT++;
		if (P.DoAdUnits) StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit] += nc;
#pragma omp critical (tot_treat)
		Cells[Hosts[ai].pcell].tot_treat++;
		if (P.OutputBitmap)
		{
			x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				j = y * bmh->width + x;
				if ((j < bmh->imagesize) && (j >= 0))
				{
#pragma omp atomic
					bmi4[j]++;
				}
			}
		}
	}
}

void DoPlaceClose(int i, int j, unsigned short int ts, int tn, int DoAnyway)
{
	int k, ai, j1, j2, l, f, m, f2;
	unsigned short trig;
#ifdef ABSENTEEISM_PLACE_CLOSURE
	unsigned short int t_old, t_new;
#endif


	f2 = 0;
	/*	if((j<0)||(j>=P.Nplace[i]))
			fprintf(stderr,"** %i %i *\n",i,j);
		else
	*/
#ifdef ABSENTEEISM_PLACE_CLOSURE
	t_new = ts / P.TimeStepsPerDay;
#endif
	trig = 0;
#pragma omp critical (closeplace)
	{
		if (!PLACE_TO_BE_CLOSED(i, j))
		{
			if ((!DoAnyway) && (Places[i][j].control_trig < USHRT_MAX - 2))
			{
#ifdef ABSENTEEISM_PLACE_CLOSURE
				t_old = Places[i][j].AbsentLastUpdateTime;
				if (t_new >= t_old + MAX_ABSENT_TIME)
					for (l = 0; l < MAX_ABSENT_TIME; l++) Places[i][j].Absent[l] = 0;
				else
					for (l = t_old; l < t_new; l++) Places[i][j].Absent[l % MAX_ABSENT_TIME] = 0;
				for (l = t_new; l < t_new + P.usCaseAbsenteeismDuration / P.TimeStepsPerDay; l++) Places[i][j].Absent[l % MAX_ABSENT_TIME]++;
				trig = Places[i][j].Absent[t_new % MAX_ABSENT_TIME];
				Places[i][j].AbsentLastUpdateTime = t_new;
				if ((P.PlaceCloseByAdminUnit) && (P.PlaceCloseAdunitPlaceTypes[i] > 0)
					&& (((double)trig) / ((double)Places[i][j].n) > P.PlaceCloseCasePropThresh))
				{
					//fprintf(stderr,"** %i %i %i %i %lg ## ",i,j,(int) Places[i][j].control_trig, (int) Places[i][j].n,P.PlaceCloseCasePropThresh);
					k = Mcells[Places[i][j].mcell].adunit;
					if (AdUnits[k].place_close_trig < USHRT_MAX - 1) AdUnits[k].place_close_trig++;
				}
#else
				trig = (++Places[i][j].control_trig);
				if ((P.PlaceCloseByAdminUnit) && (P.PlaceCloseAdunitPlaceTypes[i] > 0)
					&& (((double)Places[i][j].control_trig) / ((double)Places[i][j].n) > P.PlaceCloseCasePropThresh))
				{
					//fprintf(stderr,"** %i %i %i %i %lg ## ",i,j,(int) Places[i][j].control_trig, (int) Places[i][j].n,P.PlaceCloseCasePropThresh);
					k = Mcells[Places[i][j].mcell].adunit;
					if (AdUnits[k].place_close_trig < USHRT_MAX - 1) AdUnits[k].place_close_trig++;
				}
#endif
			}
			if (Places[i][j].control_trig < USHRT_MAX - 1)
			{
				if (P.PlaceCloseFracIncTrig > 0)
					k = (((double)trig) / ((double)Places[i][j].n) > P.PlaceCloseFracIncTrig);
				else
					k = (((int)trig) >= P.PlaceCloseIncTrig);
				if (((!P.PlaceCloseByAdminUnit) && (k)) || (DoAnyway))
				{
					if (P.DoPlaceCloseOnceOnly)
						Places[i][j].control_trig = USHRT_MAX - 1;  // Places only close once
					else
						Places[i][j].control_trig = 0;
					if ((P.PlaceCloseEffect[i] == 0) || (ranf_mt(tn) >= P.PlaceCloseEffect[i]))
					{
						Places[i][j].close_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.PlaceCloseDelayMean));
						Places[i][j].close_end_time = ts + ((unsigned short int) (P.TimeStepsPerDay * (P.PlaceCloseDelayMean + P.PlaceCloseDuration)));
						f2 = 1;
					}
					else
						Places[i][j].close_start_time = Places[i][j].close_end_time = ts + ((unsigned short int) (P.TimeStepsPerDay * (P.PlaceCloseDelayMean + P.PlaceCloseDuration)));
				}
			}
		}
	}
	if (f2)
	{
		if (P.DoRealSymptWithdrawal)
			for (k = 0; k < Places[i][j].n; k++)
			{
				ai = Places[i][j].members[k];
				Hosts[ai].absent_start_time = Places[i][j].close_start_time;
				Hosts[ai].absent_stop_time = Places[i][j].close_end_time;
				if ((HOST_AGE_YEAR(ai) >= P.CaseAbsentChildAgeCutoff) && (Hosts[ai].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] >= 0)) StateT[tn].cumAPC++;
				if ((HOST_AGE_YEAR(ai) < P.CaseAbsentChildAgeCutoff) && (!HOST_ABSENT(ai)) && (!HOST_QUARANTINED(ai)))
				{
					StateT[tn].cumAPCS++;
					if ((P.CaseAbsentChildPropAdultCarers == 1) || (ranf_mt(tn) < P.CaseAbsentChildPropAdultCarers))
					{
						j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
						if ((j1 < 0) || (j2 >= P.N)) fprintf(stderr, "++ %i %i %i (%i %i %i)##  ", ai, j1, j2, i, j, k);
						f = 0;
						for (l = j1; (l < j2) && (!f); l++)
							f = ((abs(Hosts[l].inf) != 5) && (HOST_AGE_YEAR(l) >= P.CaseAbsentChildAgeCutoff) &&
								((Hosts[l].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] < 0) || (HOST_ABSENT(l)) || (HOST_QUARANTINED(l))));
						if (!f)
						{
							for (l = j1; (l < j2) && (!f); l++)
								if ((HOST_AGE_YEAR(l) >= P.CaseAbsentChildAgeCutoff) && (abs(Hosts[l].inf) != 5)) { m = l; f = 1; }
							if (f)
							{
								Hosts[m].absent_start_time = Places[i][j].close_start_time;
								Hosts[m].absent_stop_time = Places[i][j].close_end_time;
								StateT[tn].cumAPA++;
							}
						}
					}
				}
			}
	}
}


int DoVacc(int ai, int ts, int ringflag)
{
	int j;
	double x, y;
	int age;

	age = HOST_AGE_GROUP(ai);

	if ((!P.DoDistributionVaccination) && (State.cumV >= P.VaccMaxCourses))
		return 2;
	else if ((HOST_TO_BE_VACCED(ai) && (Hosts[ai].vacc_start_time > 0)) || (Hosts[ai].inf < -1) || (Hosts[ai].inf >= 5))
		return 1;
	else
	{
		if (Hosts[ai].keyworker && Hosts[ai].vacc_start_time <= 0)
		{
			Hosts[ai].revacc = 1;
		}
		Hosts[ai].vacc_start_time = ts + ((int)(P.TimeStepsPerDay * P.VaccDelayMean));
		//if (flag == 1)
		//{
		//	//This person was in a third ring and it will take longer for efficacy - mark them as such
		//	Hosts[ai].thirdVaccRing = 1;
		//}
		if (ringflag == 1) //for ring vaccination
		{
#pragma omp critical (state_cumV)
			State.cumV++;
#pragma omp critical (state_cumVa)
			State.cumVa[age]++; //added vaccination by age: ggilani 23/02/22
#pragma omp critical (state_cumVad)
			if (P.DoAdUnits)
			{
				StateT[0].cumV_adunit[Mcells[Hosts[ai].mcell].adunit]++;
			}
#pragma omp critical (state_cumV_daily)
			if (P.VaccDosePerDay >= 0)
			{
				State.cumV_daily++;
			}
		}
		else if (ringflag == 0) //for geographic vaccination
		{
#pragma omp critical (state_cumVG)
			State.cumVG++;
			//#pragma omp critical (state_cumVa)
			//			State.cumVGa[age]++; //added vaccination by age: ggilani 23/02/22
#pragma omp critical (state_cumVGad)
			if (P.DoAdUnits)
			{
				StateT[0].cumVG_adunit[Mcells[Hosts[ai].mcell].adunit]++;
			}
#pragma omp critical (state_cumVG_daily)
			if (P.VaccGeoDosePerDay >= 0)
			{
				State.cumVG_daily++;
			}
		}

#pragma omp critical (tot_vacc)
		Cells[Hosts[ai].pcell].tot_vacc++;
		Mcells[Hosts[ai].mcell].popvacc++;
		if (P.OutputBitmap)
		{
			x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				j = y * bmh->width + x;
				if ((j < bmh->imagesize) && (j >= 0))
				{
#pragma omp atomic
					bmi4[j]++;
				}
			}
		}
	}
	return 0;
}

void DoVaccNoDelay(int ai, int ts)
{
	int j;
	double x, y;

	if ((State.cumV < P.VaccMaxCourses) && (!HOST_TO_BE_VACCED(ai)) && (Hosts[ai].inf >= -1) && (Hosts[ai].inf < 5)) //changing cumulative vaccine doses to VG - to separate ring from geographic counts
	{
		Hosts[ai].vacc_start_time = ts;
#pragma omp critical (state_cumV) //changed to VG
		State.cumV++; //changed to VG
#pragma omp critical (state_cumV_daily)
		if (P.VaccDosePerDay >= 0)
		{
			State.cumV_daily++;
		}
#pragma omp critical (state_cumVad)
		if (P.DoAdUnits)
		{
			StateT[0].cumV_adunit[Mcells[Hosts[ai].mcell].adunit]++;
		}
#pragma omp critical (tot_vacc)
		Cells[Hosts[ai].pcell].tot_vacc++;
		if (P.OutputBitmap)
		{
			x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				j = y * bmh->width + x;
				if ((j < bmh->imagesize) && (j >= 0))
				{
#pragma omp atomic
					bmi4[j]++;
				}
			}
		}
	}
}

