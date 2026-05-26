/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "model_constants.h"
#include "model_structs.h"
#include "model_fns.h"
#include "model_globals_ext.h"
#include "binio.h"



void RunModel(int run) //added run number as parameter
{
	int i, j, k, l, fs, fs2, ni, i2, nu, in;
	double ir, t, cI, lcI, t2, vaccUsed;
	unsigned short int ts, ts_prev;
	int continueEvents = 1;

	lcI = 1;
	if (P.DoLoadSnapshot)
	{
		P.ts_age = (int)(P.SnapshotLoadTime * P.TimeStepsPerDay);
		t = ((double)P.ts_age) * P.TimeStep;
	}
	else
	{
		t = 0;
		P.ts_age = 0;
	}
	fs = 1;
	fs2 = 0;
	nu = 0;

	for (ns = 1; ((ns < P.NumSamples) && (!InterruptRun)); ns++) //&&(continueEvents) <-removed this
	{
		if (continueEvents)
		{
			RecordSample(t, ns - 1);
			//update hospitalisation parameters at the beginning of every time step? ggilani - 11/03/2017
			if ((P.DoHospitalisation) && (t >= P.ETUTimeStart))
			{
				UpdateHospitals(t);
			}
			//update contact tracing parameters at the beginning of every time step? ggilani - 11/03/2017
			if ((P.DoContactTracing) && (t >= P.ContactTracingTimeStart))
			{
				UpdateContactTracing(t);
			}
			//update safe burial parameters parameters at the beginning of every time step? ggilani - 11/03/2017
			if ((P.DoFuneralTransmission) && (t >= P.FuneralControlTimeStart))
			{
				UpdateSDB(t);
			}
			//update vaccination parameters at the beginning of every time step
			if (((P.DoRingVaccination) && (t > P.VaccTimeStart)) || ((P.DoGeoVaccination) && (t > P.VaccTimeStart)))
			{
				UpdateVaccination(t, ns - 1);
			}
			//update vaccination parameters at the beginning of every time step
			if (P.DoUpdateCaseDetection)//&&(t>=P.TimeToUpdateCaseDetection))
			{
				UpdateCaseDetection(t);
			}
			if (P.VaccDosePerDay >= 0) //if constrained by total number of vaccine doses per day, reset each day
			{
				State.cumV_daily = 0;
				State.cumVG_daily = 0;
			}
		}

		//Only run to a certain number of infections: ggilani 28/10/14
		if (P.LimitNumInfections) continueEvents = (State.cumI < P.MaxNumInfections);
		fprintf(stderr, "\r    t=%lg   %i    %i|%i    %i     %i   %i (%lg %lg %lg)   %lg    ", t, State.S, State.L, State.I, State.R, State.D, State.cumD, State.cumT, State.cumV, State.cumVG, sqrt(State.maxRad2) / 1000); //added State.cumVG
		for (j = 0; ((j < P.UpdatesPerSample) && (!InterruptRun) && (continueEvents)); j++)
		{
			ts = (unsigned short int) (P.TimeStepsPerDay * t);

			//check to see if interruptions to interventions are in place
			if (P.DoInterruptIntervention)
			{
				P.InterruptIntervention = 0;
				for (in = 0; in < P.NDaysInterrupt; in++)
				{
					if (P.DaysInterruptIntervention[in] == (int)t)
					{
						P.InterruptIntervention = 1;
					}
				}
			}
			else
			{
				P.InterruptIntervention = 0;
			}

			//if we are to reset random numbers after an intervention event, specific time
			//if we are to reset random numbers after an intervention event, specific time
			if (P.ResetSeedsPostIntervention)
			{
				if ((P.ResetSeedsFlag == 0) && (ts >= (P.TimeToResetSeeds * P.TimeStepsPerDay)))
				{
					P.newseed1 = (int)(ranf() * 1e8);
					P.newseed2 = (int)(ranf() * 1e8);
					setall(P.newseed1, P.newseed2);
					P.ResetSeedsFlag = 1;
				}
			}

			if (fs)
			{
				if (P.DoAirports) TravelDepartSweep(t);
				k = (int)t;
				if (P.DurImportTimeProfile > 0)
				{
					if (k < P.DurImportTimeProfile)
						ir = P.ImportInfectionTimeProfile[k] * ((t > P.InfectionImportChangeTime) ? (P.InfectionImportRate2 / P.InfectionImportRate1) : 1.0);
					else
						ir = 0;
				}
				else
					ir = (t > P.InfectionImportChangeTime) ? P.InfectionImportRate2 : P.InfectionImportRate1;
#ifdef IMPORT_POP_SIZE
				ir *= ((double)P.N) / IMPORT_POP_SIZE;
#endif
				if (ir > 0)
				{
					ni = (int)ignpoi(P.TimeStep * ir);
					if (ni > 0)
					{
						SeedInfection(t, &ni, 1, run);
					}
				}
				if (P.FalsePositivePerCapitaIncidence > 0)
				{
					ni = (int)ignpoi(P.TimeStep * P.FalsePositivePerCapitaIncidence * ((double)P.N));
					if (ni > 0)
					{
						for (k = 0; k < ni; k++)
						{
							do
							{
								l = (int)(((double)P.N) * ranf());
							} while ((abs(Hosts[l].inf) == 5) || (ranf() > P.FalsePositiveAgeRate[HOST_AGE_GROUP(l)]));
							DoFalseCase(l, t, ts, 0);
						}
					}
				}
				InfectSweep(t, run); //adding run number as a parameter to infect sweep so we can track run number: ggilani - 15/10/14
				if (!P.DoSI) IncubRecoverySweep(t, run);
				// If we are doing hospitalisation by admin unit, process the hospitalisation queues filled during the incubation recovery sweep: ggilani - 24/11/14
				if ((P.DoHospitalisation) && (P.DoETUByAdUnit))
				{
					HospitalSweepAdunits(t);
				}
				// If doing new contact tracing, update numbers of people under contact tracing after each time step
				if ((P.DoContactTracing) && (t >= P.ContactTracingTimeStart))
				{
					ContactTracingSweep(t);
				}
				if (((P.DoRingVaccination) || (P.DoGeoVaccination)) && (t >= P.VaccTimeStart))
				{
					VaccSweep(t); //process ring vaccination queue
				}


				nu++;
				fs2 = ((P.DoDeath) || (P.DoSIS) || (State.L + State.I > 0) || (ir > 0) || (P.FalsePositivePerCapitaIncidence > 0));
				if (!TreatSweep(t))
				{
					if ((!fs2) && (State.L + State.I == 0) && (P.FalsePositivePerCapitaIncidence == 0)) { if ((ir == 0) && (((int)t) > P.DurImportTimeProfile)) fs = 0; }
				}
				if (P.DoAirports) TravelReturnSweep(t);
			}
			t += P.TimeStep;
			if (P.DoDeath) P.ts_age++;
			if ((P.DoSaveSnapshot) && (t <= P.SnapshotSaveTime) && (t + P.TimeStep > P.SnapshotSaveTime)) SaveSnapshot();
			if (t > P.TreatNewCoursesStartTime) P.TreatMaxCourses += P.TimeStep * P.TreatNewCoursesRate;
			if ((t > P.VaccNewCoursesStartTime) && (t < P.VaccNewCoursesEndTime))
			{
				if (P.DoVaccDailyReplenishment)
				{
					P.VaccMaxCourses += P.TimeStep * P.VaccNewCoursesRate;
				}
				else if (P.DoVaccBulkReplenishment)
				{
					if (t > P.VaccNewCoursesBoostStartTime)
					{

						P.VaccMaxCourses += P.VaccNewCoursesBulk;
						P.VaccNewCoursesBoostStartTime += P.VaccNewCoursesDelay;

					}
					else if (t > P.VaccNewCoursesStartTime)
					{
						if ((P.VaccNewCoursesStartTime - P.VaccNewCoursesInitDelay) < 0)
						{
							vaccUsed = 0;

						}
						else
						{
							vaccUsed = ((TimeSeries[(int)(P.VaccNewCoursesStartTime)].cumV + TimeSeries[(int)(P.VaccNewCoursesStartTime)].cumVG) - (TimeSeries[(int)(P.VaccNewCoursesStartTime - P.VaccNewCoursesInitDelay)].cumV + TimeSeries[(int)(P.VaccNewCoursesStartTime - P.VaccNewCoursesInitDelay)].cumVG));

						}
						//vaccUsed = ((TimeSeries[(int)(P.VaccNewCoursesStartTime)].cumV + TimeSeries[(int)(P.VaccNewCoursesStartTime)].cumVG) - (TimeSeries[(int)(P.VaccNewCoursesStartTime-P.VaccNewCoursesInitDelay)].cumV + TimeSeries[(int)(P.VaccNewCoursesStartTime-P.VaccNewCoursesInitDelay)].cumVG));
						P.VaccMaxCourses += min(P.VaccNewCoursesInitBulk, vaccUsed);
						P.VaccNewCoursesStartTime += P.VaccNewCoursesInitDelay;
					}
				}
			}
			cI = ((double)(State.S)) / ((double)P.N);
		}

	}
	RecordSample(t, P.NumSamples - 1);
	fprintf(stderr, "\nEnd of run\n");
	t2 = t + P.SampleTime;
	while (fs)
	{
		fs = TreatSweep(t2);
		t2 += P.SampleStep;
	}
	//	fprintf(stderr,"End RunModel\n");
	if (P.DoAirports)
	{
		t2 = t;
		for (t2 = t; t2 <= t + MAX_TRAVEL_TIME; t2 += P.TimeStep)
			TravelReturnSweep(t2);
	}
	/*	fprintf(stderr,"Checking consistency of final state...\n");
		for(i=j=k=ni=fs2=i2=0;i<P.N;i++)
			{
			if(i%1000==0) fprintf(stderr,"\r*** %i              ",i);
			if(Hosts[i].inf==0) j++;
			if((Hosts[i].pcell<P.NC)&&(Hosts[i].pcell>=0))
				{
				if(Cells[Hosts[i].pcell].susceptible[Hosts[i].listpos]!=i)
					{
					k++;
					for(l=fs=0;(l<Cells[Hosts[i].pcell].n)&&(!fs);l++)
						fs=(Cells[Hosts[i].pcell].susceptible[l]==i);
					if(!fs) ni++;
					}
				else
					{
					if((Hosts[i].listpos>Cells[Hosts[i].pcell].S-1)&&(Hosts[i].inf==0)) i2++;
					}
				}
			else
				fs2++;
			}
		fprintf(stderr,"\n*** susceptibles=%i\tincorrect listpos=%i\thosts not found in cell list=%i\tincorrect cell refs=%i\tincorrect positioning in cell susc list=%i\n",j,k,ni,fs2,i2);
	*/	RecordInfTypes();
}



void TravelReturnSweep(double t)
{
	int i, j, k, l, n, nr, ner, tn;

	if (floor(1 + t + P.TimeStep) != floor(1 + t))
	{
		nr = ner = 0;
		k = (int)floor(t);
		l = 1 + k % MAX_TRAVEL_TIME;
#pragma omp parallel for private(i,j,k,n,tn) reduction(+:nr,ner) schedule(static,1)
		for (tn = 0; tn < P.NumThreads; tn++)
		{
			for (j = tn; j < P.Nplace[HOTEL_PLACE_TYPE]; j += P.NumThreads)
			{
				n = Places[HOTEL_PLACE_TYPE][j].n;
				for (k = n - 1; k >= 0; k--)
				{
					i = Places[HOTEL_PLACE_TYPE][j].members[k];
					if (Hosts[i].Travelling == l)
					{
						n--;
						/*						if((n<0)||(Places[HOTEL_PLACE_TYPE][j].members[n]<0)||(Places[HOTEL_PLACE_TYPE][j].members[n]>=P.N))
													{fprintf(stderr,"### %i %i %i %i\n",j,k,n,Places[HOTEL_PLACE_TYPE][j].members[n]);ner++;}
												else if((k<0)||(k>n))
													{fprintf(stderr,"@ %i %i %i %i\n",j,k,n,Places[HOTEL_PLACE_TYPE][j].members[n]);ner++;}
												else
						*/
						if (k != n)
						{
							Places[HOTEL_PLACE_TYPE][j].members[k] = Places[HOTEL_PLACE_TYPE][j].members[n];
						}
						nr++;
						if (Hosts[i].PlaceLinks[HOTEL_PLACE_TYPE] != j)
						{
							ner++; fprintf(stderr, "(%i %i) ", j, Hosts[i].PlaceLinks[HOTEL_PLACE_TYPE]);
						}
						Hosts[i].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
						Hosts[i].Travelling = 0;
					}
				}
				Places[HOTEL_PLACE_TYPE][j].n = n;
			}
		}
		fprintf(stderr, " d=%i e=%i>", nr, ner);
	}
}

void TravelDepartSweep(double t)
{
	int c, i, i2, j, k, l, d, d2, m, n, f, f2, f3, mps, nld, nad, nsk, tn, bm, hp;
	double s, s2, nl;
	cell* ct;

	if (floor(1 + t - P.TimeStep) != floor(1 + t))
	{
		bm = ((P.DoBlanketMoveRestr) && (t >= P.MoveRestrTimeStart) && (t < P.MoveRestrTimeStart + P.MoveRestrDuration));
		mps = 2 * ((int)P.PlaceTypeMeanSize[HOTEL_PLACE_TYPE]) - P.NumThreads - 1;
		k = (int)floor(t);
		d = k % MAX_TRAVEL_TIME;
		nad = nld = nsk = 0;
#pragma omp parallel for private(i,i2,j,k,l,d2,m,n,s,f,f2,f3,tn,hp) reduction(+:nad,nsk) schedule(static,1)
		for (tn = 0; tn < P.NumThreads; tn++)
			for (i = tn; i < P.Nairports; i += P.NumThreads)
				if ((Airports[i].total_traffic > 0) && (Airports[i].num_mcell > 0))
				{
					s = Airports[i].total_traffic;
					if ((t > P.AirportCloseTimeStart) && (t < P.AirportCloseTimeStart + P.AirportCloseTimeStartBase))
						s *= P.AirportCloseEffectiveness;
					n = (s > 0) ? ((int)ignpoi_mt((double)s, tn)) : 0;
					f3 = 0;
					j = 0;
					while (j < n)
					{
						s = ranf_mt(tn);
						l = Airports[i].Inv_DestMcells[(int)floor(s * 1024)];
						while (Airports[i].DestMcells[l].prob < s) l++;
						l = Airports[i].DestMcells[l].id;
						k = (int)(ranf_mt(tn) * ((double)Mcells[l].n));
						i2 = Mcells[l].members[k];
						if ((abs(Hosts[i2].inf) < 2) && (Hosts[i2].inf != -2))
						{
							d2 = HOST_AGE_GROUP(i2);
							if ((P.RelativeTravelRate[d2] == 1) || (ranf_mt(tn) < P.RelativeTravelRate[d2]))
							{
								f2 = 1;
#pragma omp critical
								{
									if (Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] == -1)
									{
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -2;
										f2 = 0;
									}
								}
								if (!f2)
								{
									s = ranf_mt(tn);
									l = Airports[i].Inv_prop_traffic[(int)floor(s * 128)];
									while (Airports[i].prop_traffic[l] < s) l++;
									k = Airports[i].conn_airports[l];
									f2 = 0;
									if (bm)
									{
										if (dist2_raw(Airports[i].loc_x, Airports[i].loc_y, Airports[k].loc_x, Airports[k].loc_y) > P.MoveRestrRadius2)
										{
											if (ranf_mt(tn) > P.MoveRestrEffect)
											{
												f2 = 1;
												nsk++;
												j++;
#pragma omp critical
												Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
											}
										}
									}
									if (!f2)
									{
										f = 1;
										f2 = 0;
										do
										{
											s = ranf_mt(tn);
											m = Airports[k].Inv_DestPlaces[(int)floor(s * 1024)];
											while (Airports[k].DestPlaces[m].prob < s) m++;
											l = Airports[k].DestPlaces[m].id;
#pragma omp critical
											{
												if ((hp = Places[HOTEL_PLACE_TYPE][l].n) < mps)
												{
													f = 0;
													Places[HOTEL_PLACE_TYPE][l].n++;
												}
											}
											if (!f)
											{
												f3 = 0;
												Places[HOTEL_PLACE_TYPE][l].members[hp] = i2;
												d2 = (d + P.InvJourneyDurationDistrib[(int)(ranf_mt(tn) * 1024.0)]) % MAX_TRAVEL_TIME;
												Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = l;
												Hosts[i2].Travelling = 1 + d2;
												nad++;
												j++;
											}
											f2++;
										} while ((f) && (f2 < 300));
										if (f)
										{
#pragma omp critical
											Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
											if (++f3 > 100)
											{
												j++; nsk++;
											}
										}
									}
								}
							}
						}
						else
							j++;
					}
				}
		fprintf(stderr, "<ar=%i as=%i", nad, nsk);
		nl = ((double)P.PlaceTypeMeanSize[HOTEL_PLACE_TYPE]) * P.HotelPropLocal / P.MeanLocalJourneyTime;
		nsk = 0;
#pragma omp parallel for private(c,i,i2,j,l,d2,m,n,s,s2,ct,f,f2,f3,tn,hp) reduction(+:nld,nsk) schedule(static,1)
		for (tn = 0; tn < P.NumThreads; tn++)
			for (i = tn; i < P.Nplace[HOTEL_PLACE_TYPE]; i += P.NumThreads)
			{
				c = ((int)(Places[HOTEL_PLACE_TYPE][i].loc_x / P.cwidth)) * P.nch + ((int)(Places[HOTEL_PLACE_TYPE][i].loc_y / P.cheight));
				n = (int)ignpoi_mt(nl * Cells[c].tot_prob, tn);
				if (Places[HOTEL_PLACE_TYPE][i].n + n > mps)
				{
					nsk += (Places[HOTEL_PLACE_TYPE][i].n + n - mps);
					n = mps - Places[HOTEL_PLACE_TYPE][i].n;
				}
				for (j = 0; j < n; j++)
				{
					do
					{
						f = 0;
						s = ranf_mt(tn);
						l = Cells[c].InvCDF[(int)floor(s * 1024)];
						while (Cells[c].cum_trans[l] < s) l++;
						ct = CellLookup[l];
						m = (int)(ranf_mt(tn) * ((double)ct->S0));
						if (m < (ct->S + ct->L))
						{
							i2 = ct->susceptible[m];
							d2 = HOST_AGE_GROUP(i2);
							f3 = 0;
							if ((Hosts[i2].Travelling == 0) && ((P.RelativeTravelRate[d2] == 1) || (ranf_mt(tn) < P.RelativeTravelRate[d2])))
							{
#pragma omp critical
								{if (Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] == -1) { Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -2; f3 = 1; }}
							}
							if (f3)
							{
								s2 = dist2_raw(Households[Hosts[i2].hh].loc_x, Households[Hosts[i2].hh].loc_y, Places[HOTEL_PLACE_TYPE][i].loc_x, Places[HOTEL_PLACE_TYPE][i].loc_y);
								f2 = 1;
								if ((bm) && (s2 > P.MoveRestrRadius2))
								{
									if (ranf_mt(tn) >= P.MoveRestrEffect)
									{
#pragma omp critical
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
										nsk++;
										f2 = 0;
									}
								}
								if (f2)
								{
									s = numKernel(s2) / Cells[c].max_trans[l];
									if (ranf_mt(tn) >= s)
									{
#pragma omp critical
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
										f = 1;
									}
									else
									{
										d2 = (d + P.InvLocalJourneyDurationDistrib[(int)(ranf_mt(tn) * 1024.0)]) % MAX_TRAVEL_TIME;
										hp = Places[HOTEL_PLACE_TYPE][i].n;
										Places[HOTEL_PLACE_TYPE][i].n++;
										Places[HOTEL_PLACE_TYPE][i].members[hp] = i2;
										Hosts[i2].Travelling = 1 + d2;
										nld++;
#pragma omp critical
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = i;
									}
								}
							}
							else
								f = 1;
						}
						else
							nsk++;
					} while (f);
				}
			}
		fprintf(stderr, " l=%i ls=%i ", nld, nsk);
	}
}


/* function: HospitalSweepAdunits
 * Processes the hospital queues filled during IncubationRecoverySweep
 *
 * Parameters:
 *	double t, the current time in the simulation
 *
 * Returns: void
 *
 * Author: ggilani; Date: 24/11/14
 *
 */
void HospitalSweepAdunits(double t)
{
	int i, j, k, l, tn, ad;
	int numBeds, numBeds_inAU, numBeds_outAU;
	int numFreeBeds, numFreeBeds_inAU, numFreeBeds_outAU;
	int flag, ntrials;
	int max_trials = 50; //maximum number of times to try and fill outside of admin unit hospital beds
	int age; //added age to help track hospitalisations by age.
	int ts;

	ts = (int)t * P.TimeStepsPerDayInt;

	//combine admission lists across threads into a single list per admin unit
	for (i = 0; i < P.NumAdunits; i++)
	{
		for (j = 0; j < P.NumThreads; j++)
		{
			for (k = 0; k < StateT[j].nh_queue[i]; k++)
			{
				AdUnits[i].h_queue[k + AdUnits[i].nh_queue] = StateT[j].h_queue[i][k];
			}
			AdUnits[i].nh_queue += StateT[j].nh_queue[i];
		}
	}

	//This is for ETUs which is on an admin unit basis and will only start a certain time after the outbreak has been detected

	if (t >= P.ETUTimeStart)
	{

		//process hospitalisation discharges, admissions here, outside of main parallel loop?
		//take care of discharges first
		for (i = 0; i < P.NumAdunits; i++)
		{
			for (j = 0; j < P.NumThreads; j++)
			{
				for (k = 0; k < StateT[j].nhd_queue[i]; k++)
				{
					if (Hosts[StateT[j].hd_queue[i][k]].etu)
					{
						Hosts[StateT[j].hd_queue[i][k]].etu = 0;
						//if host was hospitalised in the same admin unit of which they are a residence, we reduce the number of currently in-use within admin unit beds in the current admin unit
						AdUnits[i].currentETUBeds--;
						StateT[j].ETU_adunit[i]--;
					}
				}
			}
		}

		//if interventions are interrupted, no new cases go into hospital
		if (!P.InterruptIntervention)
		{
			//process admissions
			//#pragma omp parallel for private(i,tn,numFreeBeds,j,numBeds,numBeds_inAU) schedule(static,1)
			for (tn = 0; tn < P.NumThreads; tn++)
			{
				for (i = tn; i < P.NumAdunits; i += P.NumThreads)
				{
					numBeds = BedsAvailablePerAdUnit(t, i);
					numFreeBeds = numBeds - (AdUnits[i].currentETUBeds);

					//add something here to check capacity of beds and change flag for that admin unit if necessary
					if ((numBeds > 0) && (AdUnits[i].nh_queue > numFreeBeds) && (StateT[tn].capETU_adunit[i] == 0))
					{
						StateT[tn].capETU_adunit[i] = 1;
					}

					if (numFreeBeds != 0)
					{
						//find number of beds and number of free beds for within admin unit cases

						if (AdUnits[i].nh_queue <= numFreeBeds)
						{
							for (j = 0; j < AdUnits[i].nh_queue; j++)
							{
								age = HOST_AGE_GROUP(AdUnits[i].h_queue[j]);
								Hosts[AdUnits[i].h_queue[j]].etu = Hosts[AdUnits[i].h_queue[j]].recovery_time;
								//set the admin unit identifier in which they are hospitalised
								AdUnits[i].currentETUBeds++;
								StateT[tn].ETU_adunit[i]++;
								StateT[tn].cumETU_adunit[i]++;
								StateT[tn].cumETU++;
								StateT[tn].cumETUa[age]++; //added hosp by age: ggilani 23/02/22
							}
							AdUnits[i].nh_queue = 0;
						}
						else
						{
							SampleWithoutReplacement(tn, numFreeBeds, AdUnits[i].nh_queue);
							for (j = 0; j < numFreeBeds; j++)
							{
								age = HOST_AGE_GROUP(AdUnits[i].h_queue[SamplingQueue[tn][j]]);
								Hosts[AdUnits[i].h_queue[SamplingQueue[tn][j]]].etu = Hosts[AdUnits[i].h_queue[SamplingQueue[tn][j]]].recovery_time;
								//Hosts[AdUnits[i].h_queue[j]].hospitalised=Hosts[AdUnits[i].h_queue[j]].recovery_time;
								AdUnits[i].currentETUBeds++;
								StateT[tn].ETU_adunit[i]++;
								StateT[tn].cumETU_adunit[i]++;
								StateT[tn].cumETU++;
								StateT[tn].cumETUa[age]++; //added hosp by age: ggilani 23/02/22
							}
						}
					}
				}
			}
		}
	}

	if ((P.DoHospitalisation) & (P.IncludeHospitalPlaceType))
	{

		//This is for individual hospital places. Before ETUs are available, and when they are full, cases will attempt to access these hospital beds.
		//process hospitalisation discharges, admissions here, outside of main parallel loop?
		//take care of discharges first
		for (i = 0; i < P.Nplace[P.HospPlaceTypeNum]; i++)
		{
			for (j = Places[P.HospPlaceTypeNum][i].nhcws; j < Places[P.HospPlaceTypeNum][i].n_current; j++)
			{
				k = Places[P.HospPlaceTypeNum][i].members[j];
				if (Hosts[k].hospitalised == ts)
				{
					Hosts[k].hospitalised = 0;
					//StateT[j].H_adunit[i]--; - not currently working as not threaded by adunit
					//if host was hospitalised in the same admin unit of which they are a residence, we reduce the number of currently in-use within admin unit beds in the current admin unit
					Places[P.HospPlaceTypeNum][i].n_current--;
					//swap people round
					if (!Hosts[k].hcw)
					{
						l = Places[P.HospPlaceTypeNum][i].members[Places[P.HospPlaceTypeNum][i].n_current];
						Places[P.HospPlaceTypeNum][i].members[Places[P.HospPlaceTypeNum][i].n_current] = j;
						Places[P.HospPlaceTypeNum][i].members[j] = l;
					}
				}
			}
		}

		//Now admit people to hospitals
		for (tn = 0; tn < P.NumThreads; tn++)
		{
			for (i = tn; i < P.NumAdunits; i += P.NumThreads)
			{
				for (j = 0; j < AdUnits[i].nh_queue; j++)
				{
					k = AdUnits[i].h_queue[j];
					if (Hosts[k].etu == 0)
					{
						if (Places[P.HospPlaceTypeNum][Hosts[k].PlaceLinks[P.HospPlaceTypeNum]].n_current < Places[P.HospPlaceTypeNum][Hosts[k].PlaceLinks[P.HospPlaceTypeNum]].maxcapacity)
						{
							if (!Hosts[k].hcw)
							{
								Places[P.HospPlaceTypeNum][Hosts[k].PlaceLinks[P.HospPlaceTypeNum]].members[Places[P.HospPlaceTypeNum][Hosts[k].PlaceLinks[P.HospPlaceTypeNum]].n_current] = k;
							}
							Places[P.HospPlaceTypeNum][Hosts[k].PlaceLinks[P.HospPlaceTypeNum]].n_current++;
							Hosts[k].hospitalised = Hosts[k].recovery_time;
							age = HOST_AGE_GROUP(AdUnits[i].h_queue[j]);
							//StateT[tn].H_adunit[i]++;
							StateT[tn].cumH_adunit[i]++;
							StateT[tn].cumH++;
							StateT[tn].cumHa[age]++;
						}
					}
				}
			}
		}
	}

}

/* function HospitalSweep
 *
 * Purpose: to sweep through all hosts currently in hospitalised and discharge them if the current time is past their recovery/death time (stored in hospitalised variable)
 *
 * Input: t, current time
 * Returns: void
 *
 * Author: ggilani, Date: 30/10/14
 */
 /*void HospitalSweep(double t)
 {
	 int tn,nDischarge,nErr,j,k,nPatients,patient;

	 //Discharges first
	 nr=ner=0;
 #pragma omp parallel for private(j,k,nPatients,patient,tn) reduction(+:nDischarge,nErr) schedule(static,1)
	 for(tn=0;tn<P.NumThreads;tn++)
	 {
		 for(j=tn;j<P.Nplace[HOSP_PLACE_TYPE];j+=P.NumThreads)
		 {
			 nPatients=Places[HOTEL_PLACE_TYPE][j].n;
			 for(k=nPatients-1;k>=0;k--)
			 {
				 patient=Places[HOSP_PLACE_TYPE][j].members[k];
				 if(Hosts[patient].hospitalised==(int)(t*P.TimeStepsPerDay))
				 {
					 nPatients--;
					 if(k!=nPatients)
						 {Places[HOSP_PLACE_TYPE][j].members[k]=Places[HOSP_PLACE_TYPE][j].members[nPatients];}
					 nDischarge++;
					 if(Hosts[patient].PlaceLinks[HOSP_PLACE_TYPE]!=j)
						 {ner++;fprintf(stderr,"(%i %i) ",j,Hosts[patient].PlaceLinks[HOTEL_PLACE_TYPE]);}
					 Hosts[patient].PlaceLinks[HOSP_PLACE_TYPE]=-1;
					 Hosts[patient].hospitalised=0;
				 }
			 }
			 Places[HOSP_PLACE_TYPE][j].n=nPatients;
		 }
	 }
	 fprintf(stderr," d=%i e=%i>",nDischarge,nErr);

	 //Admissions second

 }*/

void InfectSweep(double t, int run) //added run number as argument in order to record it in event log
{
	int i, j, k, l, m, n, i2, b, i3, f, f2, tn, cq, bm, ci;
	double seasonality, sbeta, hbeta, s, s2, s3, s4, s5, fp, contact_scale;
	cell* c, * ct;
	microcell* mi, * mt, * mp;
	unsigned short int ts;
	person* si;

	if (!P.DoSeasonality)
		seasonality = 1.0;
	else
		seasonality = P.Seasonality[((int)t) % DAYS_PER_YEAR];
	ts = (unsigned short int) (P.TimeStepsPerDay * t);
	fp = P.TimeStep / (1 - P.FalsePositiveRate);
	sbeta = seasonality * fp * P.LocalBeta;
	hbeta = (P.DoHouseholds) ? (seasonality * fp * P.HouseholdTrans) : 0;
	bm = ((P.DoBlanketMoveRestr) && (t >= P.MoveRestrTimeStart) && (t < P.MoveRestrTimeStart + P.MoveRestrDuration));
#pragma omp parallel for private(j,k,l,m,n,i2,b,i3,f,f2,s,s2,s3,s4,s5,c,ct,mi,mt,mp,cq,ci,si,contact_scale) schedule(static,1)
	for (tn = 0; tn < P.NumThreads; tn++)
		for (b = tn; b < P.NCP; b += P.NumThreads)
		{
			c = CellLookup[b];
			s5 = 0;
			for (j = 0; j < c->I; j++)
			{
				ci = c->infected[j];
				si = Hosts + ci;

				if (hbeta > 0)
				{
					if ((Households[si->hh].nh > 1) && (!(si->nc_plus_hh_disabled & HH_DISABLED)) && (!si->Travelling))
					{
						l = Households[si->hh].FirstPerson;
						m = l + Households[si->hh].nh;
						s3 = hbeta * CalcHouseInf(ci, ts);
						f = 0;
						for (i3 = l; (i3 < m) && (!f); i3++)
							for (i2 = 0; (i2 < P.PlaceTypeNum) && (!f); i2++)
								if (Hosts[i3].PlaceLinks[i2] >= 0)
								{
									f = PLACE_CLOSED(i2, Hosts[i3].PlaceLinks[i2]);
								}
						if (f) { s3 *= P.PlaceCloseHouseholdRelContact; }/* NumPCD++;}*/
						for (i3 = l; i3 < m; i3++)
						{
							if ((Hosts[i3].inf == 0) && (!(Hosts[i3].nc_plus_hh_disabled & HH_DISABLED)) && (!Hosts[i3].Travelling))
							{
								s = s3 * CalcHouseSusc(i3, ts, ci, tn);
								if (ranf_mt(tn) < s)
								{
									cq = Hosts[i3].pcell % P.NumThreads;
									if ((StateT[tn].n_queue[cq] < P.InfQueuePeakLength)) //(Hosts[i3].infector==-1)&&
									{
										if ((P.FalsePositiveRate > 0) && (ranf_mt(tn) < P.FalsePositiveRate))
											Hosts[i3].infector = Hosts[i3].infect_type = -1;
										else
										{
											Hosts[i3].infector = ci;
											Hosts[i3].infect_type = 1 + INFECT_TYPE_MASK * (1 + si->infect_type / INFECT_TYPE_MASK);
										}
										StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++] = i3;
									}
								}
							}
						}
					}
				}
				if (P.DoPlaces)
				{
					if (!HOST_ABSENT(ci))
					{
						mi = Mcells + si->mcell;
						for (k = 0; k < P.PlaceTypeNum; k++)
						{
							l = si->PlaceLinks[k];
							if (k == P.HospPlaceTypeNum && Hosts[ci].hospitalised == 0)
							{
								l = -1;
							}
							if (Hosts[ci].hospitalised)
							{
								contact_scale = 0.5;
							}
							else
							{
								contact_scale = 1.0;
							}
							if ((l >= 0) && (!PLACE_CLOSED(k, l)))
							{
								s3 = fp * seasonality * contact_scale * CalcPlaceInf(ci, k, ts); //scale contacts so that hospitalised cases split their contacts between network and staff
								mp = Mcells + Places[k][l].mcell;
								if (bm)
								{
									if ((dist2_raw(Households[si->hh].loc_x, Households[si->hh].loc_y,
										Places[k][l].loc_x, Places[k][l].loc_y) > P.MoveRestrRadius2))
										s3 *= P.MoveRestrEffect;
								}
								else if ((mi->moverest != mp->moverest) && ((mi->moverest == 2) || (mp->moverest == 2)))
									s3 *= P.MoveRestrEffect;

								if ((k != HOTEL_PLACE_TYPE) && (!si->Travelling))
								{
									i2 = (si->PlaceGroupLinks[k]);
									s4 = s3;
									//s4=s3*(1-P.PlaceTypePropBetweenGroupLinks[k]*P.PlaceTypeGroupSizeParam1[k]/((double) Places[k][l].n));
									if (s4 < 0)
									{
										fprintf(stderr, "@@@ %lg\n", s4);
										exit(1);
									}
									else if (s4 >= 1)
									{
										if (k == P.HospPlaceTypeNum)
										{
											n = Places[k][l].nhcws;
										}
										else
										{
											n = Places[k][l].group_size[i2];
										}
									}
									else
									{
										if (k == P.HospPlaceTypeNum)
										{
											n = (int)ignbin_mt((long)Places[k][l].nhcws, s4, tn);
										}
										else
										{
											n = (int)ignbin_mt((long)Places[k][l].group_size[i2], s4, tn);
										}
									}
									if (n > 0)
									{
										if (k == P.HospPlaceTypeNum)
										{
											SampleWithoutReplacement(tn, n, Places[k][l].nhcws);
										}
										else
										{
											SampleWithoutReplacement(tn, n, Places[k][l].group_size[i2]);
										}
									}
									for (m = 0; m < n; m++)
									{
										if (k == P.HospPlaceTypeNum)
										{
											i3 = Places[k][l].members[SamplingQueue[tn][m]];
										}
										else
										{
											i3 = Places[k][l].members[Places[k][l].group_start[i2] + SamplingQueue[tn][m]];
										}
										if ((Hosts[i3].inf == 0) && (!HOST_ABSENT(i3)))
										{
											mt = Mcells + Hosts[i3].mcell;
											ct = Cells + Hosts[i3].pcell;
											s = CalcPlaceSusc(i3, k, ts, ci, tn);
											if (bm)
											{
												if ((dist2_raw(Households[Hosts[i3].hh].loc_x, Households[Hosts[i3].hh].loc_y,
													Places[k][l].loc_x, Places[k][l].loc_y) > P.MoveRestrRadius2))
													s *= P.MoveRestrEffect;
											}
											else if ((mt->moverest != mp->moverest) && ((mt->moverest == 2) || (mp->moverest == 2)))
												s *= P.MoveRestrEffect;
											if ((s == 1) || (ranf_mt(tn) < s))
											{
												cq = Hosts[i3].pcell % P.NumThreads;
												if ((StateT[tn].n_queue[cq] < P.InfQueuePeakLength)) //(Hosts[i3].infector==-1)&&
												{
													if ((P.FalsePositiveRate > 0) && (ranf_mt(tn) < P.FalsePositiveRate))
														Hosts[i3].infector = Hosts[i3].infect_type = -1;
													else
													{
														Hosts[i3].infector = ci;
														Hosts[i3].infect_type = 2 + k + INFECT_TYPE_MASK * (1 + si->infect_type / INFECT_TYPE_MASK);
													}
													StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++] = i3;
												}
											}
										}
									}
								}
								if ((k == HOTEL_PLACE_TYPE) || (!si->Travelling))
								{
									s3 *= P.PlaceTypePropBetweenGroupLinks[k] * P.PlaceTypeGroupSizeParam1[k] / ((double)Places[k][l].n);
									if (s3 < 0)
									{
										fprintf(stderr, "@@@ %lg\n", s3);
										exit(1);
									}
									else if (s3 >= 1)
										n = Places[k][l].n;
									else
										n = (int)ignbin_mt((long)Places[k][l].n, s3, tn);
									if (n > 0) SampleWithoutReplacement(tn, n, Places[k][l].n);
									for (m = 0; m < n; m++)
									{
										i3 = Places[k][l].members[SamplingQueue[tn][m]];
										if ((Hosts[i3].inf == 0) && (!HOST_ABSENT(i3)))
										{
											mt = Mcells + Hosts[i3].mcell;
											ct = Cells + Hosts[i3].pcell;
											s = CalcPlaceSusc(i3, k, ts, ci, tn);
											if (bm)
											{
												if ((dist2_raw(Households[Hosts[i3].hh].loc_x, Households[Hosts[i3].hh].loc_y,
													Places[k][l].loc_x, Places[k][l].loc_y) > P.MoveRestrRadius2))
													s *= P.MoveRestrEffect;
											}
											else if ((mt->moverest != mp->moverest) && ((mt->moverest == 2) || (mp->moverest == 2)))
												s *= P.MoveRestrEffect;
											if ((s == 1) || (ranf_mt(tn) < s))
											{
												cq = Hosts[i3].pcell % P.NumThreads;
												if ((StateT[tn].n_queue[cq] < P.InfQueuePeakLength))//(Hosts[i3].infector==-1)&&
												{
													if ((P.FalsePositiveRate > 0) && (ranf_mt(tn) < P.FalsePositiveRate))
														Hosts[i3].infector = Hosts[i3].infect_type = -1;
													else
													{
														Hosts[i3].infector = ci;
														Hosts[i3].infect_type = 2 + k + NUM_PLACE_TYPES + INFECT_TYPE_MASK * (1 + si->infect_type / INFECT_TYPE_MASK);
													}
													StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++] = i3;
												}
											}
										}
									}
								}
							}
						}
					}
				}
				if (sbeta > 0)
				{
					if (si->Travelling)
					{
						s2 = 0; f = 0;
					}
					else
					{
						s2 = CalcSpatialInf(ci, ts);
						f = 0;
						if (P.DoPlaces)
							for (i3 = 0; (i3 < P.PlaceTypeNum) && (!f); i3++)
								if (si->PlaceLinks[i3] >= 0)
								{
									f = PLACE_CLOSED(i3, si->PlaceLinks[i3]);
								}
					}
					if (f)
					{
						s2 *= P.PlaceCloseSpatialRelContact;
						/* NumPCD++; */
						s5 += s2;
						StateT[tn].cell_inf[j] = -s5;
					}
					else
					{
						s5 += s2;
						StateT[tn].cell_inf[j] = s5;
					}
				}
			}
			if (s5 > 0)
			{
				n = (int)ignpoi_mt(s5 * sbeta * ((double)c->tot_prob), tn);
				i2 = c->I;
				if (n > 0)
				{
					for (j = 0; j < i2 - 1; j++) StateT[tn].cell_inf[j] /= s5;
					StateT[tn].cell_inf[i2 - 1] = (StateT[tn].cell_inf[i2 - 1] < 0) ? -1 : 1;
				}
				for (k = 0; k < n; k++)
				{
					if (i2 == 1)
						j = 0;
					else
					{
						s = ranf_mt(tn);
						j = m = i2 / 2;
						f = 1;
						do
						{
							if (m > 1) m /= 2;
							if ((j > 0) && (fabs(StateT[tn].cell_inf[j - 1]) >= s))
							{
								j -= m;
								if (j == 0)
									f = 0;
								//									else
								//										f=(fabs(StateT[tn].cell_inf[j-1])>=s);
							}
							else if ((j < i2 - 1) && (fabs(StateT[tn].cell_inf[j]) < s))
							{
								j += m;
								if (j == i2 - 1)
									f = 0;
								//									else
								//										f=(fabs(StateT[tn].cell_inf[j])<s);
							}
							else
							{
								/*									fprintf(stderr,"##%i %i %i %i: %lg %lg %lg\n",i2,j,m,f,s,fabs(StateT[tn].cell_inf[j-1]),fabs(StateT[tn].cell_inf[j]));
								*/									f = 0;
							}
						} while (f);
						/*							if(fabs(StateT[tn].cell_inf[j])<s) fprintf(stderr,"##%i %i %i %i: %lg %lg %lg\n",i2,j,m,f,s,fabs(StateT[tn].cell_inf[j-1]),fabs(StateT[tn].cell_inf[j]));
												if((j>0)&&(fabs(StateT[tn].cell_inf[j-1])>=s)) fprintf(stderr,"@@%i %i %i %i: %lg %lg %lg\n",i2,j,m,f,s,fabs(StateT[tn].cell_inf[j-1]),fabs(StateT[tn].cell_inf[j]));
						*/
					}
					f = (StateT[tn].cell_inf[j] < 0);
					ci = c->infected[j];
					si = Hosts + ci;
					do
					{
						s = ranf_mt(tn);
						l = c->InvCDF[(int)floor(s * 1024)];
						while (c->cum_trans[l] < s) l++;
						ct = CellLookup[l];
						m = (int)(ranf_mt(tn) * ((double)ct->S0));
						i3 = ct->susceptible[m];
						s2 = dist2(Hosts + i3, Hosts + ci);
						s = numKernel(s2) / c->max_trans[l];
						//alter acceptance probability for cross border effect here: testing the cross border effect - ggilani 19/01/15
						if ((AdUnits[Mcells[Hosts[ci].mcell].adunit].id / P.CountryDivisor) != (AdUnits[Mcells[Hosts[i3].mcell].adunit].id / P.CountryDivisor))
						{
							s *= P.PropCrossBorderInf;
						}
						f2 = 0;
						if ((ranf_mt(tn) >= s) || (abs(Hosts[i3].inf) == 5))
						{
							f2 = 1;
						}
						//add cross border effect - relative contact assumed to be less if infector and infectee live in different countries: ggilani 09/12/14
						//else if(((AdUnits[Mcells[Hosts[ci].mcell].adunit].id/P.CountryDivisor)!=(AdUnits[Mcells[Hosts[i3].mcell].adunit].id/P.CountryDivisor))&&(ranf_mt(tn)>=P.PropCrossBorderInf)) //checking to see if they are in the same country
						//{
						//	f2=1;
						//}
						else if (m < ct->S)
						{
							if ((!Hosts[i3].Travelling) && ((c != ct) || (Hosts[i3].hh != si->hh)))
							{
								mi = Mcells + si->mcell;
								mt = Mcells + Hosts[i3].mcell;
								s = CalcSpatialSusc(i3, ts, ci, tn);
								if (bm)
								{
									if ((dist2_raw(Households[si->hh].loc_x, Households[si->hh].loc_y,
										Households[Hosts[i3].hh].loc_x, Households[Hosts[i3].hh].loc_y) > P.MoveRestrRadius2))
										s *= P.MoveRestrEffect;
								}
								else if ((mt->moverest != mi->moverest) && ((mt->moverest == 2) || (mi->moverest == 2)))
									s *= P.MoveRestrEffect;
								if (!f)
								{
									for (m = f2 = 0; (m < P.PlaceTypeNum) && (!f2); m++)
										if (Hosts[i3].PlaceLinks[m] >= 0)
										{
											f2 = PLACE_CLOSED(m, Hosts[i3].PlaceLinks[m]);
										}
									if (f2) { s *= P.PlaceCloseSpatialRelContact; }/* NumPCD++;} */
									f2 = 0;
								}
								if ((s == 1) || (ranf_mt(tn) < s))
								{
									cq = ((int)(ct - Cells)) % P.NumThreads;
									if ((Hosts[i3].inf == 0) && (StateT[tn].n_queue[cq] < P.InfQueuePeakLength)) //Hosts[i3].infector==-1
									{
										if ((P.FalsePositiveRate > 0) && (ranf_mt(tn) < P.FalsePositiveRate))
											Hosts[i3].infector = Hosts[i3].infect_type = -1;
										else
										{
											Hosts[i3].infector = ci;
											Hosts[i3].infect_type = 2 + 2 * NUM_PLACE_TYPES + INFECT_TYPE_MASK * (1 + si->infect_type / INFECT_TYPE_MASK);
										}
										StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++] = i3;
									}
								}
							}
						}
					} while (f2);
				}
			}
		}


#pragma omp parallel for private(i,k) schedule(static,1)
	for (j = 0; j < P.NumThreads; j++)
	{
		for (k = 0; k < P.NumThreads; k++)
		{
			for (i = 0; i < StateT[k].n_queue[j]; i++)
			{
				if (Hosts[StateT[k].inf_queue[j][i]].infect_type == -1)
					DoFalseCase(StateT[k].inf_queue[j][i], t, ts, j);
				else
					DoInfect(StateT[k].inf_queue[j][i], t, j, run);
			}
			StateT[k].n_queue[j] = 0;
		}
	}
}

void IncubRecoverySweep(double t, int run)
{
	int i, j, k, l, b, tn, ci, day;
	double currInfSafeFuneral;
	cell* c;
	person* si;
	unsigned short int ts, tc;
	double currPropSafeFunerals, currRelInfFuneral;

	ts = (unsigned short int) (P.TimeStepsPerDay * t);

	//if we're doing hospitalisation by admin unit, set number of people in discharge and admission queues to zero; otherwise if doing hospitalisation by place, only set state hospitalisation queue to zero
	if (P.DoHospitalisation)
	{
		if (P.DoETUByAdUnit)
		{
			for (i = 0; i < P.NumAdunits; i++)
			{
				AdUnits[i].nh_queue = 0;
				for (j = 0; j < P.NumThreads; j++)
				{
					StateT[j].nh_queue[i] = 0;
					StateT[j].nhd_queue[i] = 0;
				}
			}
		}
		else
			//If we're doing hospitalisation by place, we'll still use the same queue to store individuals who need hospitalisation, but we'll store them all in the same queue, not by admin unit
		{
			for (j = 0; j < P.NumThreads; j++)
			{
				StateT[j].nhd_queue[0] = 0;
			}
		}
	}

#pragma omp parallel for private(j,k,l,b,c,tn,tc,ci,si,day) schedule(static,1)
	for (tn = 0; tn < P.NumThreads; tn++)
	{
		for (b = tn; b < P.NCP; b += P.NumThreads)
		{
			c = CellLookup[b];
			for (j = ((int)c->L - 1); j >= 0; j--)
				if (ts >= Hosts[c->latent[j]].latent_time) DoIncub(c->latent[j], ts, tn, run);
			StateT[tn].n_queue[0] = 0;
			for (j = c->I - 1; j >= 0; j--)
			{
				ci = c->infected[j];
				si = Hosts + ci;
				tc = si->latent_time + ((int)(P.LatentToSymptDelay / P.TimeStep));
				if ((P.DoSymptoms) && (ts == tc))
					DoCase(ci, t, ts, tn);
				if ((ts == si->detect_time) && (si->detected))
				{
					DoDetectedCase(ci, t, ts, tn);
				}

				//Now considered which of infected have reached time to hospitalisation
				if ((P.DoHospitalisation) && ((ts >= si->hospital_time) && (ts < (si->hospital_time + (int)(P.HospWaitingTime * P.TimeStepsPerDay)))) && (abs(si->inf) != 6) && (!si->hospitalised) && (!si->etu) && (si->detected))
					//A lot of conditions! To enter loop, we must be doing hospitalisation by admin unit, the host must be within their waiting hospitalisation time, not already hospitalised and not dead but infectious! Must also be detected case
				{
					//mark someone to be admitted
					if (P.DoETUByAdUnit)
					{
						StateT[tn].h_queue[Mcells[si->mcell].adunit][StateT[tn].nh_queue[Mcells[si->mcell].adunit]++] = ci;
					}
					else
					{
						StateT[tn].h_queue[0][StateT[tn].nh_queue[0]++] = ci;
					}
				}


				//Adding code to assign recovery or death when leaving the infectious class: ggilani - 22/10/14
				if (ts >= si->recovery_time)
				{
					if ((P.DoHospitalisation) && (P.DoETUByAdUnit) && ((si->hospitalised) || (si->etu)))
					{
						//if someone has reached their recovery time, regardless of whether it's a death or recovery, they will leave hospital at this point if they are hospitalised
						//mark someone to be discharged; however, if we're doing hospitalisation by place, this gets taken care of in HospitalSweep
						StateT[tn].hd_queue[Mcells[si->mcell].adunit][StateT[tn].nhd_queue[Mcells[si->mcell].adunit]++] = ci;
					}
					if ((si->to_die) && (P.DoFuneralTransmission) && (abs(si->inf) != 6)) //only go through this if block if using waiting time until leaving infectious class, and if host is alive and infectious
					{
						//set infection status flag to dead and infectious: -6
						si->inf = -6;
						//do some accounting here to make sure we count the death at the right time
						StateT[tn].cumD++;
						StateT[tn].cumDa[HOST_AGE_GROUP(ci)]++;
						if (P.DoAdUnits) StateT[tn].cumD_adunit[Mcells[si->mcell].adunit]++;

						// set recovery time to current recovery time plus length of funeral transmission duration
						si->recovery_time += (unsigned short int)(P.FuneralTransmissionDuration * P.TimeStepsPerDay);

						if ((P.DoHospitalisation) && ((si->hospitalised) || (si->etu))) //this is regardless of whether we are considering hospitalisation by admin unit or place
						{
							//if in hospital, they are definitely detected when they die
							StateT[tn].cumDD++;
							if (P.DoAdUnits) StateT[tn].cumDD_adunit[Mcells[si->mcell].adunit]++;
							if ((t >= P.FuneralControlTimeStart) && (State.cumSDB_adunit[Mcells[si->mcell].adunit] < AdUnits[Mcells[si->mcell].adunit].maxSDB))
							{
								//if someone has died in hospital, we assumed that they will have a safe burial
								si->infectiousMult = (P.RelativeInfectiousnessFuneral * P.RelInfSafeFuneral);
								//if in hospital, they definitely have a safe burial
								StateT[tn].cumSDB++;
								if (P.DoAdUnits) StateT[tn].cumSDB_adunit[Mcells[si->mcell].adunit]++;
							}
						}
						else if (si->detected)
						{
							//if a detected case, they are definitely detected when they die
							StateT[tn].cumDD++;
							if (P.DoAdUnits) StateT[tn].cumDD_adunit[Mcells[si->mcell].adunit]++;
							//alter host's infectiousness, taking into account relative reduction in infectiousness due to safe burial
							if ((t >= P.FuneralControlTimeStart) && (ranf_mt(tn) <= P.ProportionSafeFuneral) && (State.cumSDB_adunit[Mcells[si->mcell].adunit] < AdUnits[Mcells[si->mcell].adunit].maxSDB))
							{
								//if safe burials in effect, they have a safe burial with probability ProportionSafeFuneral
								si->infectiousMult = (P.RelativeInfectiousnessFuneral * P.RelInfSafeFuneral);
								StateT[tn].cumSDB++;
								if (P.DoAdUnits) StateT[tn].cumSDB_adunit[Mcells[si->mcell].adunit]++;
							}
						}
						//if undetected, they don't have a safe burial
						else
						{
							//alter host's infectiousness
							si->infectiousMult = P.RelativeInfectiousnessFuneral;
							si->safeBurial = 0;
						}
						//}

					}
					if ((!si->to_die) && (ts >= si->recovery_time))
						StateT[tn].inf_queue[0][StateT[tn].n_queue[0]++] = ci;
					else if ((si->to_die) && (ts >= si->recovery_time))
					{
						if ((!HOST_TREATED(ci)) || (abs(si->inf) == 6)) // if infectious status is 6, then host is already dead and cannot be treated!
							DoDeath(ci, tn, run);
						else if (ranf_mt(tn) >= P.EvolResistRelTreatDeathDrop[si->resist])
							StateT[tn].inf_queue[0][StateT[tn].n_queue[0]++] = ci;
						else
							DoDeath(ci, tn, run);
					}

				}

			}

			for (j = 0; j < StateT[tn].n_queue[0]; j++) DoRecover(StateT[tn].inf_queue[0][j], run, tn);
			StateT[tn].n_queue[0] = 0;
			l = (P.EvolInfectMutationRate == 0) ? 0 : ((int)ignbin_mt((long)c->I, P.TimeStep * P.EvolInfectMutationRate, tn));
			for (j = 0; j < l; j++)
			{
				k = (int)(((double)c->I) * ranf_mt(tn));
				if (Hosts[c->infected[k]].base_inf_level < P.EvolInfectMax)
					Hosts[c->infected[k]].base_inf_level += P.EvolInfectStep;
			}
		}
	}

}

/*
 * Function: Update ContactTracing
 *
 * Purpose: to count the number of people in each admin unit who are being contact traced each day
 * Parameters: none
 * Returns: none
 *
 * Author: ggilani, 08/06/17
 */
void ContactTracingSweep(double t)
{

	int i, j, k, tn;
	unsigned short int ts;
	int numCT;

	//find current time step
	ts = (unsigned short int) (P.TimeStepsPerDay * t);

	//take care of people no longer requiring contact tracing first
	for (i = 0; i < P.NumAdunits; i++)
	{
		if (AdUnits[i].contactTraceThresholdCrossed == 1)
		{
			for (j = 0; j < AdUnits[i].nct;)
			{
				//if((j==0)&&(ts==293))
				//{
				//	fprintf(stderr,"check here\n");
				//}
				if (ts >= (Hosts[AdUnits[i].ct[j]].contactTraced_end_time))
				{
					//stop contact tracing this host
					Hosts[AdUnits[i].ct[j]].contactTraced = 0;
					k = AdUnits[i].ct[j];
					AdUnits[i].ct[j] = AdUnits[i].ct[AdUnits[i].nct - 1];
					AdUnits[i].ct[AdUnits[i].nct - 1] = k;
					AdUnits[i].nct--;
				}
				else
				{
					j++;
				}
			}
		}
	}

	//now sort out people in the queue wanting to join contact tracing
	//first combine individual thread lists of queues into complete admin unit lists
	for (i = 0; i < P.NumAdunits; i++)
	{
		for (j = 0; j < P.NumThreads; j++)
		{
			for (k = 0; k < StateT[j].nct_queue[i]; k++)
			{
				AdUnits[i].ct_queue[k + AdUnits[i].nct_queue] = StateT[j].ct_queue[i][k];
			}
			AdUnits[i].nct_queue += StateT[j].nct_queue[i];
		}
	}

	//if interventions interrupted, don't contact trace new contacts
	if (!P.InterruptIntervention)
	{

		//reorder queue? Haven't done this so far, as it doesn't make sense to me that not all contacts of same infected person would be traced - should this be done? 14/06/17

		//now contact trace people for whom there is space
		for (tn = 0; tn < P.NumThreads; tn++)
		{
			for (i = tn; i < P.NumAdunits; i += P.NumThreads)
			{

				if (AdUnits[i].contactTraceThresholdCrossed == 1)
				{
					//find current available capacity for contact tracing
					numCT = AdUnits[i].contactTraceCapacity - AdUnits[i].nct;

					for (j = 0; j < AdUnits[i].nct_queue; j++)
					{
						if (((numCT > 0) && (Hosts[AdUnits[i].ct_queue[j]].contactTraced == 0)))
						{
							Hosts[AdUnits[i].ct_queue[j]].contactTraced = 1; //set that they are being contact traced
							Hosts[AdUnits[i].ct_queue[j]].contactTraced_start_time = ts; //set the time at which contact tracing starts
							//add something for lost to follow up contacts - ggilani 11/05/22
							if (ranf_mt(tn) < P.propContactLost) //if this contact is going to be lost to follow up
							{
								Hosts[AdUnits[i].ct_queue[j]].contactTraced_end_time = ts + (unsigned short int) ((P.contactTraceDuration * ranf_mt(tn)) * P.TimeStepsPerDay);
							}
							else
							{
								Hosts[AdUnits[i].ct_queue[j]].contactTraced_end_time = ts + (unsigned short int) (P.contactTraceDuration * P.TimeStepsPerDay); //set the time at which contact tracing ends

							}AdUnits[i].ct[AdUnits[i].nct] = AdUnits[i].ct_queue[j]; //store the id of the host being contact traced in the array
							AdUnits[i].nct++; //increment the number of people being incremented
							StateT[tn].cumCT_adunit[i]++;
							StateT[tn].cumCT++;
							numCT--;
						}
					}

					//if((numCT!=0)&&(AdUnits[i].nct_queue>0))
					//{
					//	//if(i==0)
					//	//{
					//	//	fprintf(stderr,"check here\n");
					//	//}
					//	if(AdUnits[i].nct_queue<=numCT)
					//	{
					//		for(j=0;j<AdUnits[i].nct_queue;j++)
					//		{
					//			Hosts[AdUnits[i].ct_queue[j]].contactTraced=1; //set that they are being contact traced
					//			Hosts[AdUnits[i].ct_queue[j]].contactTraced_start_time=ts; //set the time at which contact tracing starts
					//			Hosts[AdUnits[i].ct_queue[j]].contactTraced_end_time=ts+(unsigned short int) (P.contactTraceDuration*P.TimeStepsPerDay); //set the time at which contact tracing ends
					//			AdUnits[i].ct[AdUnits[i].nct]=AdUnits[i].ct_queue[j]; //store the id of the host being contact traced in the array
					//			AdUnits[i].nct++; //increment the number of people being incremented
					//			StateT[tn].cumCT_adunit[i]++;
					//			StateT[tn].cumCT++;
					//		}
					//	}
					//	else
					//	{
					//		//SampleWithoutReplacement(tn,numFreeBeds_inAU,AdUnits[i].nh_queue);
					//		for(j=0;j<numCT;j++)
					//		{
					//			Hosts[AdUnits[i].ct_queue[j]].contactTraced=1; //set that they are being contact traced
					//			Hosts[AdUnits[i].ct_queue[j]].contactTraced_start_time=ts; //set the time at which contact tracing starts
					//			Hosts[AdUnits[i].ct_queue[j]].contactTraced_end_time=ts+(unsigned short int) (P.contactTraceDuration*P.TimeStepsPerDay);
					//			//Hosts[AdUnits[i].h_queue[SamplingQueue[tn][j]]].hospitalised=Hosts[AdUnits[i].h_queue[SamplingQueue[tn][j]]].recovery_time;
					//			AdUnits[i].ct[AdUnits[i].nct]=AdUnits[i].ct_queue[j]; //store the id of the host being contact traced in the array
					//			AdUnits[i].nct++; //increment the number of people being incremented
					//			StateT[tn].cumCT_adunit[i]++;
					//			StateT[tn].cumCT++;
					//		}
					//	}	
					//}
				}

			}
		}


	}

	//because people don't stay in the contact tracing queue over multiple days (they are either contact traced or they aren't), reset the number of people in the queues to zero after each sweep

	for (i = 0; i < P.NumAdunits; i++)
	{
		AdUnits[i].nct_queue = 0;
		for (j = 0; j < P.NumThreads; j++)
		{
			StateT[j].nct_queue[i] = 0;
		}
	}

	/*Older code which just counted the number of people needing contact tracing - but was very slow*/

	//int tn,i,j,k;

	////loop over admin units and reset counts of people being contact traced to zero
	//for(i=0;i<P.NumAdunits;i++)
	//{
	//	AdUnits[i].contactTraceCurrent=0;
	//}

	////loop over all people in the population and if they are currently being contact traced, add them to the contact tracing list for their admin unit
	//for(i=0;i<P.N;i++)
	//{
	//	if((Hosts[i].contactTraced>0)&&(Hosts[i].contactTraced<t))
	//	{
	//		AdUnits[Mcells[Hosts[i].mcell].adunit].contactTraceCurrent++;
	//	}
	//}

}

/*
* Function: Vaccination sweep
*
* Purpose: to process the ring and geo vaccination queue
* Parameters: none
* Returns: none
*
* Author: ggilani, 21/08/19
*/
void VaccSweep(double t)
{
	int i, j, k, l, tn, maxvacc;
	unsigned short int ts;
	int numCT;

	//find current time step
	ts = (unsigned short int) (P.TimeStepsPerDay * t);

	//combine queues for individual threads into one master list
	for (i = 0; i < P.NumThreads; i++)
	{
		for (j = 0; j < StateT[i].nvacc_queue; j++) // (j = 0; j < StateT[i].nringvacc_queue; j++) //for (j = StateT[i].ringvacc_cum; j < StateT[i].nringvacc_queue; j++)
		{

			if (State.nvacc_queue < P.N) //to stop overspill from the queue: ggilani
			{
				State.vacc_queue[State.nvacc_queue] = StateT[i].vacc_queue[j];
				State.ring_queue[State.nvacc_queue] = StateT[i].ring_queue[j];
				State.nvacc_queue++;
			}
		}
		StateT[i].nvacc_queue = 0; // StateT[i].ringvacc_cum = StateT[i].nringvacc_queue; //no need to update this?
	}

	//process the ring vaccination queue
	k = State.vacc_cum;
	l = State.vacc_ind;

	//add some code to see if interventions are being interrupted - if not, do some vaccinations. if yes, no vaccinations are performed but stay in the queue.
	if (!P.InterruptIntervention)
	{
		for (i = State.vacc_ind; i < State.nvacc_queue; i++)
		{
			if ((State.cumV + State.cumVG) < P.VaccMaxCourses)
			{
				if ((State.ring_queue[i] == 1) && (State.cumV_daily < P.VaccDosePerDay)) //for ring vaccination
				{
					DoVacc(State.vacc_queue[i], ts, 1);
					k++;
					l++;
				}
				else if ((State.ring_queue[i] == 0) && (State.cumVG_daily < P.VaccGeoDosePerDay)) //for geo vaccination
				{
					DoVacc(State.vacc_queue[i], ts, 0);
					k++;
					l++;
				}
			}
		}
	}

	//if ((State.nringvacc_queue - State.ringvacc_cum) <= P.VaccDosePerDay)
	//{
	//	maxvacc = State.nringvacc_queue;
	//}
	//else
	//{
	//	maxvacc = State.ringvacc_cum + P.VaccDosePerDay;
	//}


	//for (i = State.ringvacc_cum;i<maxvacc; i++)
	//{
	//
	//	DoVacc(State.ringvacc_queue[i],ts);
	//	k++;
	//}
	State.vacc_cum = k;
	State.vacc_ind = l;

	if (!P.ResetVaccQueue)
	{
		j = 0;
		//remove those vacced from the queue
		for (i = State.vacc_ind; i < State.nvacc_queue; i++)
		{
			State.vacc_queue[j] = State.vacc_queue[i];
			j++;
		}
		State.vacc_ind = 0;
		State.nvacc_queue = j;
	}
	else
	{
		State.vacc_ind = 0;
		State.nvacc_queue = 0;
	}

}

int TreatSweep(double t)
{
	int i, j, i2, j2, k, l, m, b, rd, prd, mrd, f, f1, f2, f3, f4, tn, bs, ad, ad2, ii, jj, npropvacc;
	int minx, maxx, miny, trig_thresh, nckwp, FirstPerson, LastPerson;
	double maxdist, maxdist2, dosepercase, dosepercell, logpop;
	unsigned short int ts, tstf, tstb, tsvb, tspf, tsmb, tsmf, tssdf, tskwpf;
	int global_trig;
	double r, rad2;

	ts = (unsigned short int) (P.TimeStepsPerDay * t);
	rd = (int)ceil(P.TreatRadius / P.mcwidth);
	prd = (int)ceil(P.PlaceCloseRadius / P.mcwidth);
	mrd = (int)ceil(P.MoveRestrRadius / P.mcwidth);
	f = f1 = 0;
	if (P.DoGlobalTriggers)
	{
		if (P.DoPerCapitaTriggers)
			global_trig = (int)floor(((double)State.trigDC) * P.GlobalIncThreshPop / ((double)P.N));
		else
			global_trig = State.trigDC;
	}
	else
		global_trig = 0;
	if ((P.DoPlaces) && (t >= P.TreatTimeStart) && (t < P.TreatTimeStart + P.TreatPlaceGeogDuration) && (State.cumT < P.TreatMaxCourses))
	{
		tstf = (unsigned short int) (P.TimeStepsPerDay * (t + P.TreatDelayMean + P.TreatProphCourseLength));
#pragma omp parallel for private(i,j,k,l,m,f) reduction(+:f1) schedule(static,1)
		for (i = 0; i < P.NumThreads; i++)
		{
			for (j = 0; j < P.PlaceTypeNum; j++)
			{
				for (k = 0; k < StateT[i].np_queue[j]; k++)
				{
					l = StateT[i].p_queue[j][k];
					if (P.DoPlaceGroupTreat)
					{
						f = StateT[i].pg_queue[j][k];
						for (m = ((int)Places[j][l].group_start[f]); m < ((int)(Places[j][l].group_start[f] + Places[j][l].group_size[f])); m++)
						{
							/*							if((Places[j][l].members[m]<0)||(Places[j][l].members[m]>P.N-1))
															fprintf(stderr,"\n*** npq=%i gn=%i h=%i m=%i j=%i l=%i f=%i s=%i n=%i ***\n",
																StateT[i].np_queue[j],
																Places[j][l].n,
																Places[j][l].members[m],
																m,j,l,f,
																(int) Places[j][l].group_start[f],
																(int) Places[j][l].group_size[f]);
														else
							*/							if ((!HOST_TO_BE_TREATED(Places[j][l].members[m])) && ((P.TreatPlaceTotalProp[j] == 1) || (ranf_mt(i) < P.TreatPlaceTotalProp[j])))
	DoProph(Places[j][l].members[m], ts, i);
						}
					}
					else
					{
						if ((Places[j][l].treat) && (!PLACE_TREATED(j, l)))
						{
							f1 = 1;
							Places[j][l].treat_end_time = tstf;
							for (m = 0; m < Places[j][l].n; m++)
								if (!HOST_TO_BE_TREATED(Places[j][l].members[m]))
								{
									if ((P.TreatPlaceTotalProp[j] == 1) || (ranf_mt(i) < P.TreatPlaceTotalProp[j]))
										DoProph(Places[j][l].members[m], ts, i);
								}
						}
						Places[j][l].treat = 0;
					}
				}
				StateT[i].np_queue[j] = 0;
			}
		}
	}
	if ((P.DoMassVacc) && (t >= P.VaccTimeStart))
		for (j = 0; j < 2; j++)
		{
			m = P.VaccMaxCourses;
			if (m > State.n_mvacc) m = State.n_mvacc;
#pragma omp parallel for private(i) schedule(static,1000)
			for (i = State.mvacc_cum; i < m; i++)
				DoVacc(State.mvacc_queue[i], ts, 0);
			State.mvacc_cum = m;
		}
	if ((t >= P.TreatTimeStart) || (t >= P.VaccTimeStart) || (t >= P.PlaceCloseTimeStart) || (t >= P.MoveRestrTimeStart) || (t >= P.SocDistTimeStart) || (t >= P.KeyWorkerProphTimeStart))
	{
		tstf = (unsigned short int) (P.TimeStepsPerDay * (t + P.TreatProphCourseLength) - 1);
		tstb = (unsigned short int) (P.TimeStepsPerDay * (t + P.TreatDelayMean));
		//tsvb= (unsigned short int) (P.TimeStepsPerDay * t)+1; //modified this so that the delay gets taken care of in the vaccine queue and ring and geo vaccinations align timewise
		tsvb = (unsigned short int) (P.TimeStepsPerDay * (t + P.VaccDelayMean));
		tspf = (unsigned short int) ceil(P.TimeStepsPerDay * (t + P.PlaceCloseDelayMean + P.PlaceCloseDuration));
		tsmf = (unsigned short int) ceil(P.TimeStepsPerDay * (t + P.MoveRestrDuration));
		tsmb = (unsigned short int) floor(P.TimeStepsPerDay * (t + P.MoveDelayMean));
		tssdf = (unsigned short int) ceil(P.TimeStepsPerDay * (t + P.SocDistDuration));
		tskwpf = (unsigned short int) ceil(P.TimeStepsPerDay * (t + P.KeyWorkerProphRenewalDuration));
		nckwp = (int)ceil(P.KeyWorkerProphDuration / P.TreatProphCourseLength);
#pragma omp parallel for private(tn,i,i2,j,j2,k,l,m,b,bs,minx,maxx,miny,f2,f3,f4,trig_thresh,r,ad,ad2,maxdist,maxdist2,dosepercase,logpop,ii,jj,npropvacc,rad2,FirstPerson,LastPerson) reduction(+:f) schedule(static,1) //,ii,ll,nvacc,maxvacc,maxdist,maxdist2,dosepercase,logpop
		for (tn = 0; tn < P.NumThreads; tn++)
			for (bs = tn; bs < P.NMCP; bs += P.NumThreads)
			{
				b = (int)(McellLookup[bs] - Mcells);
				ad = (P.DoAdUnits) ? AdUnits[Mcells[b].adunit].id : 0;
				if ((!P.RestrictTreatToTarget) || (Mcells[b].country == P.TargetCountry))
				{
					if ((Mcells[b].treat == 2) && (ts >= Mcells[b].treat_end_time))
					{
						f = 1;
						Mcells[b].treat = 0;
					}
					if ((Mcells[b].treat == 1) && (ts >= Mcells[b].treat_start_time))
					{
						f = 1;
						Mcells[b].treat = 2;
						Mcells[b].treat_trig = 0;
						Mcells[b].treat_end_time = tstf;
						for (i = 0; i < Mcells[b].n; i++)
						{
							l = Mcells[b].members[i];
							if ((!HOST_TO_BE_TREATED(l)) && ((P.TreatPropRadial == 1) || (ranf_mt(tn) < P.TreatPropRadial)))
								DoProphNoDelay(l, ts, tn, 1);
						}
					}

					trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.TreatCellIncThresh)) / P.IncThreshPop)) : P.TreatCellIncThresh;
					if ((t >= P.TreatTimeStart) && (Mcells[b].treat == 0) && (((Mcells[b].treat_trig >= trig_thresh) && (!P.DoGlobalTriggers)) || (global_trig >= P.TreatCellIncThresh)) && (P.TreatRadius2 > 0))
					{
						minx = (b / P.nmch); miny = (b % P.nmch);
						k = b;
						maxx = 0;
						i = j = m = f2 = 0;
						l = f3 = 1;
						if ((!P.TreatByAdminUnit) || (ad > 0))
						{
							ad2 = ad / P.TreatAdminUnitDivisor;
							do
							{
								if ((minx >= 0) && (minx < P.nmcw) && (miny >= 0) && (miny < P.nmch))
								{
									if (P.TreatByAdminUnit)
										f4 = (AdUnits[Mcells[k].adunit].id / P.TreatAdminUnitDivisor == ad2);
									else
										f4 = ((r = dist2_mm(Mcells + b, Mcells + k)) < P.TreatRadius2);
									if (f4)
									{
										f = f2 = 1;
										if ((Mcells[k].n > 0) && (Mcells[k].treat == 0) && ((!P.RestrictTreatToTarget) || (Mcells[k].country == P.TargetCountry)))
										{
											Mcells[k].treat_start_time = tstb;
											Mcells[k].treat = 1;
											maxx += Mcells[k].n;
										}
									}
								}
								if (j == 0)
									minx = minx + 1;
								else if (j == 1)
									miny = miny - 1;
								else if (j == 2)
									minx = minx - 1;
								else if (j == 3)
									miny = miny + 1;
								m = (m + 1) % l;
								if (m == 0)
								{
									j = (j + 1) % 4;
									i = (i + 1) % 2;
									if (i == 0) l++;
									if (j == 1) { f3 = f2; f2 = 0; }
								}
								k = ((minx + P.nmcw) % P.nmcw) * P.nmch + (miny + P.nmch) % P.nmch;
							} while ((f3) && (maxx < P.TreatMaxCoursesPerCase));
						}
					}
					//if Mcells.vacc==2, it is being vaccinated as a trigger cell, if Mcells.vacc==1, it is being vaccinated as a peripheral cell, but can still trigger vaccination, if Mcells.vacc==0, it
					//is currently not scheduled to be vaccinated

					if ((P.LimitGeoVaccDosesPerCase) && (P.StopVaccinationPostThreshold))
					{
						npropvacc = Mcells[b].n * (0.9 * P.VaccProp);
						//calculate proportion vaccinated in Mcell
						if ((Mcells[b].vacc == 2) && (ts >= Mcells[b].vacc_start_time) && (Mcells[b].popvacc < npropvacc))
						{
							f = 1;
							Mcells[b].vacc = 0;
							Mcells[b].vacc_trig = 0;

						}
					}
					else
					{
						if ((Mcells[b].vacc == 2) && (ts >= Mcells[b].vacc_start_time))
						{
							f = 1;
							Mcells[b].vacc = 0;
							Mcells[b].vacc_trig = 0;

						}
					}
					trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.VaccCellIncThresh)) / P.IncThreshPop)) : P.VaccCellIncThresh;
					if ((P.LimitGeoVaccDosesPerCase) && (P.DoGeoVaccination))
					{
						if ((!P.DoMassVacc) && (P.VaccRadius2 > 0) && (t >= P.VaccTimeStart) && (Mcells[b].vacc != 2) && (((Mcells[b].vacc_trig >= trig_thresh) && (!P.DoGlobalTriggers) && (State.cumV < P.VaccMaxCourses)) || (global_trig >= P.VaccCellIncThresh))) //changed from VaccTimeStart to VaccTimeStarGeo
						{
							minx = (b / P.nmch); miny = (b % P.nmch);
							k = b;
							maxx = 0;
							//number of vaccine doses to distribute geographically in this cell
							//nvacc=maxvacc=P.VaccDosesPerCasePerCell * Mcells[b].vacc_trig;
							//max distance between cells;
							maxdist = maxdist2 = dosepercase = dosepercell = 0;
							i = j = m = f2 = 0;
							l = f3 = 1;
							if ((!P.VaccByAdminUnit) || (ad > 0))
							{
								ad2 = ad / P.VaccAdminUnitDivisor;
								do
								{
									if ((minx >= 0) && (minx < P.nmcw) && (miny >= 0) && (miny < P.nmch))
									{
										if (P.VaccByAdminUnit)
										{
											f4 = (AdUnits[Mcells[k].adunit].id / P.VaccAdminUnitDivisor == ad2);
											r = 1e20;
										}
										else
											f4 = ((r = dist2_mm(Mcells + b, Mcells + k)) < P.VaccRadius2);
										if (f4)
										{
											f = f2 = 1;
											if (r < P.VaccMinRadius2)
												Mcells[k].vacc = 3;
											else if ((Mcells[k].n > 0) && ((!P.RestrictTreatToTarget) || (Mcells[k].country == P.TargetCountry))) //(nvacc>0)&&// (Mcells[k].vacc == 0) &&
											{
												Mcells[k].vacc_start_time = tsvb;
												//set vaccination status - if this is the trigger cell, set it equal to 2
												if (k == b)
													Mcells[k].vacc = 2;
												else if (Mcells[k].vacc != 2)
													Mcells[k].vacc = 1;

												//sample without replacement
												SampleWithoutReplacement(tn, Mcells[k].nh, Mcells[k].nh);

												//add to queue
												//for (ii = 0; ii < Mcells[k].n && maxx <= (P.VaccDosesPerCasePerCell * Mcells[b].vacc_trig); ii++)
												for (ii = 0; ii < Mcells[k].nh && maxx <= (P.VaccDosesPerCasePerCell * Mcells[b].vacc_trig); ii++)
												{
													FirstPerson = Households[Mcells[k].FirstHousehold + SamplingQueue[tn][ii]].FirstPerson;
													LastPerson = FirstPerson + Households[Mcells[k].FirstHousehold + SamplingQueue[tn][ii]].nh;
													for (jj = FirstPerson; jj < LastPerson; jj++)
													{
														if (((Hosts[jj].vacc_accept < P.VaccProp)) && !(HOST_TO_BE_VACCED(jj)) && (State.cumVG < P.VaccMaxCourses) && !(Hosts[jj].inf >= 5) && (HOST_AGE_YEAR(jj) >= P.MinVaccAge))
														{
															StateT[tn].vacc_queue[StateT[tn].nvacc_queue] = jj;
															StateT[tn].ring_queue[StateT[tn].nvacc_queue] = 0; //set this to zero to indicate geographic vaccination
															StateT[tn].nvacc_queue++;
															maxx++;
														}
													}
												}

												if (r > maxdist2)
												{
													maxdist2 = r;
												}


											}
										}

									}

									if (j == 0)
										minx = minx + 1;
									else if (j == 1)
										miny = miny - 1;
									else if (j == 2)
										minx = minx - 1;
									else if (j == 3)
										miny = miny + 1;
									m = (m + 1) % l;
									if (m == 0)
									{
										j = (j + 1) % 4;
										i = (i + 1) % 2;
										if (i == 0) l++;
										if (j == 1) { f3 = f2; f2 = 0; }
									}
									k = ((minx + P.nmcw) % P.nmcw) * P.nmch + (miny + P.nmch) % P.nmch;
									if ((k < 0) || (k > P.nmcw * P.nmch))
									{
										fprintf(stderr, "### %i %i %i\n", k, minx, miny);
									}
								} while ((f3) && (maxx <= P.VaccDosesPerCasePerCell * Mcells[b].vacc_trig));


							}

							dosepercase = ((Mcells[b].vacc_trig > 0) ? ((double)maxx / (double)Mcells[b].vacc_trig) : 0);
							dosepercell = maxx;
							maxdist = sqrt(maxdist2);
							StateT[tn].vaccdose_dist[VACCDOSE_GROUP((int)dosepercase)]++;
							StateT[tn].vaccdosecell_dist[VACCDOSECELL_GROUP((int)dosepercell)]++;
							StateT[tn].vaccdistance_dist[VACCDIST_GROUP((int)maxdist)]++;
							logpop = log((double)Mcells[b].n);
							StateT[tn].vaccpop_dist[POPDIST_GROUP((int)logpop)]++;

							Mcells[b].ntriggervacc++;
							Mcells[b].totalvacc += maxx;
							if (maxdist < Mcells[b].minvaccdist)
							{
								Mcells[b].minvaccdist = maxdist;
								Mcells[b].minvaccdist_dose = maxx;
								Mcells[b].minvaccdist_t = t;
							}
							if (maxdist > Mcells[b].maxvaccdist)
							{
								Mcells[b].maxvaccdist = maxdist;
								Mcells[b].maxvaccdist_dose = maxx;
								Mcells[b].maxvaccdist_t = t;
							}


						}
					}
					else
					{
						if ((!P.DoMassVacc) && (P.VaccRadius2 > 0) && (t >= P.VaccTimeStart) && (Mcells[b].vacc == 0) && (((Mcells[b].vacc_trig >= trig_thresh) && (!P.DoGlobalTriggers) && (State.cumV < P.VaccMaxCourses)) || (global_trig >= P.VaccCellIncThresh))) //changed from VaccTimeStart to VaccTimeStarGeo
						{
							minx = (b / P.nmch); miny = (b % P.nmch);
							k = b;
							//number of vaccine doses to distribute geographically in this cell
							//nvacc=maxvacc=P.VaccDosesPerCasePerCell * Mcells[b].vacc_trig;
							//max distance between cells;
							maxx = maxdist = maxdist2 = dosepercase = dosepercell = 0;
							if (Mcells[b].n > P.PopHighDensityCell)
							{
								rad2 = P.VaccRadiusHighDensity2;
							}
							else
							{
								rad2 = P.VaccRadius2;
							}
							i = j = m = f2 = 0;
							l = f3 = 1;
							if ((!P.VaccByAdminUnit) || (ad > 0))
							{
								ad2 = ad / P.VaccAdminUnitDivisor;
								do
								{
									if ((minx >= 0) && (minx < P.nmcw) && (miny >= 0) && (miny < P.nmch))
									{
										if (P.VaccByAdminUnit)
										{
											f4 = (AdUnits[Mcells[k].adunit].id / P.VaccAdminUnitDivisor == ad2);
											r = 1e20;
										}
										else
											f4 = ((r = dist2_mm(Mcells + b, Mcells + k)) < rad2);
										if (f4)
										{
											f = f2 = 1;
											if (r < P.VaccMinRadius2)
												Mcells[k].vacc = 3;
											else if ((Mcells[k].n > 0) && (Mcells[k].vacc == 0) && ((!P.RestrictTreatToTarget) || (Mcells[k].country == P.TargetCountry))) //(nvacc>0)&&// (Mcells[k].vacc == 0) &&
											{
												Mcells[k].vacc_start_time = tsvb;
												//set vaccination status - if this is the trigger cell, set it equal to 2
												Mcells[k].vacc = 1;

												//add to queue
												for (ii = 0; ii < Mcells[k].n; ii++)
												{
													jj = Mcells[k].members[ii];
													if (((Hosts[jj].vacc_accept < P.VaccProp)) && !(HOST_TO_BE_VACCED(jj)) && (State.cumV < P.VaccMaxCourses) && !(Hosts[jj].inf >= 5) && (HOST_AGE_YEAR(jj) >= P.MinVaccAge))
													{
														StateT[tn].vacc_queue[StateT[tn].nvacc_queue] = jj;
														StateT[tn].nvacc_queue++;
														maxx++;
													}
												}

												if (r > maxdist2)
												{
													maxdist2 = r;
												}
											}
										}

									}

									if (j == 0)
										minx = minx + 1;
									else if (j == 1)
										miny = miny - 1;
									else if (j == 2)
										minx = minx - 1;
									else if (j == 3)
										miny = miny + 1;
									m = (m + 1) % l;
									if (m == 0)
									{
										j = (j + 1) % 4;
										i = (i + 1) % 2;
										if (i == 0) l++;
										if (j == 1) { f3 = f2; f2 = 0; }
									}
									k = ((minx + P.nmcw) % P.nmcw) * P.nmch + (miny + P.nmch) % P.nmch;
									if ((k < 0) || (k > P.nmcw * P.nmch))
									{
										fprintf(stderr, "### %i %i %i\n", k, minx, miny);
									}
								} while (f3);
							}

							dosepercase = ((Mcells[b].vacc_trig > 0) ? ((double)maxx / (double)Mcells[b].vacc_trig) : 0);
							dosepercell = maxx;
							maxdist = sqrt(maxdist2);
							StateT[tn].vaccdose_dist[VACCDOSE_GROUP((int)dosepercase)]++;
							StateT[tn].vaccdosecell_dist[VACCDOSECELL_GROUP((int)dosepercell)]++;
							StateT[tn].vaccdistance_dist[VACCDIST_GROUP((int)maxdist)]++;
							logpop = log((double)Mcells[b].n);
							StateT[tn].vaccpop_dist[POPDIST_GROUP((int)logpop)]++;

							Mcells[b].ntriggervacc++;
							Mcells[b].totalvacc += maxx;
							if (maxdist < Mcells[b].minvaccdist)
							{
								Mcells[b].minvaccdist = maxdist;
								Mcells[b].minvaccdist_dose = maxx;
								Mcells[b].minvaccdist_t = t;
							}
							if (maxdist > Mcells[b].maxvaccdist)
							{
								Mcells[b].maxvaccdist = maxdist;
								Mcells[b].maxvaccdist_dose = maxx;
								Mcells[b].maxvaccdist_t = t;
							}
						}
					}
					if ((Mcells[b].placeclose == 2) && (ts >= Mcells[b].place_end_time))
					{
						f = 1;
						Mcells[b].placeclose = P.DoPlaceCloseOnceOnly;
					}
					trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.PlaceCloseCellIncThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncThresh;
					if ((P.DoPlaces) && (t >= P.PlaceCloseTimeStart) && (Mcells[b].placeclose == 0))
					{
						if (((P.PlaceCloseByAdminUnit) && (AdUnits[Mcells[b].adunit].place_close_trig < USHRT_MAX - 1) && (((double)AdUnits[Mcells[b].adunit].place_close_trig) / ((double)AdUnits[Mcells[b].adunit].NP) > P.PlaceCloseAdunitPropThresh))
							|| ((!P.PlaceCloseByAdminUnit) && (((Mcells[b].place_trig >= trig_thresh) && (!P.DoGlobalTriggers)) || (global_trig >= P.PlaceCloseCellIncThresh))))
						{
							//							if(P.PlaceCloseByAdminUnit) AdUnits[Mcells[b].adunit].place_close_trig=USHRT_MAX-1; // This means schools only close once
							AdUnits[Mcells[b].adunit].place_close_trig = 0; //Schools can close multiple times
							minx = (b / P.nmch); miny = (b % P.nmch);
							k = b;
							i = j = m = f2 = 0;
							l = f3 = 1;
							if ((!P.PlaceCloseByAdminUnit) || (ad > 0))
							{
								ad2 = ad / P.PlaceCloseAdminUnitDivisor;
								do
								{
									if ((minx >= 0) && (minx < P.nmcw) && (miny >= 0) && (miny < P.nmch))
									{
										if (P.PlaceCloseByAdminUnit)
											f4 = (AdUnits[Mcells[k].adunit].id / P.PlaceCloseAdminUnitDivisor == ad2);
										else
											f4 = ((r = dist2_mm(Mcells + b, Mcells + k)) < P.PlaceCloseRadius2);
										if (f4)
										{
											f2 = f = 1;
											if ((Mcells[k].n > 0) && (Mcells[k].placeclose == 0) && ((!P.RestrictTreatToTarget) || (Mcells[k].country == P.TargetCountry)))
											{
												Mcells[k].place_end_time = tspf;
												Mcells[k].place_trig = 0;
												Mcells[k].placeclose = 2;
												for (j2 = 0; j2 < P.PlaceTypeNum; j2++)
													if (j2 != HOTEL_PLACE_TYPE)
														for (i2 = 0; i2 < Mcells[k].np[j2]; i2++)
															DoPlaceClose(j2, Mcells[k].places[j2][i2], ts, tn, 1);
											}
										}
									}
									if (j == 0)
										minx = minx + 1;
									else if (j == 1)
										miny = miny - 1;
									else if (j == 2)
										minx = minx - 1;
									else if (j == 3)
										miny = miny + 1;
									m = (m + 1) % l;
									if (m == 0)
									{
										j = (j + 1) % 4;
										i = (i + 1) % 2;
										if (i == 0) l++;
										if (j == 1) { f3 = f2; f2 = 0; }
									}
									k = ((minx + P.nmcw) % P.nmcw) * P.nmch + (miny + P.nmch) % P.nmch;
									if ((k < 0) || (k > P.nmcw * P.nmch)) fprintf(stderr, "### %i %i %i\n ", k, minx, miny);
								} while (f3);
							}
						}
					}
					if ((Mcells[b].moverest == 2) && (ts >= Mcells[b].move_end_time))
					{
						f = 1;
						Mcells[b].moverest = 0;
					}
					if ((Mcells[b].moverest == 1) && (ts >= Mcells[b].move_start_time))
					{
						f = 1;
						Mcells[b].moverest = 2;
						Mcells[b].move_trig = 0;
						Mcells[b].move_end_time = tsmf;
					}
					trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.MoveRestrCellIncThresh)) / P.IncThreshPop)) : P.MoveRestrCellIncThresh;
					if ((t >= P.MoveRestrTimeStart) && (Mcells[b].moverest == 0) && (((Mcells[b].move_trig >= trig_thresh) && (!P.DoGlobalTriggers)) || (global_trig >= P.MoveRestrCellIncThresh)))
					{
						minx = (b / P.nmch); miny = (b % P.nmch);
						k = b;
						i = j = m = f2 = 0;
						l = f3 = 1;
						if ((!P.MoveRestrByAdminUnit) || (ad > 0))
						{
							ad2 = ad / P.MoveRestrAdminUnitDivisor;
							do
							{
								if ((minx >= 0) && (minx < P.nmcw) && (miny >= 0) && (miny < P.nmch))
								{
									if (P.MoveRestrByAdminUnit)
										f4 = (AdUnits[Mcells[k].adunit].id / P.MoveRestrAdminUnitDivisor == ad2);
									else
										f4 = ((r = dist2_mm(Mcells + b, Mcells + k)) < P.MoveRestrRadius2);
									if (f4)
									{
										f = f2 = 1;
										if ((Mcells[k].n > 0) && (Mcells[k].moverest == 0) && ((!P.RestrictTreatToTarget) || (Mcells[k].country == P.TargetCountry)))
										{
											Mcells[k].move_start_time = tsmb;
											Mcells[k].moverest = 1;
										}
									}
								}
								if (j == 0)
									minx = minx + 1;
								else if (j == 1)
									miny = miny - 1;
								else if (j == 2)
									minx = minx - 1;
								else if (j == 3)
									miny = miny + 1;
								m = (m + 1) % l;
								if (m == 0)
								{
									j = (j + 1) % 4;
									i = (i + 1) % 2;
									if (i == 0) l++;
									if (j == 1) { f3 = f2; f2 = 0; }
								}
								k = ((minx + P.nmcw) % P.nmcw) * P.nmch + (miny + P.nmch) % P.nmch;
							} while (f3);
						}
					}

					if ((Mcells[b].socdist == 2) && (ts >= Mcells[b].socdist_end_time))
					{
						f = 1;
						Mcells[b].socdist = 0;
					}
					trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.SocDistCellIncThresh)) / P.IncThreshPop)) : P.SocDistCellIncThresh;
					if ((t >= P.SocDistTimeStart) && (Mcells[b].socdist == 0) && (((Mcells[b].socdist_trig >= trig_thresh) && (!P.DoGlobalTriggers)) || (global_trig >= P.SocDistCellIncThresh)))
					{
						minx = (b / P.nmch); miny = (b % P.nmch);
						k = b;
						i = j = m = f2 = 0;
						l = f3 = 1;
						do
						{
							if ((minx >= 0) && (minx < P.nmcw) && (miny >= 0) && (miny < P.nmch))
								if (dist2_mm(Mcells + b, Mcells + k) < P.SocDistRadius2)
								{
									f = f2 = 1;
									if ((Mcells[k].n > 0) && (Mcells[k].socdist == 0) && ((!P.RestrictTreatToTarget) || (Mcells[k].country == P.TargetCountry)))
									{
										Mcells[k].socdist = 2;
										Mcells[k].socdist_trig = 0;
										Mcells[k].socdist_end_time = tssdf;
									}
								}
							if (j == 0)
								minx = minx + 1;
							else if (j == 1)
								miny = miny - 1;
							else if (j == 2)
								minx = minx - 1;
							else if (j == 3)
								miny = miny + 1;
							m = (m + 1) % l;
							if (m == 0)
							{
								j = (j + 1) % 4;
								i = (i + 1) % 2;
								if (i == 0) l++;
								if (j == 1) { f3 = f2; f2 = 0; }
							}
							k = ((minx + P.nmcw) % P.nmcw) * P.nmch + (miny + P.nmch) % P.nmch;
						} while (f3);
					}

					if ((Mcells[b].keyworkerproph == 2) && (ts >= Mcells[b].keyworkerproph_end_time))
					{
						f = 1;
						Mcells[b].keyworkerproph = 0;
					}
					trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.KeyWorkerProphCellIncThresh)) / P.IncThreshPop)) : P.KeyWorkerProphCellIncThresh;
					if ((P.DoPlaces) && (t >= P.KeyWorkerProphTimeStart) && (Mcells[b].keyworkerproph == 0) && (((Mcells[b].keyworkerproph_trig >= trig_thresh) && (!P.DoGlobalTriggers)) || (global_trig >= P.KeyWorkerProphCellIncThresh)))
					{
						minx = (b / P.nmch); miny = (b % P.nmch);
						k = b;
						i = j = m = f2 = 0;
						l = f3 = 1;
						do
						{
							if ((minx >= 0) && (minx < P.nmcw) && (miny >= 0) && (miny < P.nmch))
								if (dist2_mm(Mcells + b, Mcells + k) < P.KeyWorkerProphRadius2)
								{
									f = f2 = 1;
									if ((Mcells[k].n > 0) && (Mcells[k].keyworkerproph == 0) && ((!P.RestrictTreatToTarget) || (Mcells[k].country == P.TargetCountry)))
									{
										Mcells[k].keyworkerproph = 2;
										Mcells[k].keyworkerproph_trig = 0;
										Mcells[k].keyworkerproph_end_time = tskwpf;
										for (i2 = 0; i2 < Mcells[k].n; i2++)
										{
											j2 = Mcells[k].members[i2];
											if ((Hosts[j2].keyworker) && (!HOST_TO_BE_TREATED(j2)))
												DoProphNoDelay(j2, ts, tn, nckwp);
										}
									}
								}
							if (j == 0)
								minx = minx + 1;
							else if (j == 1)
								miny = miny - 1;
							else if (j == 2)
								minx = minx - 1;
							else if (j == 3)
								miny = miny + 1;
							m = (m + 1) % l;
							if (m == 0)
							{
								j = (j + 1) % 4;
								i = (i + 1) % 2;
								if (i == 0) l++;
								if (j == 1) { f3 = f2; f2 = 0; }
							}
							k = ((minx + P.nmcw) % P.nmcw) * P.nmch + (miny + P.nmch) % P.nmch;
						} while (f3);
					}

				}
			}
		for (i = 0; i < P.NumThreads; i++)
		{
			State.cumT += StateT[i].cumT;
			State.cumTP += StateT[i].cumTP;
			State.cumUT += StateT[i].cumUT;
			//State.cumV+=StateT[i].cumV;
			StateT[i].cumT = StateT[i].cumUT = StateT[i].cumTP = 0; //StateT[i].cumV=0;
		}
	}
	f += f1;


	//if (t >= P.GeoVaccTimeStart)
	//{
	//	//added this - process geographical vaccine queue - 26/11/19
	//	//combine queues for individual threads into one master list
	//	for (i = 0; i < P.NumThreads; i++)
	//	{
	//		for (j = StateT[i].geovacc_cum; j < StateT[i].ngeovacc_queue; j++)
	//		{
	//			State.geovacc_queue[State.ngeovacc_queue] = StateT[i].geovacc_queue[j];
	//			State.ngeovacc_queue++;
	//		}
	//		StateT[i].geovacc_cum = StateT[i].ngeovacc_queue;
	//	}

	//	//process the ring vaccination queue
	//	k = State.geovacc_cum;
	//	l = State.geovacc_ind;

	//	for (i = State.geovacc_ind; i < State.ngeovacc_queue; i++)
	//	{
	//		if (((P.VaccDosePerDay >= 0) ? ((State.cumV_daily) < P.VaccDosePerDay) : 1) && ((State.cumV) < P.VaccMaxCourses))
	//		{
	//			if (!HOST_TO_BE_VACCED(State.geovacc_queue[i]))
	//			{

	//				DoVacc(State.geovacc_queue[i], Mcells[Hosts[State.geovacc_queue[i]].mcell].vacc_start_time);
	//				k++;
	//			}
	//			l++;
	//		}
	//	}
	//	State.geovacc_cum = k;
	//	State.geovacc_ind = l;
	//}
	for (j = 0; j < P.NumThreads; j++)
	{
		for (i = 0; i < NUM_VACCDIST_GROUPS; i++)
		{
			vaccdistance_dist[i] = vaccdistance_dist[i] + StateT[j].vaccdistance_dist[i];
			StateT[j].vaccdistance_dist[i] = 0;
		}
		for (i = 0; i < NUM_VACCDOSE_GROUPS; i++)
		{
			vaccdose_dist[i] = vaccdose_dist[i] + StateT[j].vaccdose_dist[i];
			StateT[j].vaccdose_dist[i] = 0;
		}
		for (i = 0; i < NUM_VACCDOSE_GROUPS; i++)
		{
			vaccdosering_dist[i] = vaccdosering_dist[i] + StateT[j].vaccdosering_dist[i];
			StateT[j].vaccdosering_dist[i] = 0;
		}
		for (i = 0; i < NUM_VACCDOSECELL_GROUPS; i++)
		{
			vaccdosecell_dist[i] = vaccdosecell_dist[i] + StateT[j].vaccdosecell_dist[i];
			StateT[j].vaccdosecell_dist[i] = 0;
		}
		for (i = 0; i < NUM_POP_GROUPS; i++)
		{
			vaccpop_dist[i] = vaccpop_dist[i] + StateT[j].vaccpop_dist[i];
			StateT[j].vaccpop_dist[i] = 0;
		}
	}



	return (f > 0);
}
