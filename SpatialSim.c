/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "model_constants.h"
#include "model_structs.h"
#include "model_fns.h"
#include "model_globals.h"
#include "binio.h"



int main(int argc,char *argv[])
{
	char ParamFile[1024],DensityFile[1024],NetworkFile[1024],AirTravelFile[1024],SchoolFile[1024],RegDemogFile[1024],InterventionFile[MAXINTFILE][1024],PreParamFile[1024],buf[2048],*sep;
	int i,j,k,GotP,GotPP,GotO,GotD,GotL,GotS,GotA,GotScF,GotIF,Perr,cl;
	double t,t2;

	Perr=0;
	fprintf(stderr,"sizeof(int)=%zi sizeof(long)=%zi sizeof(float)=%zi sizeof(double)=%zi sizeof(unsigned short int)=%zi sizeof(int *)=%zi\n",sizeof(int),sizeof(long),sizeof(float),sizeof(double),sizeof(unsigned short int),sizeof(int *));
	cl=clock();
	if(argc<7)
		Perr=1;
	else
		{
		i=argc-4;
		sscanf(argv[i],"%li",&P.seed1);
		sscanf(argv[i+1],"%li",&P.seed2);
		sscanf(argv[i+2],"%li",&P.seed3);
		sscanf(argv[i+3],"%li",&P.seed4);
		GotP=GotO=GotD=GotL=GotS=GotA=GotScF=GotIF=GotPP=0;
		P.PlaceCloseIndepThresh=P.LoadSaveNetwork=P.DoHeteroDensity=P.DoPeriodicBoundaries=P.DoSchoolFile=P.DoAdunitDemog=P.OutputDensFile=P.MaxNumThreads=P.DoInterventionFile=0;
		P.PreControlClusterIdCaseThreshold=0;
		P.R0scale=1.0;
		P.KernelOffsetScale=P.KernelPowerScale=1.0; //added this so that kernel parameters are only changed if input from the command line: ggilani - 15/10/2014
		P.CT_scale1=P.CT_scale2=P.CT_thresh1=P.CT_thresh2=P.CTinc_scale1=P.CTinc_scale2=1; //added this so that contact tracing capacity is only scaled if input on the command line
		P.BC_scale=1.0; //added this to scale cross border contact from the command line: ggilani 17/12/14
		P.RR1=P.RR2=P.RR3=1.0; //added this to scale reporting rates from command line: ggilani 03/02/15
		P.RRAlt=1.0; //added this to scale reporting rate for Conakry from the command line: ggilani 18/06/2015
		P.AltAU1=P.AltAU2=0; //this is going to be equal to the adunit code which gets the separate capital city reporting rate - if we don't read RRAlt in from the command line it will remain 0 and won't have any effect
		P.CapitalCityAddEffect=0.0; //initially have the increased additive connectivity between the capital and other adunits set to zero, i.e. no increase
		P.RelativeSusceptibilityGuinea=1.0; //relative susceptibility in Guinea - used to scale now R0 in Guinea.
		P.RelativeSusceptibilityLiberia=1.0; //relative susceptibility in Liberia - used to scale R0 in Liberia
		P.VaccRingScale=1;
		P.VaccCaseScale=1;
		P.VaccPropScale=1.0; //added scaling factors for ring vaccination so that they can be set from the command line.
		P.VaccEffTimeScale=1.0;
		P.VaccDelayScale=1.0;
		P.DoSaveSnapshot=P.DoLoadSnapshot=P.DoRadiationMobility=P.DoRoadNetwork=0; //also set radiation mobility and road network to zero initially: ggilani 09/02/15
		for(i=1;i<argc-4;i++)
			{
			if((argv[i][0]!='/') && ((argv[i][2]!=':')&&(argv[i][3]!=':'))) Perr=1;
			if(argv[i][1]=='P' && argv[i][2]==':')
				{
				GotP=1;
				sscanf(&argv[i][3],"%s",ParamFile);
				}
			else if(argv[i][1]=='O' && argv[i][2]==':')
				{
				GotO=1;
				sscanf(&argv[i][3],"%s",OutFileBase);
				}
			else if(argv[i][1]=='D' && argv[i][2]==':')
				{
				GotD=1;
				sscanf(&argv[i][3],"%s",DensityFile);
				P.DoHeteroDensity=1;
				P.DoPeriodicBoundaries=0;
				}
			else if(argv[i][1]=='L' && argv[i][2]==':')
				{
				GotL=1;
				P.LoadSaveNetwork=1;
				sscanf(&argv[i][3],"%s",NetworkFile);
				}
			else if(argv[i][1]=='S' && argv[i][2]==':')
				{
				P.LoadSaveNetwork=2;
				GotS=1;
				sscanf(&argv[i][3],"%s",NetworkFile);
				}
			else if(argv[i][1]=='R' && argv[i][2]==':')
				{
				sscanf(&argv[i][3],"%lg",&P.R0scale);
				}
			else if(argv[i][1]=='K' && argv[i][2]=='P' && argv[i][3]==':') //added Kernel Power and Offset scaling so that it can easily be altered from the command line in order to vary the kernel quickly: ggilani - 15/10/14
				{
					sscanf(&argv[i][4],"%lg",&P.KernelPowerScale);
				}
			else if(argv[i][1]=='K' && argv[i][2]=='O' && argv[i][3]==':')
				{
					sscanf(&argv[i][4],"%lg",&P.KernelOffsetScale);
				}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '1' && argv[i][5] == ':') // generic command line specified param - matched to #1 in param file
			{
				sscanf(&argv[i][6], "%lf", &P.clP1);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '2' && argv[i][5] == ':') // generic command line specified param - matched to #2 in param file
			{
				sscanf(&argv[i][6], "%lf", &P.clP2);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '3' && argv[i][5] == ':') // generic command line specified param - matched to #3 in param file
			{
				sscanf(&argv[i][6], "%lf", &P.clP3);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '4' && argv[i][5] == ':') // generic command line specified param - matched to #4 in param file
			{
				sscanf(&argv[i][6], "%lf", &P.clP4);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '5' && argv[i][5] == ':') // generic command line specified param - matched to #5 in param file
			{
				sscanf(&argv[i][6], "%lf", &P.clP5);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '6' && argv[i][5] == ':') // generic command line specified param - matched to #6 in param file
			{
				sscanf(&argv[i][6], "%lf", &P.clP6);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '7' && argv[i][5] == ':') // generic command line specified param - matched to #7 in param file
			{
				sscanf(&argv[i][6], "%lf", &P.clP7);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '8' && argv[i][5] == ':') // generic command line specified param - matched to #8 in param file
			{
				sscanf(&argv[i][6], "%lf", &P.clP8);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '9' && argv[i][5] == ':') // generic command line specified param - matched to #8 in param file
			{
				sscanf(&argv[i][6], "%lf", &P.clP9);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '1' && argv[i][5] == '0' && argv[i][6] == ':') // generic command line specified param - matched to #8 in param file
			{
				sscanf(&argv[i][7], "%lf", &P.clP10);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '1' && argv[i][5] == '1' && argv[i][6] == ':') // generic command line specified param - matched to #8 in param file
			{
				sscanf(&argv[i][7], "%lf", &P.clP11);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '1' && argv[i][5] == '2' && argv[i][6] == ':') // generic command line specified param - matched to #8 in param file
			{
				sscanf(&argv[i][7], "%lf", &P.clP12);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '1' && argv[i][5] == '3' && argv[i][6] == ':') // generic command line specified param - matched to #8 in param file
			{
				sscanf(&argv[i][7], "%lf", &P.clP13);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '1' && argv[i][5] == '4' && argv[i][6] == ':') // generic command line specified param - matched to #8 in param file
			{
				sscanf(&argv[i][7], "%lf", &P.clP14);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '1' && argv[i][5] == '5' && argv[i][6] == ':') // generic command line specified param - matched to #8 in param file
			{
				sscanf(&argv[i][7], "%lf", &P.clP15);
			}


			else if(argv[i][1]=='C' && argv[i][2]=='C' && argv[i][3]=='1' && argv[i][4]==':') //added contact tracing capacity scaling
				{
					sscanf(&argv[i][5],"%i",&P.CT_scale1);
			}
			else if(argv[i][1]=='C' && argv[i][2]=='C' && argv[i][3]=='2' && argv[i][4]==':') //added contact tracing capacity scaling
			{
				sscanf(&argv[i][5],"%i",&P.CT_scale2);
			}
			else if(argv[i][1]=='C' && argv[i][2]=='I' && argv[i][3]=='1' && argv[i][4]==':') //added contact tracing capacity scaling
				{
					sscanf(&argv[i][5],"%i",&P.CTinc_scale1);
			}
			else if(argv[i][1]=='C' && argv[i][2]=='I' && argv[i][3]=='2' && argv[i][4]==':') //added contact tracing capacity scaling
			{
				sscanf(&argv[i][5],"%i",&P.CTinc_scale2);
			}
			else if(argv[i][1]=='C' && argv[i][2]=='T' && argv[i][3]=='1' && argv[i][4]==':') //added contact tracing threshold scaling
				{
					sscanf(&argv[i][5],"%i",&P.CT_thresh1);
				}
			else if(argv[i][1]=='C' && argv[i][2]=='T' && argv[i][3]=='2' && argv[i][4]==':') //added contact tracing threshold scaling
				{
					sscanf(&argv[i][5],"%i",&P.CT_thresh2);
				}
			else if(argv[i][1]=='B' && argv[i][2]=='C' && argv[i][3]==':') //added border control scaling
			{
				sscanf(&argv[i][4],"%lg",&P.BC_scale);
			}
			else if(argv[i][1]=='V' && argv[i][2]=='R' && argv[i][3]==':') //added number of rings scaling
			{
				sscanf(&argv[i][4],"%i",&P.VaccRingScale);
			}
			else if(argv[i][1]=='V' && argv[i][2]=='C' && argv[i][3]==':') //added number of threshold cases
			{
				sscanf(&argv[i][4],"%i",&P.VaccCaseScale);
			}
			else if(argv[i][1]=='V' && argv[i][2]=='P' && argv[i][3]==':') //added vaccination proportion scaling
			{
				sscanf(&argv[i][4],"%lg",&P.VaccPropScale);
			}
			else if(argv[i][1]=='V' && argv[i][2]=='D' && argv[i][3]==':') //added vaccination delay scaling
			{
				sscanf(&argv[i][4],"%lg",&P.VaccDelayScale);
			}
			else if(argv[i][1]=='V' && argv[i][2]=='E' && argv[i][3]==':') //added vaccination time to efficacy scaling
			{
				sscanf(&argv[i][4],"%lg",&P.VaccEffTimeScale);
			}
			else if(argv[i][1]=='R' && argv[i][2]=='R' && argv[i][3]=='1' && argv[i][4]==':') //case detection rate 1, formerly case detection scaling for guinea: ggilani 15/02/22
			{
				sscanf(&argv[i][5],"%lg",&P.RR1);
			}
			else if(argv[i][1]=='R' && argv[i][2]=='R' && argv[i][3]=='2' && argv[i][4]==':') //case detection rate 2, formerly case detection scaling for liberia: ggilani 15/02/22
			{
				sscanf(&argv[i][5],"%lg",&P.RR2);
			}
			else if(argv[i][1]=='R' && argv[i][2]=='R' && argv[i][3]=='3' && argv[i][4]==':') //case detection rate 3, formerly case detection scaling for sierra leone: ggilani 15/02/22
			{
				sscanf(&argv[i][5],"%lg",&P.RR3);
			}
			else if(argv[i][1]=='R' && argv[i][2]=='R' && argv[i][3]=='A' && argv[i][4]==':') //added case detection scaling for Conakry specifically
			{
				sscanf(&argv[i][5],"%lg",&P.RRAlt);
			}
			else if(argv[i][1]=='R' && argv[i][2]=='R' && argv[i][3]=='A' && argv[i][4]=='1' && argv[i][5]==':') //added case detection scaling for Conakry specifically - specifies which adunit is conakry
			{
				sscanf(&argv[i][6],"%i",&P.AltAU1);
			}
			else if(argv[i][1]=='R' && argv[i][2]=='R' && argv[i][3]=='A' && argv[i][4]=='2' && argv[i][5]==':') //added case detection scaling for Conakry specifically - specifies which adunit is conakry
			{
				sscanf(&argv[i][6],"%i",&P.AltAU2);
			}
			else if(argv[i][1]=='I' && argv[i][2]=='C' && argv[i][3]==':') //Added connectivity between admin units and capital cities
			{
				sscanf(&argv[i][4],"%lg",&P.CapitalCityAddEffect);
			}
			else if(argv[i][1]=='R' && argv[i][2]=='S' && argv[i][3]=='G' && argv[i][4]==':') //Relative susceptibility in Guinea
			{
				sscanf(&argv[i][5],"%lg",&P.RelativeSusceptibilityGuinea);
			}
			else if(argv[i][1]=='R' && argv[i][2]=='S' && argv[i][3]=='L' && argv[i][4]==':') //Relative susceptibility in Liberia
			{
				sscanf(&argv[i][5],"%lg",&P.RelativeSusceptibilityLiberia);
			}
			else if(argv[i][1]=='A' && argv[i][2]==':')
				{
				GotA=1;
				sscanf(&argv[i][3],"%s",AirTravelFile);
				}
			else if(argv[i][1]=='s' && argv[i][2]==':')
				{
				GotScF=1;
				sscanf(&argv[i][3],"%s",SchoolFile);
				}
			else if(argv[i][1]=='T' && argv[i][2]==':')
				{
				sscanf(&argv[i][3],"%i",&P.PreControlClusterIdCaseThreshold);
				}
			else if(argv[i][1]=='C' && argv[i][2]==':')
				{
				sscanf(&argv[i][3],"%i",&P.PlaceCloseIndepThresh);
				}
			else if(argv[i][1]=='d' && argv[i][2]==':')
				{
				P.DoAdunitDemog=1;
				sscanf(&argv[i][3],"%s",RegDemogFile);
				}
			else if(argv[i][1]=='c' && argv[i][2]==':')
				{
				sscanf(&argv[i][3],"%i",&P.MaxNumThreads);
				}
			//STB New Distribution Network 
			else if(argv[i][1]=='M' && argv[i][2]==':')
				{
				P.OutputDensFile=1;
				sscanf(&argv[i][3],"%s",OutDensFile);
				}
			else if(argv[i][1]=='I' && argv[i][2]==':')
				{
				GotIF=1;
				sscanf(&argv[i][3],"%s",InterventionFile[P.DoInterventionFile]);
				P.DoInterventionFile++;
				}
			else if(argv[i][1]=='L' && argv[i][2]=='S' && argv[i][3]==':')
				{
				sscanf(&argv[i][4],"%s",SnapshotLoadFile);
				P.DoLoadSnapshot=1;
				}
			else if(argv[i][1]=='P' && argv[i][2]=='P' && argv[i][3]==':')
				{
				sscanf(&argv[i][4],"%s",PreParamFile);
				GotPP=1;
				}
			else if(argv[i][1]=='S' && argv[i][2]=='S' && argv[i][3]==':')
				{
				sscanf(&argv[i][4],"%s",buf);
				fprintf(stderr,"### %s\n",buf);
				sep=strchr(buf,',');
				if(!sep) 
					Perr=1;
				else
					{
					P.DoSaveSnapshot=1;
					*sep=' ';
					sscanf(buf,"%lg %s",&(P.SnapshotSaveTime),SnapshotSaveFile);
					}
				}
			else if(argv[i][1]=='N' && argv[i][2]==':')
				{
				fprintf(stderr,"This Version of SpatialSim doesn't have FRESSCA compiled into it\n");
				}
			//added this for radiation mobility model: ggilani 09/02/15
			else if(argv[i][1]=='r' && argv[i][2]==':')
				{
				P.DoRadiationMobility=1;
				sscanf(&argv[i][3],"%s",RadiationFile);
				}

			//added this for transport network file: ggilani 12/02/15
			else if(argv[i][1]=='R' && argv[i][2]=='N' && argv[i][3]==':')
				{
				P.DoRoadNetwork=1;
				sscanf(&argv[i][4],"%s",RoadNetworkFile);
				}
			}
		if(((GotS)&&(GotL))||(!GotP)||(!GotO)) Perr=1;
		}
	sprintf(OutFile,"%s",OutFileBase);
	fprintf(stderr,"Param=%s\nOut=%s\nDens=%s\n",ParamFile,OutFile,DensityFile);
	if(Perr) ERR_CRITICAL("Syntax:\nSpatialSim /P:ParamFile /O:OutputFile [/A:AirTravelFile] [/s:SchoolFile] [/D:DensityFile] [/L:NetworkFileToLoad | /S:NetworkFileToSave] [/R:R0scaling] Seed1 Seed2 Seed3 Seed4\n");
#ifdef DO_OMP_PARALLEL
	P.NumThreads=omp_get_max_threads();
	if((P.MaxNumThreads>0)&&(P.MaxNumThreads<P.NumThreads)) P.NumThreads=P.MaxNumThreads;
	if(P.NumThreads>MAX_NUM_THREADS)
		{
		fprintf(stderr,"Assigned number of threads > MAX_NUM_THREADS\n");
		omp_set_num_threads(MAX_NUM_THREADS);
		}
	else
		omp_set_num_threads(P.NumThreads);
#pragma omp parallel default(shared)
		{
		fprintf(stderr,"Thread %i initialised\n",omp_get_thread_num());
		}
	/* fprintf(stderr,"int=%i\tfloat=%i\tdouble=%i\tint *=%i\n",(int) sizeof(int),(int) sizeof(float),(int) sizeof(double),(int) sizeof(int *));
	*/
#else
	P.NumThreads=1;
#endif
	if(!GotPP)
		{
#ifdef UNIX
		sprintf(PreParamFile,"../Pre_%s",ParamFile);
#else
		sprintf(PreParamFile,"..\\Pre_%s",ParamFile);
#endif
		}
	ReadParams(ParamFile,PreParamFile);
	if(GotScF) P.DoSchoolFile=1;
	if(P.DoAirports)
		{
		if(!GotA) ERR_CRITICAL("Syntax:\nSpatialSim /P:ParamFile /O:OutputFile /A:AirTravelFile [/s:SchoolFile] [/D:DensityFile] [/L:NetworkFileToLoad | /S:NetworkFileToSave] [/R:R0scaling] Seed1 Seed2 Seed3 Seed4\n");
		ReadAirTravel(AirTravelFile);
		}
	
	//initialise random number counts to zero
	count_ranf=0;
	count_ranfmt=0;

	SetupModel(DensityFile,NetworkFile,SchoolFile,RegDemogFile);


//	signal(SIGABRT,HandleBreak);
//	signal(SIGINT,HandleBreak);
//	signal(SIGTERM,HandleBreak);
	for(i=0;i<MAX_ADUNITS;i++) AdUnits[i].NI=0;
	if(P.DoInterventionFile>0)
		for(i=0;i<P.DoInterventionFile;i++)
			ReadInterventions(InterventionFile[i]);

	fprintf(stderr,"Model setup in %lg seconds\n",((double) (clock()-cl))/CLOCKS_PER_SEC);

	//add some code to generate allow for an array of seeds to be generated which will be used for each run
	if ((P.ResetSeeds) && (P.DoFixedSeeds))
	{
		//error check for number of runs
		if (P.NR > MAX_FIXED_SEEDS)
		{
			ERR_CRITICAL("Number of runs is greater than the size of the seed array!\n");
		}

		for (i = 0; (i < P.NR); i++)
		{
			P.FixedSeeds[i][0]	= (int)(ranf() * 1e8);
			P.FixedSeeds[i][1]	= (int)(ranf() * 1e8);
		}
	}

	if(!P.ResetSeeds)
	{
		setall(P.seed3,P.seed4);
	}

	P.NRactE=P.NRactNE=0;
	for(i=0;(i<P.NR)&&(P.NRactNE<P.NRN)&&(!InterruptRun);i++)
		{
		if(P.NR>1)
			{
			sprintf(OutFile,"%s.%i",OutFileBase,i);
			fprintf(stderr,"Realisation %i   (time=%lg nr_ne=%i)\n",i+1,((double) (clock()-cl))/CLOCKS_PER_SEC,P.NRactNE);
			}
		if(P.ResetSeeds)
		{
			if (P.KeepSameSeeds) 
			{
				P.newseed1 = P.seed3;
				P.newseed2 = P.seed4;
			}
			else if (P.DoFixedSeeds)
			{
				P.newseed1 = P.FixedSeeds[i][0];
				P.newseed2 = P.FixedSeeds[i][1];
			}
			else
			{
				if (i == 0)
				{
					P.newseed1 = P.seed3;
					P.newseed2 = P.seed4;
				}
				else
				{
					//sample 4 new seeds to use as random number generators
					P.newseed1 = (int) (ranf()*1e8);
					P.newseed2 = (int) (ranf()*1e8);
				}
			}
			//save these seeds to file
			SaveRandomSeeds();
			//reset seeds
			setall(P.newseed1,P.newseed2);
			//fprintf(stderr, "%i, %i\n", P.newseed1,P.newseed2);
			//fprintf(stderr, "%f\n", ranf());
		}

		InitModel(i); //passing run number into RunModel so we can save run number in the infection event log: ggilani - 15/10/2014
		//fprintf(stderr, "%f\n", ranf());
		//fprintf(stderr, "%f\n", ranf_mt(5));
		if(P.DoLoadSnapshot) LoadSnapshot();
		RunModel(i); //passing run number into RunModel so we can save run number in the infection event log: ggilani - 15/10/2014
		if(((!TimeSeries[P.NumSamples-1].extinct)||(!P.OutputNonExtinct))&&(P.OutputAll)) 
			{
				SaveResults(); 
			}
		if ((P.DoRecordInfEvents) && (P.RecordInfEventsPerRun == 1))
		{
			SaveEvents();
		}
		}
	sprintf(OutFile,"%s",OutFileBase);
	
	//Calculate origin destination matrix if needed
	if((P.DoAdUnits)&&(P.DoOriginDestinationMatrix))
	{
		CalcOriginDestMatrix_adunit();
		SaveOriginDestMatrix();
	}

	P.NRactual=P.NRactNE;
	TSMean=TSMeanNE;TSVar=TSVarNE;
	if((P.DoRecordInfEvents)&&(P.RecordInfEventsPerRun==0))
	{
		SaveEvents();
	}
	if (P.DoSummaryOutput) //added flag here to choose whether to save summary results or not: ggilani 29/03/22
	{
		sprintf(OutFile, "%s.avNE", OutFileBase);
		SaveSummaryResults();
		P.NRactual = P.NRactE;
		TSMean = TSMeanE; TSVar = TSVarE;
		sprintf(OutFile, "%s.avE", OutFileBase);
		SaveSummaryResults();
	}

#ifdef WIN32_BM
	Gdiplus::GdiplusShutdown(m_gdiplusToken);
#endif
	fprintf(stderr,"Extinction in %i out of %i runs\n",P.NRactE,P.NRactNE+P.NRactE);
	fprintf(stderr,"Model ran in %lg seconds\n",((double) (clock()-cl))/CLOCKS_PER_SEC);
	fprintf(stderr,"Model finished\n");
}



void HandleBreak(int signum)
{
InterruptRun++;
fprintf(stderr,"Please wait whilst current run finishes...\n");
if(InterruptRun==2)
	{
	signal(SIGABRT,SIG_DFL);
	signal(SIGINT,SIG_DFL);
	signal(SIGTERM,SIG_DFL);
	}
}


