/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "model_constants.h"
#include "model_structs.h"
#include "model_fns.h"
#include "model_globals_ext.h"
#include "binio.h"


void SetupModel(char *DensityFile,char *NetworkFile,char *SchoolFile, char *RegDemogFile)
{
	int i,j,k,l,m,i2,j2,l2,m2,tn,BedCapacity; //added tn as variable for multi-threaded loops: 28/11/14
	unsigned int rn;
	double t,s,s2,s3,x,y,t2,t3,d,q,contact_scale,inf_period,ProbHosp,probMort;
	char buf[2048];
	FILE *dat,*dat2;


	if(!(Xcg1=(long *) malloc(MAX_NUM_THREADS*CACHE_LINE_SIZE*sizeof(long)))) ERR_CRITICAL("Unable to allocate ranf storage\n");
	if(!(Xcg2=(long *) malloc(MAX_NUM_THREADS*CACHE_LINE_SIZE*sizeof(long)))) ERR_CRITICAL("Unable to allocate ranf storage\n");
	setall(P.seed1,P.seed2);

	P.DoBin=-1;
	if(P.DoHeteroDensity)
		{
		if(P.DoAdunitBoundaries)
			{
			fprintf(stderr,"Scanning population density file\n");
			if(!(dat=fopen(DensityFile,"rb"))) ERR_CRITICAL("Unable to open density file\n");
			fread_big(&(P.BinFileLen),sizeof(unsigned int),1,dat);
			if(P.BinFileLen==0xf0f0f0f0) //code for first 4 bytes of binary file ## NOTE - SHOULD BE LONG LONG TO COPE WITH BIGGER POPULATIONS
				{
				P.DoBin=1;
				fread_big(&(P.BinFileLen),sizeof(unsigned int),1,dat);
				if(!(BinFileBuf=(void *) malloc(P.BinFileLen*sizeof(bin_file)))) ERR_CRITICAL("Unable to allocate binary file buffer\n");
				fread_big(BinFileBuf,sizeof(bin_file),(size_t) P.BinFileLen,dat);
				BF=(bin_file *) BinFileBuf;
				fclose(dat);
				}
			else
				{
				P.DoBin=0;
				fclose(dat);
				if(!(dat=fopen(DensityFile,"r"))) ERR_CRITICAL("Unable to open density file\n");
				P.BinFileLen=UINT_MAX-1;
				}
			P.SpatialBoundingBox[0]=P.SpatialBoundingBox[1]=1e10;
			P.SpatialBoundingBox[2]=P.SpatialBoundingBox[3]=-1e10;
			s2=0;
			if(P.DoBin==0)
				{
				fgets(buf,2047,dat);
				if(feof(dat)) rn=P.BinFileLen;
				}
			for(rn=0;rn<P.BinFileLen;rn++)
				{
				if(P.DoBin==0)
					{
					sscanf(buf,"%lg\t%lg\t%lg\t%i\t%i",&x,&y,&t,&i2,&l);
					if(l/P.CountryDivisor!=i2) //temporarily changed this to 10000 from 100 to work with new admin codes - ggilani 30/05/2018, now changed to CountryDivisor - ggilani 13/05/2019
					{
						//fprintf(stderr,"# %lg %lg %lg %i %i\n",x,y,t,i2,l);
					}
					//if(i2==P.TargetCountry)
					//{
					//	fprintf(stderr,"# %lg %lg %lg %i %i\n",x,y,t,i2,l);
					//}
					fgets(buf,2047,dat);
					if(feof(dat)) rn=P.BinFileLen;
					}
				else
					{
					x=BF[rn].x;
					y=BF[rn].y;
					t=BF[rn].pop;
					i2=BF[rn].cnt;
					l=BF[rn].ad;
//					fprintf(stderr,"# %lg %lg %lg %i\t",x,y,t,l);

					}
				m=(l%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor;
				if(P.AdunitLevel1Lookup[m]>=0)
					if(AdUnits[P.AdunitLevel1Lookup[m]].id/P.AdunitLevel1Mask==l/P.AdunitLevel1Mask)
						{
						AdUnits[P.AdunitLevel1Lookup[m]].cnt_id=i2;
						s2+=t;
						if(x<P.SpatialBoundingBox[0]) P.SpatialBoundingBox[0]=x;
						if(x>=P.SpatialBoundingBox[2]) P.SpatialBoundingBox[2]=x+1e-6;
						if(y<P.SpatialBoundingBox[1]) P.SpatialBoundingBox[1]=y;
						if(y>=P.SpatialBoundingBox[3]) P.SpatialBoundingBox[3]=y+1e-6;
						}
				}
			if(!P.DoSpecifyPop) P.N=(int) s2;
			if(P.DoBin==0) fclose(dat);
			}

#ifdef COUNTRY_THAILAND
		P.width=P.SpatialBoundingBox[2]-P.SpatialBoundingBox[0];
		P.height=P.SpatialBoundingBox[3]-P.SpatialBoundingBox[1];
		P.ncw=(int) (P.width/P.cwidth);
		P.nch=(int) (P.height/P.cwidth);
		P.cwidth=P.width/((double) P.ncw);
		P.cheight=P.height/((double) P.nch);
#else
		P.cheight=P.cwidth;
		P.SpatialBoundingBox[0]=floor(P.SpatialBoundingBox[0]/P.cwidth)*P.cwidth;
		P.SpatialBoundingBox[1]=floor(P.SpatialBoundingBox[1]/P.cheight)*P.cheight;
		P.SpatialBoundingBox[2]=ceil(P.SpatialBoundingBox[2]/P.cwidth)*P.cwidth;
		P.SpatialBoundingBox[3]=ceil(P.SpatialBoundingBox[3]/P.cheight)*P.cheight;
		P.width=P.SpatialBoundingBox[2]-P.SpatialBoundingBox[0];
		P.height=P.SpatialBoundingBox[3]-P.SpatialBoundingBox[1];
		P.ncw=4*((int) ceil(P.width/P.cwidth/4));
		P.nch=4*((int) ceil(P.height/P.cheight/4));
		P.width=((double) P.ncw)*P.cwidth;
		P.height=((double) P.nch)*P.cheight;
		P.SpatialBoundingBox[2]=P.SpatialBoundingBox[0]+P.width;
		P.SpatialBoundingBox[3]=P.SpatialBoundingBox[1]+P.height;
#endif
		P.NC=P.ncw*P.nch;
		fprintf(stderr,"Adjusted bounding box = (%lg, %lg)- (%lg, %lg)\n",P.SpatialBoundingBox[0],P.SpatialBoundingBox[1],P.SpatialBoundingBox[2],P.SpatialBoundingBox[3]);
		fprintf(stderr,"Number of cells = %i (%i x %i)\n",P.NC,P.ncw,P.nch);
		fprintf(stderr,"Population size = %i \n",P.N);
		s=1;
		P.DoPeriodicBoundaries=0;
		}
	else
		{
		P.ncw=P.nch=(int) sqrt((double) P.NC);
		P.NC=P.ncw*P.nch;
		fprintf(stderr,"Number of cells adjusted to be %i (%i^2)\n",P.NC,P.ncw);
		s=floor(sqrt((double) P.N));
		P.SpatialBoundingBox[0]=P.SpatialBoundingBox[1]=0;
		P.SpatialBoundingBox[2]=P.SpatialBoundingBox[3]=s;
		P.N=(int) (s*s);
		fprintf(stderr,"Population size adjusted to be %i (%lg^2)\n",P.N,s);
		P.width=P.height=s;
		P.cwidth=P.width/((double) P.ncw);
		P.cheight=P.height/((double) P.nch);
		}
	P.NMC=P.NMCL*P.NMCL*P.NC;
	P.nmcw=P.ncw*P.NMCL;
	P.nmch=P.nch*P.NMCL;
	fprintf(stderr,"Number of micro-cells = %i\n",P.NMC);
	P.scalex=P.BitmapScale;
	P.scaley=P.BitmapAspectScale*P.BitmapScale;
	P.bwidth=(int) (P.width*(P.BoundingBox[2]-P.BoundingBox[0])*P.scalex);
	P.bwidth=(P.bwidth+3)/4;
	P.bwidth*=4;
	P.bheight=(int) (P.height*(P.BoundingBox[3]-P.BoundingBox[1])*P.scaley);
	P.bheight+=(4-P.bheight%4)%4;
	P.bheight2=P.bheight+20; // space for colour legend
	fprintf(stderr,"Bitmap width = %i\nBitmap height = %i\n",P.bwidth,P.bheight);
	P.bminx=(int) (P.width*P.BoundingBox[0]*P.scalex);
	P.bminy=(int) (P.height*P.BoundingBox[1]*P.scaley);
	P.mcwidth=P.cwidth/((double) P.NMCL);
	P.mcheight=P.cheight/((double) P.NMCL);
	for(i=0;i<P.NumSeedLocations;i++)
		{
		P.LocationInitialInfection[i][0]-=P.SpatialBoundingBox[0];
		P.LocationInitialInfection[i][1]-=P.SpatialBoundingBox[1];
		}
	//if we are importing cases to a specific location, adjust for spatial bounding box: ggilani 01/07/2015
	if(P.DoImportToSpecLocation)
	{
		//if we haven't actually specified it, coordinates will be set to -1001 and we can see them to the middle of the spatial bounding box
		if((P.ImportLocation[0]==-1000)&&(P.ImportLocation[1]==-1000))
		{
			P.ImportLocation[0]=0.5*(P.SpatialBoundingBox[0]+P.SpatialBoundingBox[2]); //set location for importation to the centre of the bounding box
			P.ImportLocation[1]=0.5*(P.SpatialBoundingBox[1]+P.SpatialBoundingBox[3]);
		}
		else
		{
			//adjust the given location to account for the bounding box
			P.ImportLocation[0]-=P.SpatialBoundingBox[0];
			P.ImportLocation[1]-=P.SpatialBoundingBox[1];
		}
	}
	t=dist2_raw(0,0,P.width,P.height);
	if(P.DoPeriodicBoundaries) t*=0.25;
	if(!(nKernel=(double *) calloc(NKR+1,sizeof(double)))) ERR_CRITICAL("Unable to allocate kernel storage\n");
	if(!(nKernelHR=(double *) calloc(NKR+1,sizeof(double)))) ERR_CRITICAL("Unable to allocate kernel storage\n");
	P.KernelDelta=t/NKR;
//	fprintf(stderr,"** %i %lg %lg %lg %lg | %lg %lg %lg %lg \n",P.DoUTM_coords,P.SpatialBoundingBox[0],P.SpatialBoundingBox[1],P.SpatialBoundingBox[2],P.SpatialBoundingBox[3],P.width,P.height,t,P.KernelDelta);
	fprintf(stderr,"Coords xmcell=%lg m   ymcell = %lg m\n",sqrt(dist2_raw(P.width/2,P.height/2,P.width/2+P.mcwidth,P.height/2)),sqrt(dist2_raw(P.width/2,P.height/2,P.width/2,P.height/2+P.mcheight)));
	P.KernelShape=P.MoveKernelShape;
	P.KernelScale=P.MoveKernelScale;
	P.KernelP3=P.MoveKernelP3;
	P.KernelP4=P.MoveKernelP4;
	P.KernelType=P.MoveKernelType;
	t2=0.0;
	SetupPopulation(DensityFile,SchoolFile,RegDemogFile);
	if(!(TimeSeries=(results *) calloc(P.NumSamples,sizeof(results)))) ERR_CRITICAL("Unable to allocate results storage\n");
	if(!(TSMeanE=(results *) calloc(P.NumSamples,sizeof(results)))) ERR_CRITICAL("Unable to allocate results storage\n");
	if(!(TSVarE=(results *) calloc(P.NumSamples,sizeof(results)))) ERR_CRITICAL("Unable to allocate results storage\n");
	if(!(TSMeanNE=(results *) calloc(P.NumSamples,sizeof(results)))) ERR_CRITICAL("Unable to allocate results storage\n");
	if(!(TSVarNE=(results *) calloc(P.NumSamples,sizeof(results)))) ERR_CRITICAL("Unable to allocate results storage\n");
	TSMean=TSMeanE;TSVar=TSVarE;
	for(l=0;l<2;l++)
		{
		for(i=0;i<P.NumSamples;i++)
			{
			TSMean[i].S=TSMean[i].I=TSMean[i].R=TSMean[i].D=TSMean[i].L=
				TSMean[i].incL=TSMean[i].incI=TSMean[i].incR=TSMean[i].incC=TSMean[i].incFI=TSMean[i].incDC=TSMean[i].cumDC=TSMean[i].cumDD=TSMean[i].cumSDB=
				TSMean[i].incTC=TSMean[i].cumT=TSMean[i].cumTP=TSMean[i].cumUT=TSMean[i].cumV=TSMean[i].ETU=TSMean[i].incETU=TSMean[i].incCT=TSMean[i].incCC= //added contact tracing, cases who are contacts
				TSMean[i].cumTmax=TSMean[i].cumVmax=TSMean[i].incD=TSMean[i].incHQ=TSMean[i].incAC=
				TSMean[i].incAH=TSMean[i].incAA=TSMean[i].incACS=TSMean[i].incAPC=
				TSMean[i].incAPA=TSMean[i].incAPCS=TSMean[i].Rdenom=0;
			TSVar[i].S=TSVar[i].I=TSVar[i].R=TSVar[i].D=TSVar[i].L=
				TSVar[i].incL=TSVar[i].incI=TSVar[i].incR=TSVar[i].incC=TSVar[i].incTC=TSVar[i].incD=0;
			for(j=0;j<NUM_PLACE_TYPES;j++) TSMean[i].PropPlacesClosed[j]=TSVar[i].PropPlacesClosed[j]=0;
			for(j=0;j<INFECT_TYPE_MASK;j++) TSMean[i].incItype[j]=TSMean[i].Rtype[j]=0;
			for(j=0;j<NUM_AGE_GROUPS;j++) TSMean[i].incCa[j]=TSMean[i].incIa[j]=TSMean[i].incDa[j]= TSMean[i].incDCa[j]= TSMean[i].incETUa[j]= TSMean[i].incVa[j]=TSMean[i].Rage[j]=0;
			for(j=0;j<P.EvolResistNumTypes;j++)
					TSMean[i].incI_resist[j]=TSVar[i].incI_resist[j]=
					TSMean[i].incC_resist[j]=TSVar[i].incC_resist[j]=
					TSMean[i].cumT_resist[j]=TSVar[i].cumT_resist[j]=0;
			for(j=0;j<2;j++)
					TSMean[i].incI_keyworker[j]=TSVar[i].incI_keyworker[j]=
					TSMean[i].incC_keyworker[j]=TSVar[i].incC_keyworker[j]=
					TSMean[i].cumT_keyworker[j]=TSVar[i].cumT_keyworker[j]=0;
			if(P.DoAdUnits)
				for(j=0;j<=P.NumAdunits;j++)
					TSMean[i].incI_adunit[j]=TSVar[i].incI_adunit[j]=
					TSMean[i].incC_adunit[j]=TSVar[i].incC_adunit[j]=
					TSMean[i].incDC_adunit[j]=TSVar[i].incDC_adunit[j]=//added detected cases here: ggilani 03/02/15
					TSMean[i].incETU_adunit[j]=TSVar[i].incETU_adunit[j]=
					TSMean[i].incCT_adunit[j]=TSVar[i].incCT_adunit[j]= //added contact tracing
					TSMean[i].incCC_adunit[j]=TSVar[i].incCC_adunit[j]= //added cases who are contacts: ggilani 28/05/2019
					TSMean[i].cumT_adunit[j]=TSVar[i].cumT_adunit[j]=0;
			}
		TSMean=TSMeanNE;TSVar=TSVarNE;
		}

	//added memory allocation and initialisation of infection event log, if DoRecordInfEvents is set to 1: ggilani - 10/10/2014
	if(P.DoRecordInfEvents)
	{
		if(!(InfEventLog=(events *) calloc(P.MaxInfEvents,sizeof(events)))) ERR_CRITICAL("Unable to allocate events storage\n");
		if(!(nEvents=(int *) calloc(1,sizeof(int)))) ERR_CRITICAL("Unable to allocate events storage\n");
	}

	P.CellPop2=((double) P.N)*((double) P.N)/(((double) P.NC)*((double) P.NC));
	//SaveAgeDistrib();

	if((P.DoCapitalCityEffect)&&(P.DoAdUnits)) //if including a capital city effect, determine whether a cell contains a capital city or not before initialising the kernel
	{
		DetermineCellsWithCapitalCities();
	}

	fprintf(stderr,"Initialising kernel...\n");
	InitKernel(0,1.0);	fprintf(stderr,"Initialising places...\n");
	if(P.DoPlaces)
		{
		if(P.LoadSaveNetwork==1)
			LoadPeopleToPlaces(NetworkFile);
		else
			AssignPeopleToPlaces();
		}
	if((P.DoPlaces)&&(P.LoadSaveNetwork==2))
		SavePeopleToPlaces(NetworkFile);
	//SaveDistribs();

	//setall(P.seed3,P.seed4);

	for (i = 0; i < P.N; i++)
	{
		Hosts[i].keyworker = 0;
		Hosts[i].hcw = 0;
		Hosts[i].flw = 0;
	}

	StratifyPlaces(); //set hcws within StratifyPlaces()
	for(i=0;i<P.NC;i++) 
		{
		Cells[i].S=Cells[i].n;
		Cells[i].L=Cells[i].I=Cells[i].R=0;
		//Cells[i].susceptible=Cells[i].members; //added this line
		}

	P.KeyWorkerNum=P.KeyWorkerIncHouseNum=m=l=0;
	if(P.DoPlaces)
		{
		while((m<P.KeyWorkerPopNum)&&(l<1000))
			{
			i=(int) (((double) P.N)*ranf_mt(0));
			if(Hosts[i].keyworker)
				l++;
			else
				{
				Hosts[i].keyworker=1;
				m++;
				P.KeyWorkerNum++;
				P.KeyWorkerIncHouseNum++;
				l=0;
				if(ranf_mt(0)<P.KeyWorkerHouseProp)
					{
					l2=Households[Hosts[i].hh].FirstPerson;
					m2=l2+Households[Hosts[i].hh].nh;
					for(j2=l2;j2<m2;j2++)
						if(!Hosts[j2].keyworker)
							{
							Hosts[j2].keyworker=1;
							P.KeyWorkerIncHouseNum++;
							}
					}
				}
			}
		for(j=0;j<NUM_PLACE_TYPES_NOAIR;j++)
			{
			m=l=0;
			while((m<P.KeyWorkerPlaceNum[j])&&(l<1000))
				{
				k=(int) (((double) P.Nplace[j])*ranf_mt(0));
				for(i2=0;(m<P.KeyWorkerPlaceNum[j])&&(i2<Places[j][k].n);i2++)
					{
					i=Places[j][k].members[i2];
					if((i<0)||(i>=P.N)) fprintf(stderr,"## %i # ",i);
					if((Hosts[i].keyworker)||(ranf_mt(0)>=P.KeyWorkerPropInKeyPlaces[j]))
						l++;
					else
						{
						Hosts[i].keyworker=1;
						m++;
						P.KeyWorkerNum++;
						P.KeyWorkerIncHouseNum++;
						l=0;
						l2=Households[Hosts[i].hh].FirstPerson;
						m2=l2+Households[Hosts[i].hh].nh;
						for(j2=l2;j2<m2;j2++)
							if((!Hosts[j2].keyworker)&&(ranf_mt(0)<P.KeyWorkerHouseProp))
								{
								Hosts[j2].keyworker=1;
								P.KeyWorkerIncHouseNum++;
								}
						}
					}
				}
			}
		if(P.KeyWorkerNum>0) fprintf(stderr,"%i key workers selected in total\n",P.KeyWorkerNum);
		if(P.DoAdUnits)
			{
			for(i=0;i<P.NumAdunits;i++) AdUnits[i].NP=0;
				for(j=0;j<P.PlaceTypeNum;j++)
					if(P.PlaceCloseAdunitPlaceTypes[j]>0)
						{
						for(k=0;k<P.Nplace[j];k++)
							AdUnits[Mcells[Places[j][k].mcell].adunit].NP++;
						}
			}
		}
	fprintf(stderr,"Places intialised.\n");
	UpdateProbs(0);
	if(P.DoAirports) SetupAirports();
	if(P.R0scale!=1.0)
		{
		P.HouseholdTrans*=P.R0scale;
		P.R0*=P.R0scale;
		for(j=0;j<P.PlaceTypeNum;j++)
			P.PlaceTypeTrans[j]*=P.R0scale;
		fprintf(stderr,"Rescaled transmission coefficients by factor of %lg\n",P.R0scale);
		}
	t=s=t2=0;
	for(i=0;i<MAX_HOUSEHOLD_SIZE;i++)
		{
		t+=((double) (i+1))*(P.HouseholdSizeDistrib[0][i]-t2);
		t2=P.HouseholdSizeDistrib[0][i];
		}
	t2=s=0;
	s3=1.0;
#pragma omp parallel for private(i,s2,j,k,q,l,d,y,m,tn) schedule(static,1) reduction(+:s,t2) //schedule(static,1000)
	for(tn=0;tn<P.NumThreads;tn++) //changed this looping to allow for multi-threaded random numbers
	{
		for(i=tn;i<P.N;i+=P.NumThreads)
		{
		if(P.InfectiousnessSD==0)
			Hosts[i].infectiousness=P.AgeInfectiousness[HOST_AGE_GROUP(i)];
		else
			Hosts[i].infectiousness=P.AgeInfectiousness[HOST_AGE_GROUP(i)]*gen_gamma_mt(P.InfectiousnessGamA,P.InfectiousnessGamR,tn); //made this multi-threaded: 28/11/14
			//Hosts[i].infectiousness=P.AgeInfectiousness[HOST_AGE_GROUP(i)]*gen_beta_mt(P.InfectiousnessBetaA,P.InfectiousnessBetaB,tn);
		q=P.ProportionSymptomatic[HOST_AGE_GROUP(i)];
		if(ranf_mt(tn)<q) //made this multi-threaded: 28/11/14
			Hosts[i].infectiousness=-P.SymptInfectiousness*Hosts[i].infectiousness;
		j=(int) floor((q=ranf_mt(tn)*CDF_RES)); //made this multi-threaded: 28/11/14
		q-=((double) j);
		Hosts[i].recovery_time=(unsigned short int) floor(0.5-(P.InfectiousPeriod*log(q*P.infectious_icdf[j+1]+(1.0-q)*P.infectious_icdf[j])/P.TimeStep));
		//adding a step here to see, if we are doing funeral transmission, the host dies: ggilani 14/11/14
		if(P.DoFuneralTransmission)
		{
			if(P.DoMortality)
			{
				if (P.DoEventMortality)
				{
					if (ranf_mt(tn) <= P.RecoveryProb[(int)ceil(Hosts[i].recovery_time * P.TimeStep)]) Hosts[i].to_die = 0; //made this multi-threaded: 28/11/14
					else Hosts[i].to_die = 1;
				}
				else
				{
					if (P.DoAgeMortality)
					{
						probMort = P.AgeMortality[HOST_AGE_GROUP(i)];
					}
					else
					{
						probMort = P.DiseaseMortality;
					}
					if (ranf_mt(tn) < probMort)
					{
						Hosts[i].to_die = 1;
						Hosts[i].recovery_time = (unsigned short int)(P.LethalInfectiousPeriod * (double)Hosts[i].recovery_time);//lethal infectious period
					}
					else Hosts[i].to_die = 0;
				}
			}
		}
		//adding this to allow for hospital transmission
		if (P.IncludeHospitalPlaceType)
		{
			//calculate number of hospital beds to be utilised over the simulation, which will allow us to estimate the number of people who actually seek healthcare are able to access hospital care
			//due to capacity restraints
			if (P.DoEventMortality)
			{
				inf_period = P.InfectiousPeriod;
			}
			else
			{
				inf_period = P.DiseaseMortality * P.LethalInfectiousPeriod * P.InfectiousPeriod + (1.0 - P.DiseaseMortality) * P.InfectiousPeriod;
			}
			BedCapacity = (int)((double)P.HospCaseCapacity * double(P.Nplace[P.HospPlaceTypeNum] * ((double)P.SampleTime / inf_period)));
			ProbHosp = ((double)BedCapacity) / (P.PropHospSeekPreOutbreak * (double)P.N);
			Hosts[i].hospitalised = (ranf_mt(tn) < P.PropHospSeekPreOutbreak);
			//Hosts[i].hospitalised = (ranf_mt(tn) < ProbHosp);
		}
		if(P.DoHouseholds)
			{
			s2=P.TimeStep*P.HouseholdTrans*fabs(Hosts[i].infectiousness)*P.HouseholdDenomLookup[Households[Hosts[i].hh].nhr-1];
			d=1.0;l=(int) Hosts[i].recovery_time;
			for(k=0;k<l;k++) {y=1.0-s2*P.infectiousness[k];d*=((y<0)?0:y);}
			//adding effect of funeral transmission to overall household transmission: ggilani - 14/11/14
			if((P.DoFuneralTransmission)&&(Hosts[i].to_die))
			{
				//add funeral transmission time to host recovery time parameter l
				l+=(int)(P.FuneralTransmissionDuration/P.TimeStep);
				//include infectiousness of funeral transmission period, allowing for increased infectiousness of funerals
				for(;k<l;k++) {y=1.0-s2*P.RelativeInfectiousnessFuneral*P.infectiousness[k];d*=((y<0)?0:y);}
			}
			if(!((Hosts[i].nc_plus_hh_disabled & HH_DISABLED)))
				{
				l=Households[Hosts[i].hh].FirstPerson;
				m=l+Households[Hosts[i].hh].nh;
				for(k=l;k<m;k++) if((Hosts[k].inf==0)&&(k!=i)) s+=(1-d)*P.AgeSusceptibility[HOST_AGE_GROUP(i)];
				}
			}
		q=(P.LatentToSymptDelay>Hosts[i].recovery_time*P.TimeStep)?Hosts[i].recovery_time*P.TimeStep:P.LatentToSymptDelay;
		s2=fabs(Hosts[i].infectiousness)*P.RelativeSpatialContact[HOST_AGE_GROUP(i)]*P.TimeStep;
		l=(int) (q/P.TimeStep);
		for(k=0;k<l;k++) t2+=s2*P.infectiousness[k];
		s2*=((Hosts[i].infectiousness<0)?P.SymptSpatialContactRate:1);
		l=(int) Hosts[i].recovery_time;
		for(;k<l;k++) t2+=s2*P.infectiousness[k];
		//adding effect of funeral transmission spatially here: ggilani - 14/11/14
		if((P.DoFuneralTransmission)&&(Hosts[i].to_die))
		{
			//add funeral transmission time to host recovery time parameter l
			l+=(int)(P.FuneralTransmissionDuration/P.TimeStep);
			//include infectiousness of funeral transmission period, allowing for increased infectiousness of funerals
			for(;k<l;k++) t2+=s2*P.RelativeInfectiousnessFuneral*P.infectiousness[k];
		}
		}
	}
	t2*=(s3/((double) P.N));
	s/=((double) P.N);
	fprintf(stderr,"Household mean size=%lg\nHousehold R0=%lg\n",t,P.R0household=s);
	t=x=y=0;
	if(P.DoPlaces)
	  for(j=0;j<P.PlaceTypeNum;j++)
		if((j!=HOTEL_PLACE_TYPE))//&&(j!=P.HospPlaceTypeNum)) //added this so that don't double count places. Or we can do this another way and figure out who is going to hospital
			{
			if (j == 1) //check this
			{
				k = 0;
			}
//#pragma omp parallel for private(i,k,d,q,s2,s3,t3,l,m,x,y) schedule(static,1000) reduction(+:t)
			for(i=0;i<P.N;i++)
				{
				k=Hosts[i].PlaceLinks[j];
				if ((P.IncludeHospitalPlaceType)&&((j==P.HospPlaceTypeNum)&&(Hosts[i].hospitalised==0)))
				{
					k = PERSON_NOT_IN_PLACE_AT_THAT_TIME;
				}
				if (Hosts[i].hospitalised)
				{
					contact_scale = 0.5;
				}
				else
				{
					contact_scale = 1.0;
				}
				if(k>=0)
					{
					q=(P.LatentToSymptDelay>Hosts[i].recovery_time*P.TimeStep)?Hosts[i].recovery_time*P.TimeStep:P.LatentToSymptDelay;
					s2=fabs(Hosts[i].infectiousness)*P.TimeStep*P.PlaceTypeTrans[j]*contact_scale;
					if ((P.IncludeHospitalPlaceType) && (j==P.HospPlaceTypeNum))
					{
						x= s2/Places[j][k].nhcws;
					}
					else
					{
						x = s2 / P.PlaceTypeGroupSizeParam1[j];
					}

					d=1.0;l=(int) (q/P.TimeStep);
					for(m=0;m<l;m++) {y=1.0-x*P.infectiousness[m];d*=((y<0)?0:y);}

					if ((P.IncludeHospitalPlaceType) && (j == P.HospPlaceTypeNum))
					{
						s3 = ((double)(Places[j][k].nhcws - 1));
					}
					else
					{
						s3 = ((double)(Places[j][k].group_size[Hosts[i].PlaceGroupLinks[j]] - 1));
					}

					x*=((Hosts[i].infectiousness<0)?(P.SymptPlaceTypeContactRate[j]*(1-P.SymptPlaceTypeWithdrawalProp[j])):1);
					l=(int) Hosts[i].recovery_time;
					for(;m<l;m++) {y=1.0-x*P.infectiousness[m];d*=((y<0)?0:y);}
					//adding effect of funeral transmission here for within group places: ggilani 14/11/14
					if((P.DoFuneralTransmission)&&(Hosts[i].to_die))
					{
						//add funeral transmission time to host recovery time parameter l
						l+=(int)(P.FuneralTransmissionDuration/P.TimeStep);
						//include infectiousness of funeral transmission period, allowing for increased infectiousness of funerals
						for(;m<l;m++) {y=1.0-x*P.RelativeInfectiousnessFuneral*P.infectiousness[m];d*=((y<0)?0:y);}
					}
					t3=d;
					x=P.PlaceTypePropBetweenGroupLinks[j]*s2/((double) Places[j][k].n);
					d=1.0;l=(int) (q/P.TimeStep);
					for(m=0;m<l;m++) {y=1.0-x*P.infectiousness[m];d*=((y<0)?0:y);}
					x*=((Hosts[i].infectiousness<0)?(P.SymptPlaceTypeContactRate[j]*(1-P.SymptPlaceTypeWithdrawalProp[j])):1);
					l=(int) Hosts[i].recovery_time;
					for(;m<l;m++) {y=1.0-x*P.infectiousness[m];d*=((y<0)?0:y);} //check this
					//adding effect of funeral transmission here for between groups within places: ggilani 14/11/14
					if((P.DoFuneralTransmission)&&(Hosts[i].to_die))
					{
						//add funeral transmission time to host recovery time parameter l
						l+=(int)(P.FuneralTransmissionDuration/P.TimeStep);
						//include infectiousness of funeral transmission period, allowing for increased infectiousness of funerals
						for(;m<l;m++) {y=1.0-x*P.RelativeInfectiousnessFuneral*P.infectiousness[m];d*=((y<0)?0:y);}
					}
					if ((P.IncludeHospitalPlaceType) && (j == P.HospPlaceTypeNum))
					{
						t += (1 - t3 * d) * s3;// +(1 - d) * (((double)(Places[j][k].nhcws - 1)) - s3); //could get rid of this
					}
					else
					{
						t += (1 - t3 * d) * s3 + (1 - d) * (((double)(Places[j][k].n - 1)) - s3);
					}
					}
				}
			fprintf(stderr,"%lg  ",t/((double) P.N));
			}
#pragma omp parallel for private(i) schedule(static,500) reduction(+:x,y)
	for(i=0;i<P.N;i++)
		{
		x+=Hosts[i].recovery_time*P.TimeStep;
		y+=Hosts[i].recovery_time;
		Hosts[i].recovery_time=0;
		Hosts[i].to_die=0;
		Hosts[i].hospitalised = 0;
		}
	t/=((double) P.N);
	x/=((double) P.N);
	y/=((double) P.N);
	fprintf(stderr,"R0 for places = %lg\nR0 for random spatial = %lg\nOverall R0=%lg\n",P.R0places=t,P.R0spatial=P.R0-s-t,P.R0);
	fprintf(stderr,"Mean infectious period (sampled) = %lg (%lg)\n",x,y);
	if(P.DoSI)
		P.LocalBeta=(P.R0/t2-s-t);
	else
		P.LocalBeta=(P.R0-s-t)/t2;
	if((P.LocalBeta<0)||(!P.DoSpatial))
		{
		P.LocalBeta=P.R0spatial=0;
		fprintf(stderr,"Reset spatial R0 to 0\n");
		}
	fprintf(stderr,"LocalBeta = %lg\n",P.LocalBeta);
	TSMean=TSMeanNE;TSVar=TSVarNE;
	fprintf(stderr,"Calculated approx cell probabilities\n");
	for(i=0;i<INFECT_TYPE_MASK;i++) inftype_av[i]=0;
	for(i=0;i<MAX_COUNTRIES;i++) infcountry_av[i]=infcountry_num[i]=0;
	for(i=0;i<MAX_SEC_REC;i++) 
		for(j=0;j<MAX_GEN_REC;j++)
			indivR0_av[i][j]=0;
	for(i=0;i<=MAX_HOUSEHOLD_SIZE;i++)
		for(j=0;j<=MAX_HOUSEHOLD_SIZE;j++)
			inf_household_av[i][j]=case_household_av[i][j]=0;
	fprintf(stderr,"Model configuration complete.\n");
	DoInitUpdateProbs=1;
	for(i=0;i<P.NC;i++)	Cells[i].tot_treat=1;  //This makes sure InitModel intialises the cells.
	P.NRactE=P.NRactNE=0;
	if(P.ESocProportionCompliant>0)
		{
		for(i=0;i<P.N;i++)
			Hosts[i].esocdist_comply=(ranf()<P.ESocProportionCompliant)?1:0;
		}
//	else
//		{
//		for(i=0;i<P.N;i++) Hosts[i].esocdist_comply=0;
//		}
	if(P.OutputBitmap) InitBMHead();
	if(P.DoMassVacc)
		{
		if(!(State.mvacc_queue=(int *) calloc(P.N,sizeof(int)))) ERR_CRITICAL("Unable to allocate host storage\n");
		for(i=j=0;i<P.N;i++)
			{
			if((HOST_AGE_YEAR(i)>=P.VaccPriorityGroupAge[0])&&(HOST_AGE_YEAR(i)<=P.VaccPriorityGroupAge[1]))
				{
				if(Hosts[i].vacc_accept<P.VaccProp)
					State.mvacc_queue[j++]=i;
				}
			}
		k=j;
		for(i=0;i<P.N;i++)
			{
			if((HOST_AGE_YEAR(i)<P.VaccPriorityGroupAge[0])||(HOST_AGE_YEAR(i)>P.VaccPriorityGroupAge[1]))
				{
				if(Hosts[i].vacc_accept<P.VaccProp)
					State.mvacc_queue[j++]=i;
				}
			}
		State.n_mvacc=j;
		fprintf(stderr,"Number to be vaccinated=%i\n",State.n_mvacc);
		for(i=0;i<2;i++)
			{
			for(j=0;j<k;j++)
				{
				l=(int) (ranf()*((double) k));
				m=State.mvacc_queue[j];
				State.mvacc_queue[j]=State.mvacc_queue[l];
				State.mvacc_queue[l]=m;
				}
			for(j=k;j<State.n_mvacc;j++)
				{
				l=k+((int) (ranf()*((double) (State.n_mvacc-k))));
				m=State.mvacc_queue[j];
				State.mvacc_queue[j]=State.mvacc_queue[l];
				State.mvacc_queue[l]=m;
				}
			}
		fprintf(stderr,"Configured mass vaccination queue.\n");
		}

	PeakHeightSum=PeakHeightSS=PeakTimeSum=PeakTimeSS=0;
	i=(P.ncw/2)*P.nch+P.nch/2;
	j=(P.ncw/2+2)*P.nch+P.nch/2;
	fprintf(stderr,"UTM dist horiz=%lg %lg\n",sqrt(dist2_cc(Cells+i,Cells+j)),sqrt(dist2_cc(Cells+j,Cells+i)));
	j=(P.ncw/2)*P.nch+P.nch/2+2;
	fprintf(stderr,"UTM dist vert=%lg %lg\n",sqrt(dist2_cc(Cells+i,Cells+j)),sqrt(dist2_cc(Cells+j,Cells+i)));
	j=(P.ncw/2+2)*P.nch+P.nch/2+2;
	fprintf(stderr,"UTM dist diag=%lg %lg\n",sqrt(dist2_cc(Cells+i,Cells+j)),sqrt(dist2_cc(Cells+j,Cells+i)));

	//if(P.OutputBitmap)
	//{
	//	CaptureBitmap(0,0);
	//	OutputBitmap(0.0,0);
	//}

}

void SetupPopulation(char *DensityFile,char *SchoolFile, char *RegDemogFile)
	{
	int i,j,k,l,m,i2,j2,last_i,mr,ad,tn,*mcl;
	unsigned int rn,rn2;
	double t,s,x,y,xh,yh,maxd,s2,CumAgeDist[NUM_AGE_GROUPS+1],*income_distrib;
	char buf[4096],*col;
	const char delimiters[] = " \t,";
	FILE *dat,*dat2;
	bin_file rec,*BinFileOutBuf,*BFO;
	double temp_rep_rate,temp_vacc_accept; //added this to add case detection rate per household

	if(!(Cells=(cell *) calloc(P.NC,sizeof(cell)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if(!(Mcells=(microcell *) calloc(P.NMC,sizeof(microcell)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if(!(mcell_num=(int *) malloc(P.NMC*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if(!(mcell_dens=(double *) malloc(P.NMC*sizeof(double)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if(!(mcell_country=(int *) malloc(P.NMC*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if(!(mcell_adunits=(int *) malloc(P.NMC*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if(P.DoPlaces)
		if(!(Places=(place **) malloc(P.PlaceTypeNum*sizeof(place *)))) ERR_CRITICAL("Unable to allocate place storage\n");
	for(j=0;j<P.NMC;j++)
		{
		Mcells[j].n=0;
		mcell_adunits[j]=-1;
		mcell_dens[j]=0;
		mcell_num[j]=mcell_country[j]=0;
		//also set road connection to zero: ggilani 12/02/15
		Mcells[j].road_connection=0;
		}
	if(P.DoAdUnits)
		for(i=0;i<MAX_ADUNITS;i++)
			P.PopByAdunit[i][0]=P.PopByAdunit[i][1]=0;
	if(P.DoHeteroDensity)
		{
		if(!P.DoAdunitBoundaries) P.NumAdunits=0;
		if(P.DoBin==-1) 
			{
			if(!(dat=fopen(DensityFile,"rb"))) ERR_CRITICAL("Unable to open density file\n");
			fread_big(&(P.BinFileLen),sizeof(unsigned int),1,dat);
			if(P.BinFileLen==0xf0f0f0f0) //code for first 4 bytes of binary file
				{
				fprintf(stderr,"Reading binary population density file...\n");
				P.DoBin=1;
				fread_big(&(P.BinFileLen),sizeof(unsigned int),1,dat);
				if(!(BinFileBuf=(void *) malloc(((size_t) P.BinFileLen)*sizeof(bin_file)))) ERR_CRITICAL("Unable to allocate binary file buffer\n");
				fread_big(BinFileBuf,sizeof(bin_file),(size_t) P.BinFileLen,dat);
				BF=(bin_file *) BinFileBuf;
				fclose(dat);
				}
			else
				{
				fclose(dat);
				P.DoBin=0;
				}
			}
		if(P.DoBin==0)
			{
			P.BinFileLen=UINT_MAX-1;
			fprintf(stderr,"Reading ASCII population density file...\n");
			if(!(dat=fopen(DensityFile,"r"))) ERR_CRITICAL("Unable to open density file\n");
			}
//		if(!(dat2=fopen("EnvTest.txt","w"))) ERR_CRITICAL("Unable to open test file\n");
		if(P.DoBin==1)
			fprintf(stderr,"Binary density file contains %i cells.\n",(int) P.BinFileLen);
		else
			{
			fgets(buf,2047,dat);
			if(feof(dat)) rn=P.BinFileLen;
			}
		for(rn=rn2=mr=0;rn<P.BinFileLen;rn++)
			{
			if(P.DoAdUnits)
				{
				if(P.DoBin==1)
					{
					x=BF[rn].x;y=BF[rn].y;t=BF[rn].pop;i2=BF[rn].cnt;j2=BF[rn].ad; //changed from i to rn to loop over indices properly
					rec=BF[rn];
					}
				else
					sscanf(buf,"%lg\t%lg\t%lg\t%i\t%i",&x,&y,&t,&i2,&j2);
				m=(j2%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor;
				if(P.DoAdunitBoundaries)
					{
					if(P.AdunitLevel1Lookup[m]>=0)
						{
						if(j2/P.AdunitLevel1Mask==AdUnits[P.AdunitLevel1Lookup[m]].id/P.AdunitLevel1Mask) 
							{
							k=1;
							AdUnits[P.AdunitLevel1Lookup[m]].cnt_id=i2;
							}
						}
					else
						k=0;
					}
				else
					{
					k=1;
					if(P.AdunitLevel1Lookup[m]<0) 
						{
						P.AdunitLevel1Lookup[m]=P.NumAdunits;
						AdUnits[P.NumAdunits].id=j2;
						AdUnits[P.NumAdunits].cnt_id=i2;
						P.NumAdunits++;
						if(P.NumAdunits>=MAX_ADUNITS) ERR_CRITICAL("Total number of administrative units exceeds MAX_ADUNITS\n");
						}
					else
						{
						AdUnits[P.AdunitLevel1Lookup[m]].cnt_id=i2;
						}
					}
				}
			else
				{
				k=1;
				if(P.DoBin==1)
					{
					x=BF[i].x;y=BF[i].y;t=BF[i].pop;i2=BF[i].cnt;j2=BF[i].ad;
					rec=BF[rn];
					}
				else
					{
					sscanf(buf,"%lg\t%lg\t%lg\t%i",&x,&y,&t,&i2);
					j2=0;
					rec.x=x;rec.y=y;rec.pop=t;rec.cnt=i2;rec.ad=j2;
					}
				}
			if(P.DoBin==0)
				{
				fgets(buf,2047,dat);
				if(feof(dat)) rn=P.BinFileLen;
				}
			if((k)&&(x>=P.SpatialBoundingBox[0])&&(y>=P.SpatialBoundingBox[1])&&(x<P.SpatialBoundingBox[2])&&(y<P.SpatialBoundingBox[3]))
				{
				j=(int) floor((x-P.SpatialBoundingBox[0])/P.mcwidth+0.1);
				k=(int) floor((y-P.SpatialBoundingBox[1])/P.mcheight+0.1);
				l=j*P.nmch+k;
				if(l<P.NMC)
					{
					mr++;
					mcell_dens[l]+=t;
					mcell_country[l]=i2;
					//fprintf(stderr,"mcell %i, country %i, pop %lg\n",l,i2,t);
					mcell_num[l]++;
					if(P.DoAdUnits)
						{
						mcell_adunits[l]=P.AdunitLevel1Lookup[m];
						if(mcell_adunits[l]<0) fprintf(stderr,"Cell %i has adunits<0\n",l);
						P.PopByAdunit[P.AdunitLevel1Lookup[m]][0]+=t;
						}
					else
						mcell_adunits[l]=0;
					if((P.OutputDensFile)&&(P.DoBin)&&(mcell_adunits[l]>=0))
						{
						if(rn2<rn) BF[rn2]=rec;
						rn2++;
						}
					}
				}
			}
//		fclose(dat2);
		fprintf(stderr,"%i valid cells read from density file.\n",mr);
		if((P.OutputDensFile)&&(P.DoBin)) P.BinFileLen=rn2;
		if(P.DoBin==0)
			{
			fclose(dat);
			if(P.OutputDensFile)
				{
				P.DoBin=1;
				P.BinFileLen=0;
				for(l=0;l<P.NMC;l++)
					if(mcell_adunits[l]>=0) P.BinFileLen++;
				if(!(BinFileBuf=(void *) malloc(P.BinFileLen*sizeof(bin_file)))) ERR_CRITICAL("Unable to allocate binary file buffer\n");
				BF=(bin_file *) BinFileBuf;
				fprintf(stderr,"Binary density file should contain %i cells.\n",(int) P.BinFileLen);
				rn=0;
				for(l=0;l<P.NMC;l++)
					if(mcell_adunits[l]>=0)
						{
						BF[rn].x=(double) (P.mcwidth*(((double) (l/P.nmch))+0.5))+P.SpatialBoundingBox[0]; //x
						BF[rn].y=(double) (P.mcheight*(((double) (l%P.nmch))+0.5))+P.SpatialBoundingBox[1]; //y
						BF[rn].ad=(P.DoAdUnits)?(AdUnits[mcell_adunits[l]].id):0;
						BF[rn].pop=mcell_dens[l];
						BF[rn].cnt=mcell_country[l];
						rn++;
						}
				}
			}

		if(P.OutputDensFile)
			{
			if(!(dat2=fopen(OutDensFile,"wb"))) ERR_CRITICAL("Unable to open output density file\n");
			rn=0xf0f0f0f0;
			fwrite_big((void *) &rn,sizeof(unsigned int),1,dat2);
			fprintf(stderr,"Saving population density file with NC=%i...\n",(int) P.BinFileLen);
			fwrite_big((void *) &(P.BinFileLen),sizeof(unsigned int),1,dat2);
			fwrite_big(BinFileBuf,sizeof(bin_file),(size_t) P.BinFileLen,dat2);
			fclose(dat2);
			}
		if(P.DoBin==1) free(BinFileBuf);
		fprintf(stderr,"Population files read.\n");
		maxd=0;
		for(i=0;i<P.NMC;i++)
			{
			if(mcell_num[i]>0)
				{
				mcell_dens[i]/=((double) mcell_num[i]);
				Mcells[i].country=(unsigned short) mcell_country[i];
				if(P.DoAdUnits)
					Mcells[i].adunit=mcell_adunits[i];
				else
					Mcells[i].adunit=0;
				}
			else
				Mcells[i].adunit=-1;
			maxd+=mcell_dens[i];
			}
		}
	else
		{
		for(i=0;i<P.NMC;i++)
			{
			mcell_dens[i]=1.0;
			Mcells[i].country=1;
			}
		maxd=((double) P.NMC);
		}
	if(!P.DoAdUnits) P.NumAdunits=1;
	if((P.DoAdUnits)&&(P.DoAdunitDemog))
		{
		if(!(State.InvAgeDist=(int **) malloc(P.NumAdunits*sizeof(int *)))) ERR_CRITICAL("Unable to allocate InvAgeDist storage\n");
		for(i=0;i<P.NumAdunits;i++) 
			if(!(State.InvAgeDist[i]=(int *) malloc(1000*sizeof(int)))) ERR_CRITICAL("Unable to allocate InvAgeDist storage\n");
		if(!(dat=fopen(RegDemogFile,"r"))) ERR_CRITICAL("Unable to open regional demography file\n");
		for(k=0;k<P.NumAdunits;k++)
			{
			for(i=0;i<NUM_AGE_GROUPS;i++)
				P.PropAgeGroup[k][i]=0;
			for(i=0;i<MAX_HOUSEHOLD_SIZE;i++)
				P.HouseholdSizeDistrib[k][i]=0;
			P.PopByAdunit[k][1]=0;
			}
		while(!feof(dat))
			{
			fgets(buf,2047,dat);
			col=strtok(buf,delimiters);
			sscanf(col,"%i",&l);
			m=(l%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor;
			k=P.AdunitLevel1Lookup[m];
			if(k>=0)
				if(l/P.AdunitLevel1Mask==AdUnits[k].id/P.AdunitLevel1Mask)
					{
					col=strtok(NULL,delimiters);
					sscanf(col,"%lg",&x);
					P.PopByAdunit[k][1]+=x;
					t=0;
					for(i=0;i<NUM_AGE_GROUPS;i++)
						{
						col=strtok(NULL,delimiters);
						sscanf(col,"%lg",&s);
						P.PropAgeGroup[k][i]+=s;
						}
					col=strtok(NULL,delimiters);
					if(P.DoHouseholds)
						{
						sscanf(col,"%lg",&y);
						for(i=0;i<MAX_HOUSEHOLD_SIZE;i++)
							{
							col=strtok(NULL,delimiters);
							sscanf(col,"%lg",&s);
							P.HouseholdSizeDistrib[k][i]+=y*s;
							}
						}
					}
			}
		fclose(dat);
		for(k=0;k<P.NumAdunits;k++)
			{
			t=0;
			for(i=0;i<NUM_AGE_GROUPS;i++)
				t+=P.PropAgeGroup[k][i];
			CumAgeDist[0]=0;
			for(i=1;i<=NUM_AGE_GROUPS;i++)
				{
				P.PropAgeGroup[k][i-1]/=t;
				CumAgeDist[i]=CumAgeDist[i-1]+P.PropAgeGroup[k][i-1];
				}
			for(i=j=0;i<1000;i++)
				{
				t=((double) i)/1000;
				while(t>=CumAgeDist[j+1]) j++; 
				t=AGE_GROUP_WIDTH*(((double) j)+(t-CumAgeDist[j])/(CumAgeDist[j+1]-CumAgeDist[j]));
				State.InvAgeDist[k][i]=(int) t;
				}
			State.InvAgeDist[k][1000-1]=NUM_AGE_GROUPS*AGE_GROUP_WIDTH-1;
			if(P.DoHouseholds)
				{
				t=0;
				for(i=0;i<MAX_HOUSEHOLD_SIZE;i++)
					t+=P.HouseholdSizeDistrib[k][i];
				P.HouseholdSizeDistrib[k][0]/=t;
				for(i=1;i<MAX_HOUSEHOLD_SIZE-1;i++)
					P.HouseholdSizeDistrib[k][i]=P.HouseholdSizeDistrib[k][i]/t+P.HouseholdSizeDistrib[k][i-1];
				P.HouseholdSizeDistrib[k][MAX_HOUSEHOLD_SIZE-1]=1.0;
				}
			else
				{
				for(i=0;i<MAX_HOUSEHOLD_SIZE-1;i++)
					P.HouseholdSizeDistrib[k][i]=1.0;
				}
			}
		}
	else
		{
		if(!(State.InvAgeDist=(int **) malloc(sizeof(int *)))) ERR_CRITICAL("Unable to allocate InvAgeDist storage\n");
		if(!(State.InvAgeDist[0]=(int *) malloc(1000*sizeof(int)))) ERR_CRITICAL("Unable to allocate InvAgeDist storage\n");
		CumAgeDist[0]=0;
		for(i=1;i<=NUM_AGE_GROUPS;i++)
			CumAgeDist[i]=CumAgeDist[i-1]+P.PropAgeGroup[0][i-1];
		for(i=j=0;i<1000;i++)
			{
			t=((double) i)/1000;
			if(t>=CumAgeDist[j+1]) j++; 
			t=AGE_GROUP_WIDTH*(((double) j)+(t-CumAgeDist[j])/(CumAgeDist[j+1]-CumAgeDist[j]));
			State.InvAgeDist[0][i]=(int) t;
			}
		State.InvAgeDist[0][1000-1]=NUM_AGE_GROUPS*AGE_GROUP_WIDTH-1;
		}
	if((P.DoAdUnits)&&(P.DoAdunitDemog)&&(P.DoCorrectAdunitPop))
		{
		for(i=0;i<P.NumAdunits;i++)
			fprintf(stderr,"%i\t%i\t%lg\t%lg\n",i,(AdUnits[i].id%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor,P.PropAgeGroup[i][0],P.HouseholdSizeDistrib[i][0]);
		maxd=0;
		for(i=0;i<P.NMC;i++)
			{
			if(mcell_num[i]>0)
				mcell_dens[i]*=P.PopByAdunit[mcell_adunits[i]][1]/(1e-10+P.PopByAdunit[mcell_adunits[i]][0]);
			maxd+=mcell_dens[i];
			}
		t=0;
		for(i=0;i<P.NumAdunits;i++)
			t+=P.PopByAdunit[i][1];
		i=P.N;
		P.N=(int) t;
		fprintf(stderr,"Population size reset from %i to %i\n",i,P.N);
		}
	t=1.0;m=0;
	P.NMCP=0;
	for(i=0;i<(P.NMC-1);i++)
		{
		s=mcell_dens[i]/maxd/t;
		if(s>1.0) s=1.0;
		m+=(Mcells[i].n=(int) ignbin_mt((long) (P.N-m),s,0));
		t-=mcell_dens[i]/maxd;
		if(Mcells[i].n>0) P.NMCP++;
		}
	Mcells[P.NMC-1].n=P.N-m;
	if(Mcells[P.NMC-1].n>0) P.NMCP++;
	free(mcell_num);
	free(mcell_country);
	free(mcell_adunits);




	if(!(McellLookup=(microcell **) malloc(P.NMCP*sizeof(microcell *)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if(!(mcl=(int *) malloc(P.N*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	State.CellMemberArray=mcl;
	P.NCP=0;
	for(i=i2=j2=0;i<P.NC;i++)
		{
		Cells[i].n=0;
		k=(i/P.nch)*P.NMCL*P.nmch+(i%P.nch)*P.NMCL;
		Cells[i].members=mcl+j2;
		for(l=0;l<P.NMCL;l++)
			for(m=0;m<P.NMCL;m++)
				{
				j=k+m+l*P.nmch;
				if(Mcells[j].n>0)
					{
					Mcells[j].members=mcl+j2;
					//if(!(Mcells[j].members=(int *) calloc(Mcells[j].n,sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n"); //replaced line above with this to ensure members don't get mixed across microcells
					McellLookup[i2++]=Mcells+j;
					Cells[i].n+=Mcells[j].n;
					j2+=Mcells[j].n;
					}
				}
		if(Cells[i].n>0) P.NCP++;
		}
	fprintf(stderr,"Number of hosts assigned = %i\n",j2);
	if(!P.DoAdUnits) P.AdunitLevel1Lookup[0]=0;
	fprintf(stderr,"Number of cells with non-zero population = %i\n",P.NCP);
	fprintf(stderr,"Number of microcells with non-zero population = %i\n",P.NMCP);

	if(!(CellLookup=(cell **) malloc(P.NCP*sizeof(cell *)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if(!(mcl=(int *) malloc(P.N*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	State.CellSuscMemberArray=mcl;
	i2=k=0;
	for(j=0;j<P.NC;j++)
		if(Cells[j].n>0)
			{
			CellLookup[i2++]=Cells+j;
			Cells[j].susceptible=mcl+k;
			k+=Cells[j].n;
			}
	if(i2>P.NCP) fprintf(stderr,"######## Over-run on CellLookup array NCP=%i i2=%i ###########\n",P.NCP,i2);
	if(!(RevCellLookup=(int *) malloc(P.NC*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	i2=0;
	for(j=0;j<P.NC;j++)
		if(Cells[j].n>0)
			RevCellLookup[j]=i2++;
		else
			RevCellLookup[j]=-1;

	if(!(Hosts=(person *) malloc(P.N*sizeof(person)))) ERR_CRITICAL("Unable to allocate host storage\n");
	fprintf(stderr,"sizeof(person)=%i\n",(int) sizeof(person));
	for(i=0;i<P.NCP;i++)
		{
		j=(int) (CellLookup[i]-Cells);
		if(Cells[j].n>0)
			{
			if(!(Cells[j].InvCDF=(int *) malloc(1025*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
			if(!(Cells[j].max_trans=(float *) malloc(P.NCP*sizeof(float)))) ERR_CRITICAL("Unable to allocate cell storage\n");
			if(!(Cells[j].cum_trans=(float *) malloc(P.NCP*sizeof(float)))) ERR_CRITICAL("Unable to allocate cell storage\n");
			}
		}
	for(i=0;i<P.NC;i++)
		{
		Cells[i].cumTC=0;
		for(j=0;j<Cells[i].n;j++) Cells[i].members[j]=-1;
		}
	fprintf(stderr,"Cells assigned\n");
	for(i=0;i<=MAX_HOUSEHOLD_SIZE;i++) denom_household[i]=0;
	P.NH=0;
	for(i=j2=0;j2<P.NMCP;j2++)
		{
		j=(int) (McellLookup[j2]-Mcells);
		l=((j/P.nmch)/P.NMCL)*P.nch+((j%P.nmch)/P.NMCL);
		ad=((P.DoAdunitDemog)&&(P.DoAdUnits))?Mcells[j].adunit:0;
		Mcells[j].FirstHousehold = P.NH; //added first household so we can link households to microcells
		for(k=0;k<Mcells[j].n;)
			{
			m=1;
			if(P.DoHouseholds)
				{
				s=ranf_mt(0);
				while((s>P.HouseholdSizeDistrib[ad][m-1])&&(k+m<Mcells[j].n)&&(m<MAX_HOUSEHOLD_SIZE)) m++;
				}
			denom_household[m]++;
			//Case detection by household
			if (P.DoClusterVaccAccept||P.DoClusterCaseDetection)
			{
				temp_rep_rate = ranf();
				temp_vacc_accept = ranf();
			}
			for(i2=0;i2<m;i2++)
			{
//				fprintf(stderr,"%i ",i+i2);
				Hosts[i+i2].listpos=m; //used temporarily to store household size
				Mcells[j].members[k+i2]=i+i2;
				Cells[l].susceptible[Cells[l].cumTC]=i+i2;
				Cells[l].members[Cells[l].cumTC++]=i+i2;
				Hosts[i+i2].pcell=l;
				Hosts[i+i2].mcell=j;
				Hosts[i+i2].hh=P.NH;
				//add vaccination acceptance for each host here
				if (P.DoClusterVaccAccept)
				{
					Hosts[i + i2].vacc_accept = temp_vacc_accept;
				}
				else
				{
					Hosts[i + i2].vacc_accept = ranf();
				}
				if (P.DoClusterCaseDetection)
				{
					Hosts[i + i2].rep_rate = temp_rep_rate;
				}
				else
				{
					Hosts[i + i2].rep_rate = ranf();
				}
				//add case detection here: ggilani 13/06/19
				//if (P.DoCaseDetection)
				//{
				//	if (P.DoClusterCaseDetection)
				//	{
				//		Hosts[i + i2].rep_rate = temp_rep_rate;
				//	}
				//	else
				//	{
				//		Hosts[i + i2].rep_rate = ranf_mt(0);
				//	}
				//}
			}
			P.NH++;
			i+=m;
			k+=m;
			Mcells[j].nh++;
			}
		}
	if(!(Households=(household *) malloc(P.NH*sizeof(household)))) ERR_CRITICAL("Unable to allocate household storage\n");
	for(j=0;j<NUM_AGE_GROUPS;j++) AgeDist[j]=AgeDist2[j]=0;
	if(P.DoHouseholds) fprintf(stderr,"Household sizes assigned to %i people\n",i);
#pragma omp parallel for private(tn,j2,j,i,k,x,y,xh,yh,i2,m) schedule(static,1)
	for(tn=0;tn<P.NumThreads;tn++)
	  for(j2=tn;j2<P.NMCP;j2+=P.NumThreads)
		{
		j=(int) (McellLookup[j2]-Mcells);
		x=(double) (j/P.nmch);
		y=(double) (j%P.nmch);
		i=Mcells[j].members[0];
		if(j%100==0) 
			fprintf(stderr,"%i=%i (%i %i)            \r",j,Mcells[j].n,Mcells[j].adunit,(AdUnits[Mcells[j].adunit].id%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor);
		for(k=0;k<Mcells[j].n;)
			{
			m=Hosts[i].listpos;
			xh=P.mcwidth*(ranf_mt(tn)+x);
			yh=P.mcheight*(ranf_mt(tn)+y);
			AssignHouseholdAges(m,i,tn);
			for(i2=0;i2<m;i2++) Hosts[i+i2].listpos=0;
			if(P.DoHouseholds)
			{
				for(i2=0;i2<m;i2++){
					Hosts[i+i2].nc_plus_hh_disabled=0; //added this so that households are included
					Hosts[i+i2].inf=0; //added this so that infection status is set to zero and household r0 is correctly calculated
				}
			}
			Households[Hosts[i].hh].FirstPerson=i;
			Households[Hosts[i].hh].nh=m;
			Households[Hosts[i].hh].nhr=m;
			Households[Hosts[i].hh].loc_x=xh;
			Households[Hosts[i].hh].loc_y=yh;
			i+=m;
			k+=m;
			}
		}
	fprintf(stderr,"Ages/households assigned\n");


	//adding some code to call set up road networks, if including the roads network file: ggilani 12/02/15
	if(P.DoRoadNetwork)
	{
		SetupRoads();
		fprintf(stderr,"Roads assigned to microcells\n");
	}


	if(!P.DoRandomInitialInfectionLoc)
		{
		k=(int) (P.LocationInitialInfection[0][0]/P.mcwidth);
		l=(int) (P.LocationInitialInfection[0][1]/P.mcheight);
		j=k*P.nmch+l;
		
		double rand_r=0.0; //added these variables so that if initial infection location is empty we can search the 10km neighbourhood to find a suitable cell
		double rand_theta=0.0;
		int counter=0;
		if(Mcells[j].n<P.NumInitialInfections[0])
		{
			while(Mcells[j].n<P.NumInitialInfections[0]&&counter<100)
			{
				rand_r=ranf(); rand_theta=ranf(); 
				rand_r=0.083*sqrt(rand_r); rand_theta=2*PI*rand_theta; //rand_r is multiplied by 0.083 as this is roughly equal to 10km in decimal degrees
				k=(int) ((P.LocationInitialInfection[0][0]+rand_r*cos(rand_theta))/P.mcwidth);
				l=(int) ((P.LocationInitialInfection[0][1]+rand_r*sin(rand_theta))/P.mcheight);
				j=k*P.nmch+l;
				counter++;
			}
			if(counter<100)
			{
				P.LocationInitialInfection[0][0]=P.LocationInitialInfection[0][0]+rand_r*cos(rand_theta); //set LocationInitialInfection to actual one used
				P.LocationInitialInfection[0][1]=P.LocationInitialInfection[0][1]+rand_r*sin(rand_theta);
			}
		}
		if(Mcells[j].n<P.NumInitialInfections[0])
			ERR_CRITICAL("Too few people in seed microcell to start epidemic with required number of initial infectionz.\n");
		}
	fprintf(stderr,"Checking cells...\n");
	maxd=0;last_i=0;
	for(i=0;i<P.NMC;i++)
		{
		l=((i/P.nmch)/P.NMCL)*P.nch+((i%P.nmch)/P.NMCL);
		if(Cells[l].n==0)
			mcell_dens[i]=0;
		else
			last_i=i;
		maxd+=mcell_dens[i];
		}
	fprintf(stderr,"Allocating place/age groups...\n");
	for(k=0;k<NUM_AGE_GROUPS*AGE_GROUP_WIDTH;k++)
		{
		for(l=0;l<P.PlaceTypeNum;l++)
			{
			PropPlaces[k][l]=PropPlacesC[k][l]=0.0;
			if((k<P.PlaceTypeAgeMax[l])&&(k>=P.PlaceTypeAgeMin[l]))
				PropPlaces[k][l]+=P.PlaceTypePropAgeGroup[l];
			if((k<P.PlaceTypeAgeMax2[l])&&(k>=P.PlaceTypeAgeMin2[l]))
				PropPlaces[k][l]+=P.PlaceTypePropAgeGroup2[l];
			if((k<P.PlaceTypeAgeMax3[l])&&(k>=P.PlaceTypeAgeMin3[l]))
				PropPlaces[k][l]+=P.PlaceTypePropAgeGroup3[l];
			if(l==HOTEL_PLACE_TYPE)
				PropPlacesC[k][l]=((l>0)?PropPlacesC[k][l-1]:0);
			else
				PropPlacesC[k][l]=PropPlaces[k][l]+((l>0)?PropPlacesC[k][l-1]:0);
			}
		}
/*
	for(l=0;l<P.PlaceTypeNum;l++)
		{
		for(k=0;k<NUM_AGE_GROUPS*AGE_GROUP_WIDTH;k++)
			fprintf(stderr, "%i:%lg ",k,PropPlaces[k][l]);
		fprintf(stderr,"\n");
		}
*/
/*	if((P.DoAdUnits)&&(P.DoAdunitDemog))
		{for(i=0;i<P.NumAdunits;i++) free(State.InvAgeDist[i]);}
	else
		free(State.InvAgeDist[0]);
	free(State.InvAgeDist);
*/	P.nsp=0;
	if((P.DoSchoolFile)&&(P.DoPlaces))
		{
		fprintf(stderr,"Reading school file\n");
		if(!(dat=fopen(SchoolFile,"r"))) ERR_CRITICAL("Unable to open school file\n");
		fscanf(dat,"%i",&P.nsp);
		for(j=0;j<P.nsp;j++)
			{
			fscanf(dat,"%i %i",&m,&(P.PlaceTypeMaxAgeRead[j]));
			if(!(Places[j]=(place *) malloc(m*sizeof(place)))) ERR_CRITICAL("Unable to allocate place storage\n");
			for(i=0;i<m;i++)
				if(!(Places[j][i].AvailByAge =(unsigned short int *) malloc(P.PlaceTypeMaxAgeRead[j]*sizeof(unsigned short int)))) ERR_CRITICAL("Unable to allocate place storage\n");
			P.Nplace[j]=0;
			for(i=0;i<P.NMC;i++) Mcells[i].np[j]=0;
			}
		mr=0;
		while(!feof(dat))
			{
			fscanf(dat,"%lg %lg %i %i",&x,&y,&j,&m);
			for(i=0;i<P.PlaceTypeMaxAgeRead[j];i++) fscanf(dat,"%hu",&(Places[j][P.Nplace[j]].AvailByAge[i]));
			Places[j][P.Nplace[j]].loc_x=(float) (x-P.SpatialBoundingBox[0]);
			Places[j][P.Nplace[j]].loc_y=(float) (y-P.SpatialBoundingBox[1]);
			if((x>=P.SpatialBoundingBox[0])&&(x<P.SpatialBoundingBox[2])&&(y>=P.SpatialBoundingBox[1])&&(y<P.SpatialBoundingBox[3]))
				{
				i=P.nch*((int) (Places[j][P.Nplace[j]].loc_x/P.cwidth))+((int) (Places[j][P.Nplace[j]].loc_y/P.cheight));
				if(Cells[i].n==0) mr++;					
				Places[j][P.Nplace[j]].n=m;
				i=(int) (Places[j][P.Nplace[j]].loc_x/P.mcwidth);
				k=(int) (Places[j][P.Nplace[j]].loc_y/P.mcheight);
				j2=i*P.nmch+k;
				Mcells[j2].np[j]++;
				Places[j][P.Nplace[j]].mcell=j2;
				P.Nplace[j]++;
				if(P.Nplace[j]%1000==0) fprintf(stderr,"%i read    \r",P.Nplace[j]);
				}
			}
		fclose(dat);
		fprintf(stderr,"%i schools read (%i in empty cells)      \n",P.Nplace[j],mr);
		for(i=0;i<P.NMC;i++)
			for(j=0;j<P.nsp;j++)
				if(Mcells[i].np[j]>0)
					{
					if(!(Mcells[i].places[j]=(int *) malloc(Mcells[i].np[j]*sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
					Mcells[i].np[j]=0;
					}
		for(j=0;j<P.nsp;j++)
			{
			t=s=0;
			for(i=0;i<P.N;i++)
				t+=PropPlaces[HOST_AGE_YEAR(i)][j];
			for(i=0;i<P.Nplace[j];i++)
				{
				k=Places[j][i].mcell;
				Mcells[k].places[j][Mcells[k].np[j]++]=i;
				s+=(double) Places[j][i].n;
				}
			fprintf(stderr,"School type %i: capacity=%lg demand=%lg\n",j,s,t);
			t/=s;
			for(i=0;i<P.Nplace[j];i++)
				Places[j][i].n=(int) ceil(((double) Places[j][i].n)*t);
			}
		}
	if(P.DoPlaces)
		{
		fprintf(stderr,"Configuring places...\n");
#pragma omp parallel for private(tn,j2,i,j,k,t,m,s,x,y,xh,yh) schedule(static,1)
		for(tn=0;tn<P.NumThreads;tn++)
		  for(j2=P.nsp+tn;j2<P.PlaceTypeNum;j2+=P.NumThreads)
			{
			t=0;
			P.PlaceTypeMaxAgeRead[j2]=0;
			for(i=0;i<P.N;i++)
				t+=PropPlaces[HOST_AGE_YEAR(i)][j2];
			P.Nplace[j2]=(int) ceil(t/P.PlaceTypeMeanSize[j2]);
			fprintf(stderr,"[%i:%i %g] ",j2,P.Nplace[j2],t);
			if(!(Places[j2]=(place *) malloc(P.Nplace[j2]*sizeof(place)))) ERR_CRITICAL("Unable to allocate place storage\n");
			t=1.0;
			for(m=i=k=0;i<P.NMC;i++)
				{
				s=mcell_dens[i]/maxd/t;
				if(s>1.0) s=1.0;
				if(i==last_i)
					m+=(Mcells[last_i].np[j2]=P.Nplace[j2]-m);
				else
					m+=(Mcells[i].np[j2]=(int) ignbin_mt((long) (P.Nplace[j2]-m),s,tn));
				t-=mcell_dens[i]/maxd;
				if(Mcells[i].np[j2]>0)
					{
					if(!(Mcells[i].places[j2]=(int *) malloc(Mcells[i].np[j2]*sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
					x=(double) (i/P.nmch);
					y=(double) (i%P.nmch);
					for(j=0;j<Mcells[i].np[j2];j++)
						{
						s=ranf_mt(tn);
						xh=P.mcwidth*(ranf_mt(tn)+x);
						yh=P.mcheight*(ranf_mt(tn)+y);
						Places[j2][k].loc_x=xh;
						Places[j2][k].loc_y=yh;
						Places[j2][k].n=0;
						Places[j2][k].mcell=i;
						Places[j2][k].country=Mcells[i].country;
						Mcells[i].places[j2][j]=k;
						k++;
						}
					}
				}
			}
		for(k=0;k<NUM_AGE_GROUPS*AGE_GROUP_WIDTH;k++)
			for(l=1;l<P.PlaceTypeNum;l++)
				if(l!=HOTEL_PLACE_TYPE)
					{
					if(PropPlacesC[k][l-1]<1)
						PropPlaces[k][l]/=(1-PropPlacesC[k][l-1]);
					else if (PropPlaces[k][l]!=0)
						PropPlaces[k][l]=1.0;
					}
		fprintf(stderr,"Places assigned\n");
		}
	free(mcell_dens);
	l=0;
	for(j=0;j<P.NC;j++)
		if(l<Cells[j].n) l=Cells[j].n;
	if(!(SamplingQueue=(int **) malloc(P.NumThreads*sizeof(int *)))) ERR_CRITICAL("Unable to allocate state storage\n");
	P.InfQueuePeakLength=P.N/P.NumThreads/INF_QUEUE_SCALE;
#pragma omp parallel for private(i,k) schedule(static,1)
	for(i=0;i<P.NumThreads;i++)
		{
		if(!(SamplingQueue[i]=(int *) malloc(2*(MAX_PLACE_SIZE+CACHE_LINE_SIZE)*sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
		for(k=0;k<P.NumThreads;k++)
			if(!(StateT[i].inf_queue[k]=(int *) malloc(P.InfQueuePeakLength*sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
		if(!(StateT[i].cell_inf=(float *) malloc((l+1)*sizeof(float)))) ERR_CRITICAL("Unable to allocate state storage\n");
		if(!(StateT[i].inv_cell_inf=(int *) malloc(1025*sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
		}

	if((P.DoHospitalisation)&&(P.DoETUByAdUnit))
	{
		for(i=0;i<P.NumAdunits;i++)
		{
			if(!(AdUnits[i].h_queue=(int *) malloc(P.InfQueuePeakLength*sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			for(j=0;j<P.NumThreads;j++)
			{
				if(!(StateT[j].h_queue[i]=(int *) malloc(P.InfQueuePeakLength*sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
				if(!(StateT[j].hd_queue[i]=(int *) malloc(P.InfQueuePeakLength*sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			}
		}
	}


	// set up queues for contact tracing, similar to queues for hospitalisation above: ggilani 12/06/17
//	if((P.DoNewContactTracing)&&(P.DoAdUnits))
	if (P.DoAdUnits)
	{
		for(i=0;i<P.NumAdunits;i++)
		{
			if(!(AdUnits[i].ct_queue=(int *) malloc(P.InfQueuePeakLength*sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			if(!(AdUnits[i].ct=(int *) malloc(P.InfQueuePeakLength*sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			for(j=0;j<P.NumThreads;j++)
			{
				if(!(StateT[j].ct_queue[i]=(int *) malloc(P.InfQueuePeakLength*sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			}
		}
	}

	//set up queues for vaccination (both ring and geographic)
	if ((P.DoRingVaccination)||(P.DoGeoVaccination))
	{
		for (j = 0; j < P.NumThreads; j++)
		{
			if (!(StateT[j].vacc_queue = (int*)malloc(P.N * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			if (!(StateT[j].ring_queue = (int*)malloc(P.N * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n"); // to store the ring a person is in - potentially useful for dosage
			
		}
		if (!(State.vacc_queue = (int*)malloc(P.N * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
		if (!(State.ring_queue = (int*)malloc(P.N * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n"); // to store the ring a person is in - potentially useful for dosage

	}
	fprintf(stderr, "Population: %i\n",P.N);
	//if (P.DoGeoVaccination)
	//{
	//	for (j = 0; j < P.NumThreads; j++)
	//	{
	//		if (!(StateT[j].geovacc_queue = (int*)malloc(P.N * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
	//	}
	//	if (!(State.geovacc_queue = (int*)malloc(P.N * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
	//}

	//add storage to store each ring, for both contact tracing and ring vaccination - ggilani 21/10/19
	if ((P.DoRingVaccination) || (P.DoContactTracing))
	{
		for (j = 0; j < P.NumThreads; j++)
		{
			if (!(StateT[j].ringvacclist = (int*)malloc(MAX_RING_SIZE * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			if (!(StateT[j].ringlist = (int*)malloc(MAX_RING_SIZE * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
		}
	}

	//If outputting origin-destination matrix, set up storage for flow between admin units
	if((P.DoAdUnits)&&(P.DoOriginDestinationMatrix))
	{
		for(i=0;i<P.NumAdunits;i++)
		{
			if(!(AdUnits[i].origin_dest=(double *) malloc(MAX_ADUNITS*sizeof(double)))) ERR_CRITICAL("Unable to allocate storage for origin destination matrix\n");
			for(j=0;j<P.NumThreads;j++)
			{
				if(!(StateT[j].origin_dest[i]=(double *) calloc(MAX_ADUNITS,sizeof(double)))) ERR_CRITICAL("Unable to allocate state origin destination matrix storage\n");
			}
			//initialise to zero
			for(j=0;j<P.NumAdunits;j++)
			{
				AdUnits[i].origin_dest[j]=0.0;
			}
		}
	}

	for(i=0;i<P.NC;i++)
		{
		Cells[i].cumTC=0;
		Cells[i].S=Cells[i].n;
		Cells[i].L=Cells[i].I=0;
		}
	fprintf(stderr,"Allocated cell and host memory\n");
	fprintf(stderr,"Assigned hosts to cells\n");

	if(!(income_distrib=(double *) malloc(10000*sizeof(double)))) ERR_CRITICAL("Unable to allocate storage\n");

	if(P.DoAdunitBoundaries)
		{
		for(i=0;i<P.NumAdunits;i++) P.HouseholdsByAdunit[i]=P.PrivateStockByAdunit[i]=0;
		for(i=0;i<P.NH;i++)
			{
			j=Mcells[Hosts[Households[i].FirstPerson].mcell].adunit;
			Households[i].income=P.MedianIncomeByAdunit[j]*pow(-1.442695041*log(ranf()),1/P.IncomeWeibullPower);
			if(Households[i].income>=1000000)
				income_distrib[9999]++;
			else
				income_distrib[(int) (Households[i].income/100)]++;
			P.HouseholdsByAdunit[j]++;
			}
		}
	else
		{
		for(i=0;i<P.NH;i++)
			{
			Households[i].income=P.MedianIncome*pow(-1.442695041*log(ranf()),1/P.IncomeWeibullPower);
			if(Households[i].income>=1000000)
				income_distrib[9999]++;
			else
				income_distrib[(int) (Households[i].income/100)]++;
			}
		}
	for(i=1;i<10000;i++)
		income_distrib[i]+=income_distrib[i-1];
	for(i=0;i<9999;i++)
		income_distrib[i]/=income_distrib[9999];
	income_distrib[9999]=1;
	fprintf(stderr,"Assigned household incomes\n");

	if(P.PropPrivateStockpile==0)
		{
		k=0;
		for(i=0;i<P.NH;i++)
			Households[i].stockpile=0;
		}
	else
		{
		if(P.PrivateStockpileOrderByIncome==1)
			{
			for(j=9999;(j>0)&&((1-income_distrib[j-1])<=P.PropPrivateStockpile);) j--;
			t=(j<9999)?(100.0*((double) (j+1))):1e20;
			k=0;
			for(i=0;i<P.NH;i++)
				k+=(Households[i].stockpile=(Households[i].income>t)?1:0);
			}
		else if(P.PrivateStockpileOrderByIncome==-1)
			{
			for(j=0;(j<10000)&&(income_distrib[j]<=P.PropPrivateStockpile);) j++;
			t=(j<9999)?(100.0*((double) j)):1e20;
			k=0;
			for(i=0;i<P.NH;i++)
				k+=(Households[i].stockpile=(Households[i].income<t)?1:0);
			}
		else
			{
			k=0;
			for(i=0;i<P.NH;i++)
				k+=(Households[i].stockpile=((ranf()<P.PropPrivateStockpile)?1:0));
			}
		if(P.DoAdunitBoundaries)
			{
			for(i=0;i<P.NH;i++)
				if(Households[i].stockpile!=0) P.PrivateStockByAdunit[Mcells[Hosts[Households[i].FirstPerson].mcell].adunit]++;
			}
		fprintf(stderr,"Proportion of households with private antiviral stockpiles = %lg\n",((double) k)/((double) P.NH));
		}
	free(income_distrib);

}


void SetupAirports(void)
{
	int i,j,k,l,m;
	double x,y,t,tmin;
	indexlist *base,*cur;

	fprintf(stderr,"Assigning airports to microcells\n");
	P.KernelType=P.AirportKernelType;
	P.KernelScale=P.AirportKernelScale;
	P.KernelShape=P.AirportKernelShape;
	P.KernelP3=P.AirportKernelP3;
	P.KernelP4=P.AirportKernelP4;
	InitKernel(1,1.0);
	if(!(Airports[0].DestMcells =(indexlist *) calloc(P.NMCP*NNA,sizeof(indexlist)))) ERR_CRITICAL("Unable to allocate airport storage\n");
	if(!(base =(indexlist *) calloc(P.NMCP*NNA,sizeof(indexlist)))) ERR_CRITICAL("Unable to allocate airport storage\n");
	for(i=0;i<P.Nairports;i++) Airports[i].num_mcell=0;
	cur=base;
	for(i=0;i<P.NMC;i++)
		if(Mcells[i].n>0)
			{
			Mcells[i].AirportList=cur;
			cur+=NNA;
			}
#pragma omp parallel for private(i,j,k,l,x,y,t,tmin) schedule(static,10000)
	for(i=0;i<P.NMC;i++)
		if(Mcells[i].n>0)
			{
			if(i%10000==0) fprintf(stderr,"\n%i           ",i);
			x=(((double) (i/P.nmch))+0.5)*P.mcwidth;
			y=(((double) (i%P.nmch))+0.5)*P.mcheight;
			k=l=0;
			tmin=1e20;
			for(j=0;j<P.Nairports;j++)
				if(Airports[j].total_traffic>0)
					{
					t=numKernel(dist2_raw(x,y,Airports[j].loc_x,Airports[j].loc_y))*Airports[j].total_traffic;
					if(k<NNA)
						{
						Mcells[i].AirportList[k].id=j;
						Mcells[i].AirportList[k].prob=t;
						if(t<tmin) {tmin=t;l=k;}
						k++;
						}
					else if(t>tmin)
						{
						Mcells[i].AirportList[l].id=j;
						Mcells[i].AirportList[l].prob=t;
						tmin=1e20;
						for(k=0;k<NNA;k++)
							if(Mcells[i].AirportList[k].prob<tmin)
								{
								tmin=Mcells[i].AirportList[k].prob;
								l=k;
								}
						}
					}
			for(j=0;j<NNA;j++)
				Airports[Mcells[i].AirportList[j].id].num_mcell++;
			}
	cur=Airports[0].DestMcells;
	fprintf(stderr,"Microcell airport lists collated.\n");
	for(i=0;i<P.Nairports;i++)
		{
		Airports[i].DestMcells=cur;
		cur+=Airports[i].num_mcell;
		Airports[i].num_mcell=0;
		}
#pragma omp parallel for private(i,j,k,l,t,tmin) schedule(static,10000)
	for(i=0;i<P.NMC;i++)
		if(Mcells[i].n>0)
			{
			if(i%10000==0) fprintf(stderr,"\n%i           ",i);
			t=0;
			for(j=0;j<NNA;j++)
				{
				t+=Mcells[i].AirportList[j].prob;
				k=Mcells[i].AirportList[j].id;
#pragma omp critical (airport)
				l=(Airports[k].num_mcell++);
				Airports[k].DestMcells[l].id=i;
				Airports[k].DestMcells[l].prob=Mcells[i].AirportList[j].prob*((float) Mcells[i].n);
				}
			tmin=0;
			for(j=0;j<NNA;j++)
				{
				Mcells[i].AirportList[j].prob=tmin+Mcells[i].AirportList[j].prob/t;
				tmin=Mcells[i].AirportList[j].prob;
				}
			}
	fprintf(stderr,"Airport microcell lists collated.\n");
	for(i=0;i<P.Nairports;i++)
		if(Airports[i].total_traffic>0)
			{
			for(j=1;j<Airports[i].num_mcell;j++)
				Airports[i].DestMcells[j].prob+=Airports[i].DestMcells[j-1].prob;
			t=Airports[i].DestMcells[Airports[i].num_mcell-1].prob;
			if(t==0) t=1.0;
			for(j=0;j<Airports[i].num_mcell-1;j++)
				Airports[i].DestMcells[j].prob/=t;
			if(Airports[i].num_mcell>0) Airports[i].DestMcells[Airports[i].num_mcell-1].prob=1.0;
			for(j=l=0;l<=1024;l++)
				{
				t=((double) l)/1024.0;
				while(Airports[i].DestMcells[j].prob<t) j++;
				Airports[i].Inv_DestMcells[l]=j;
				}
			l=0;
			for(j=0;j<Airports[i].num_mcell;j++)
				l+=Mcells[Airports[i].DestMcells[j].id].np[HOTEL_PLACE_TYPE];
			if(l<10)
				{
				fprintf(stderr,"(%i ",l);
				l=0;
				for(j=0;j<Airports[i].num_mcell;j++)
					l+=Mcells[Airports[i].DestMcells[j].id].n;
				fprintf(stderr,"%i %i) ",Airports[i].num_mcell,l);
				}
			}
	fprintf(stderr,"\nInitialising hotel to airport lookup tables\n");
	free(base);
#pragma omp parallel for private(i,j,l,m,t,tmin) schedule(static,1)
	for(i=0;i<P.Nairports;i++)
		if(Airports[i].total_traffic>0)
			{
			m=(int) (Airports[i].total_traffic/HOTELS_PER_1000PASSENGER/1000);
			if(m<MIN_HOTELS_PER_AIRPORT) m=MIN_HOTELS_PER_AIRPORT;
			fprintf(stderr,"\n%i    ",i);
			tmin=MAX_DIST_AIRPORT_TO_HOTEL*MAX_DIST_AIRPORT_TO_HOTEL*0.75;
			do 
				{
				tmin+=0.25*MAX_DIST_AIRPORT_TO_HOTEL*MAX_DIST_AIRPORT_TO_HOTEL;
				Airports[i].num_place=0;
				for(j=0;j<P.Nplace[HOTEL_PLACE_TYPE];j++)
					if(dist2_raw(Airports[i].loc_x,Airports[i].loc_y,
						Places[HOTEL_PLACE_TYPE][j].loc_x,Places[HOTEL_PLACE_TYPE][j].loc_y)<tmin)
						Airports[i].num_place++;	
				} 
			while(Airports[i].num_place<m);
			if(tmin>MAX_DIST_AIRPORT_TO_HOTEL*MAX_DIST_AIRPORT_TO_HOTEL) fprintf(stderr,"*** %i : %lg %i ***\n",i,sqrt(tmin),Airports[i].num_place);
			if(!(Airports[i].DestPlaces =(indexlist *) calloc(Airports[i].num_place,sizeof(indexlist)))) ERR_CRITICAL("Unable to allocate airport storage\n");
			Airports[i].num_place=0;
			for(j=0;j<P.Nplace[HOTEL_PLACE_TYPE];j++)
				if((t=dist2_raw(Airports[i].loc_x,Airports[i].loc_y,
					Places[HOTEL_PLACE_TYPE][j].loc_x,Places[HOTEL_PLACE_TYPE][j].loc_y))<tmin)
					{
					Airports[i].DestPlaces[Airports[i].num_place].prob=numKernel(t);
					Airports[i].DestPlaces[Airports[i].num_place].id=j;
					Airports[i].num_place++;
					}
			t=0;
			for(j=0;j<Airports[i].num_place;j++)
				{
				Airports[i].DestPlaces[j].prob=t+Airports[i].DestPlaces[j].prob;
				t=Airports[i].DestPlaces[j].prob;
				}
			for(j=0;j<Airports[i].num_place-1;j++)
				Airports[i].DestPlaces[j].prob/=t;
			if(Airports[i].num_place>0) Airports[i].DestPlaces[Airports[i].num_place-1].prob=1.0;
			for(j=l=0;l<=1024;l++)
				{
				t=((double) l)/1024.0;
				while(Airports[i].DestPlaces[j].prob<t) j++;
				Airports[i].Inv_DestPlaces[l]=j;
				}
			}
	P.KernelType=P.MoveKernelType;
	P.KernelScale=P.MoveKernelScale;
	P.KernelShape=P.MoveKernelShape;
	P.KernelP3=P.MoveKernelP3;
	P.KernelP4=P.MoveKernelP4;
	for(i=0;i<P.Nplace[HOTEL_PLACE_TYPE];i++) Places[HOTEL_PLACE_TYPE][i].n=0;
	InitKernel(0,1.0);
	fprintf(stderr,"\nAirport initialisation completed successfully\n");
}


void InitKernel(int DoPlaces,double norm)
{
	int i,j,l,m,im,mcl_from,mcl_to,mcl_ind,k,p,mc,mc_ind;
	double t,t2,d;
	int same_country;
	//Adding some road accessibility parameters for scaling distances if the road network is used
	double accessFrom,accessTo;

	if(P.KernelType==1)
		Kernel=ExpKernel;
	else if(P.KernelType==2)
		Kernel=PowerKernel;
	else if(P.KernelType==3)
		Kernel=GaussianKernel;
	else if(P.KernelType==4)
		Kernel=StepKernel;
	else if(P.KernelType==5)
		Kernel=PowerKernelB;
	else if(P.KernelType==6)
		Kernel=PowerKernelUS;
	else if(P.KernelType==7)
		Kernel=PowerExpKernel;
	t2=0;
#pragma omp parallel for private(i) schedule(static,500) //added private i
	for(i=0;i<=NKR;i++)
		{
		nKernel[i]=(*Kernel)(((double) i)*P.KernelDelta)/norm;
		nKernelHR[i]=(*Kernel)(((double) i)*P.KernelDelta/NK_HR)/norm;
		}
#pragma omp parallel for schedule(static,500) private(i,j,l,m,im,same_country,mc,mc_ind,k,p)
	for(i=0;i<P.NCP;i++)
		{
		l=CellLookup[i]-Cells;
		Cells[l].tot_prob=0;
		//accessFrom=0.0;//1.0;
		/*if(P.DoRoadNetwork)
		{
			mcl_ind=(l/P.nch)*P.NMCL*P.nmch+(l%P.nch)*P.NMCL;
			//loop over microcells in these cells to find populations in each admin unit and so flows
			for(k=0;k<P.NMCL;k++)
			{
				for(p=0;p<P.NMCL;p++)
				{
					//get index of microcell
					mcl_from=mcl_ind+p+k*P.nmch;
					if((Mcells[mcl_from].road_connection>0)&&(Mcells[mcl_from].road_connection<=P.MaxRoadType))
					{
						accessFrom=1;//P.RoadAccessDistance;
					}
				}
			}
		}*/
		for(j=im=0;j<P.NCP;j++)
			{
			m=CellLookup[j]-Cells;
			//accessTo=0.0;//1.0;
			/*if(P.DoRoadNetwork)
			{
				mcl_ind=(m/P.nch)*P.NMCL*P.nmch+(m%P.nch)*P.NMCL;
				//loop over microcells in these cells to find populations in each admin unit and so flows
				for(k=0;k<P.NMCL;k++)
				{
					for(p=0;p<P.NMCL;p++)
					{
						//get index of microcell
						mcl_to=mcl_ind+p+k*P.nmch;
						if((Mcells[mcl_to].road_connection>0)&&(Mcells[mcl_to].road_connection<=P.MaxRoadType))
						{
							accessTo=1;//P.RoadAccessDistance;
						}
					}
				}
			}*/
			if((P.DoRoadNetwork)&&(P.DoRoadDistanceEffect)&&((Cells[m].road_connection)||(Cells[l].road_connection)))
				{
					Cells[l].tot_prob+=(Cells[l].max_trans[j]=numKernel(dist2_cc_min(Cells+l,Cells+m))+Cells[l].road_access+Cells[m].road_access)*Cells[m].n; //scaling distance according to accessibility by road - was accessTo*accessFrom //Cells[l].road_access*Cells[m].road_access numKernel(P.RoadAccessDistance*dist2_cc_min(Cells+l,Cells+m)))*Cells[m].n; numKernel(dist2_cc_min(Cells+l,Cells+m))+P.RoadAccessDistance)*Cells[m].n;
				}
			else if((P.DoCapitalCityEffect)&&(P.DoCapitalCityDistanceEffect)&&(P.DoAdUnits)&&((Cells[l].capital_city)^(Cells[m].capital_city))) //if we're doing the capital city effect and either cell l or cell m (but not both) are in capital cities
			{
				same_country=0;
				if(Cells[l].capital_city)
				{
					//if cell m contains a capital city district, check to see whether any of cell m is in the same country
					mc=(m/P.nch)*P.NMCL*P.nmch+(m%P.nch)*P.NMCL;
					for(k=0;k<P.NMCL;k++)
					{
						for(p=0;p<P.NMCL;p++)
						{
							//get index of microcell
							mc_ind=mc+p+k*P.nmch;
							if((AdUnits[Mcells[mc_ind].adunit].id/P.CountryDivisor)==(Cells[l].capital_city/P.CountryDivisor))
							{
								same_country=1;//P.RoadAccessDistance;
							}
						}
					}
				}
				else if (Cells[m].capital_city)
				{
					//if cell m contains a capital city district, check to see whether any of cell l is in the same country
					mc=(l/P.nch)*P.NMCL*P.nmch+(l%P.nch)*P.NMCL;
					for(k=0;k<P.NMCL;k++)
					{
						for(p=0;p<P.NMCL;p++)
						{
							//get index of microcell
							mc_ind=mc+p+k*P.nmch;
							if((AdUnits[Mcells[mc_ind].adunit].id/P.CountryDivisor)==(Cells[m].capital_city/P.CountryDivisor))
							{
								same_country=1;//P.RoadAccessDistance;
							}
						}
					}
				}
				if(same_country)
				{
					Cells[l].tot_prob+=(Cells[l].max_trans[j]=numKernel(P.CapitalCityDistanceEffect*dist2_cc_min(Cells+l,Cells+m)))*Cells[m].n;  //numKernel(P.CapitalCityDistanceEffect*dist2_cc_min(Cells+l,Cells+m)))*Cells[m].n; numKernel(dist2_cc_min(Cells+l,Cells+m))+P.CapitalCityDistanceEffect)*Cells[m].n;
				}
				else
				{
					Cells[l].tot_prob+=(Cells[l].max_trans[j]=numKernel(dist2_cc_min(Cells+l,Cells+m)))*Cells[m].n;
				}
				}
			else if((P.DoCapitalCityEffect)&&(P.DoCapitalCityAddEffect)&&(P.DoAdUnits)&&((Cells[l].capital_city)||(Cells[m].capital_city))) //if we're doing the capital city effect and either cell l or cell m are in capital cities
			{
				same_country=0;
				if(Cells[l].capital_city)
				{
					//if cell m contains a capital city district, check to see whether any of cell m is in the same country
					mc=(m/P.nch)*P.NMCL*P.nmch+(m%P.nch)*P.NMCL;
					for(k=0;k<P.NMCL;k++)
					{
						for(p=0;p<P.NMCL;p++)
						{
							//get index of microcell
							mc_ind=mc+p+k*P.nmch;
							if((AdUnits[Mcells[mc_ind].adunit].id/P.CountryDivisor)==(Cells[l].capital_city/P.CountryDivisor))
							{
								same_country=1;
							}
						}
					}
				}
				else if (Cells[m].capital_city)
				{
					//if cell m contains a capital city district, check to see whether any of cell l is in the same country
					mc=(l/P.nch)*P.NMCL*P.nmch+(l%P.nch)*P.NMCL;
					for(k=0;k<P.NMCL;k++)
					{
						for(p=0;p<P.NMCL;p++)
						{
							//get index of microcell
							mc_ind=mc+p+k*P.nmch;
							if((AdUnits[Mcells[mc_ind].adunit].id/P.CountryDivisor)==(Cells[m].capital_city/P.CountryDivisor))
							{
								same_country=1;
							}
						}
					}
				}
				if(same_country)
				{
					Cells[l].tot_prob+=(Cells[l].max_trans[j]=numKernel(dist2_cc_min(Cells+l,Cells+m))+P.CapitalCityAddEffect)*Cells[m].n;  //numKernel(P.CapitalCityDistanceEffect*dist2_cc_min(Cells+l,Cells+m)))*Cells[m].n; numKernel(dist2_cc_min(Cells+l,Cells+m))+P.CapitalCityDistanceEffect)*Cells[m].n;
				}
				else
				{
					Cells[l].tot_prob+=(Cells[l].max_trans[j]=numKernel(dist2_cc_min(Cells+l,Cells+m)))*Cells[m].n;
				}
			}
			else
				{
					Cells[l].tot_prob+=(Cells[l].max_trans[j]=numKernel(dist2_cc_min(Cells+l,Cells+m)))*Cells[m].n;
				}
			}
		}
}


/** function: SetupRoads(void)
 *
 * purpose: to read in the road network file and assign all microcells which contain a road a number which corresponds to the road type
 *
 * parameters: none
 *
 * returns: none
 *
 * author: ggilani, 12/02/15
 */
void SetupRoads(void)
{
	FILE *dat;
	int nRoadCells,i,road,j,k,l,mc_ind,c_ind,cell_x,cell_y,ind;
	double x,y;
	int *RoadType;
	double *X_Coords, *Y_Coords;
	double maxDist=sqrt(2*((double) P.MaxRoadNeighbour+1)*((double) P.MaxRoadNeighbour+1)); //maximum neighbourhood area to explore, and max distance (in cells) before road accessibility has no effect
	double accessGrad; //gradient in accessibility as cells get further from the road
	double dist,access;

	if(P.DoRoadDistanceEffect)
	{
		accessGrad=(0-P.RoadAccessDistance)/maxDist;
	}
	else if(P.DoRoadPopEffect)
	{
		accessGrad=(1-P.RoadAccessPop)/maxDist;
	}

	if(!(dat=fopen(RoadNetworkFile,"r"))) ERR_CRITICAL("Unable to open road network file!\n");
	i=0;
	//allocate temporary storage to hold coordinates and road types: can't be more than the number of microcells
	if(!(X_Coords=(double *) malloc(P.NMC*sizeof(double)))) ERR_CRITICAL("Unable to allocate storage for road coordinates\n");
	if(!(Y_Coords=(double *) malloc(P.NMC*sizeof(double)))) ERR_CRITICAL("Unable to allocate storage for road coordinates\n");
	if(!(RoadType=(int *) malloc(P.NMC*sizeof(int)))) ERR_CRITICAL("Unable to allocate storage for road types\n");

	while(!feof(dat))
	{
		fscanf(dat,"%lg %lg %i",&x,&y,&road);

		X_Coords[i]=x;
		Y_Coords[i]=y;
		RoadType[i]=road; 
		i++;

	}
	nRoadCells=i;
	fclose(dat);

	//go through cells first to set up initial states
	for(i=0;i<P.NC;i++)
	{
		Cells[i].road_connection=0;
		Cells[i].road_access=0.0;
	}

	//now go through list of road coordinates, find the associated microcell and cell and assign it a road type
	for(i=0;i<nRoadCells;i++)
	{
		//check to see whether segment of road is with the simulation area, and if it's within the range of road types we are considering
		if((X_Coords[i]>=P.SpatialBoundingBox[0])&&(Y_Coords[i]>=P.SpatialBoundingBox[1])&&(X_Coords[i]<P.SpatialBoundingBox[2])&&(Y_Coords[i]<P.SpatialBoundingBox[3])&&(RoadType[i]>0)&&(RoadType[i]<=P.MaxRoadType))
		{
			j=(int)floor((X_Coords[i]-P.SpatialBoundingBox[0])/P.mcwidth);
			k=(int)floor((Y_Coords[i]-P.SpatialBoundingBox[1])/P.mcheight);
			//actual index of microcell
			mc_ind=j*P.nmch+k;

			if((AdUnits[Mcells[mc_ind].adunit].id/P.CountryDivisor)==P.TargetCountry)
			{

			//mark microcell as containing a road if it doesn't contain a road already
			Mcells[mc_ind].road_connection=1;

			cell_x=(int)floor((X_Coords[i]-P.SpatialBoundingBox[0])/P.cwidth);
			cell_y=(int)floor((Y_Coords[i]-P.SpatialBoundingBox[1])/P.cheight);
			//actual index of cell
			c_ind=cell_x*P.nch+cell_y;

			//now look as cells instead of microcells: marking cells as 1 if they directly contain a road, 2 if they are in the neighbourhood of a road. If a cell is already marked
			//with a 1, we've already checked it and it's neighbourhood, so we can ignore it. But if it previously didn't have a road in it or was adjacent to a road, we want to update it
			if(Cells[c_ind].road_connection!=1)
			{
				//mark microcell as having having a road
				Cells[c_ind].road_connection=1;
				//set road accessibility to set value
				if(P.DoRoadDistanceEffect)
				{
					Cells[c_ind].road_access=P.RoadAccessDistance;
				}
				else if(P.DoRoadPopEffect)
				{
					Cells[c_ind].road_access=P.RoadAccessPop;
				}

				//Now we want to look in the neighbourhood of the cells
				for(j=-P.MaxRoadNeighbour;j<=P.MaxRoadNeighbour;j++)
				{
					for(k=-P.MaxRoadNeighbour;k<=P.MaxRoadNeighbour;k++)
					{
						//we've already set values for the actual cell itself
						if((j!=0)||(k!=0))
						{
							//get index of new cell
							ind=(cell_x+j)*P.nch+(cell_y+k);
							//we only want to update cells that don't have roads going right through them
							if(Cells[ind].road_connection==0)
							{
								Cells[ind].road_connection=2; //set road connection
								//find distance (in terms of cells) away from main cell
								dist=sqrt((double) (j*j+k*k));
								//set road accessibility based on this
								if(P.DoRoadDistanceEffect)
								{
									Cells[ind].road_access=accessGrad*dist+P.RoadAccessDistance;
								}
								else if(P.DoRoadPopEffect)
								{
									Cells[ind].road_access=accessGrad*dist+P.RoadAccessPop;
								}
							}
							else if(Cells[ind].road_connection==2)
							{
								//find distance from cell and accessibility
								dist=sqrt((double) (j*j+k*k));
								//if accessibility is lower (i.e. it is now closer to the road), update accessibility
								if(P.DoRoadDistanceEffect)
								{
									access=accessGrad*dist+P.RoadAccessDistance;
								}
								else if(P.DoRoadPopEffect)
								{
									access=accessGrad*dist+P.RoadAccessPop;
								}
								
								if(access<Cells[ind].road_access)
								{
									Cells[ind].road_access=access;
								}
							}
						}
					}
				}
			}

			} //cell matches target country

		}
	}

	//free memory used to temporarily store road properties
	free(X_Coords);
	free(Y_Coords);
	free(RoadType);

}

/**
 * function: DetermineCellsWithCapitalCities
 *
 * parameters: none
 * returns: none
 * function: to assign cells that overlap with admin units containing capital cities a capital city flag
 * author: ggilani, 26/02/15
 */
void DetermineCellsWithCapitalCities(void)
{

	int i,j,k,mc,mc_ind;

	for(i=0;i<P.NC;i++)
	{
		Cells[i].capital_city=0;
		//loop through microcells to find out what admin units overlap with the cell
		mc=(i/P.nch)*P.NMCL*P.nmch+(i%P.nch)*P.NMCL;
		for(j=0;j<P.NMCL;j++)
		{
			for(k=0;k<P.NMCL;k++)
			{
				mc_ind=mc+k+j*P.nmch;
				//test against capital city adunits
				if(AdUnits[Mcells[mc_ind].adunit].id==P.CapitalCityAdunit)
				{
					if(Cells[i].capital_city==0)
					{
						Cells[i].capital_city=P.CapitalCityAdunit;
					}
				}
				else if(AdUnits[Mcells[mc_ind].adunit].id==P.CapitalCityAdunit2)
				{
					if(Cells[i].capital_city==0)
					{
						Cells[i].capital_city=P.CapitalCityAdunit2;
					}
				}
				else if(AdUnits[Mcells[mc_ind].adunit].id==P.CapitalCityAdunit3)
				{
					if(Cells[i].capital_city==0)
					{
						Cells[i].capital_city=P.CapitalCityAdunit3;
					}
				}

			}
		}
	}
}

/**
 * function: initRadiationModel(void)
 *
 * purpose: to initialise a radiation mobility model for spatial transmission
 *
 *
 */
//void initRadiationModel(void)
//{
//	//permutations of moving around cells
//	const int ord[24][4][2]={	{{1,1},{-1,-1},{1,-1},{-1,1}},
//								{{-1,-1},{1,1},{1,-1},{-1,1}},
//								{{1,1},{1,-1},{-1,1},{-1,-1}},
//								{{1,1},{1,-1},{-1,-1},{-1,1}},
//								{{1,1},{-1,1},{1,-1},{-1,-1}},
//								{{1,1},{-1,1},{-1,-1},{1,-1}},
//								{{1,1},{-1,-1},{-1,1},{1,-1}},
//								{{1,-1},{1,1},{-1,1},{-1,-1}},
//								{{1,-1},{1,1},{-1,-1},{-1,1}},
//								{{1,-1},{-1,1},{1,1},{-1,-1}},
//								{{1,-1},{-1,1},{-1,-1},{1,1}},
//								{{1,-1},{-1,-1},{1,1},{-1,1}},
//								{{1,-1},{-1,-1},{-1,1},{1,1}},
//								{{-1,1},{1,-1},{1,1},{-1,-1}},
//								{{-1,1},{1,-1},{-1,-1},{1,1}},
//								{{-1,1},{1,1},{1,-1},{-1,-1}},
//								{{-1,1},{1,1},{-1,-1},{1,-1}},
//								{{-1,1},{-1,-1},{1,-1},{1,1}},
//								{{-1,1},{-1,-1},{1,1},{1,-1}},
//								{{-1,-1},{1,-1},{-1,1},{1,1}},
//								{{-1,-1},{1,-1},{1,1},{-1,1}},
//								{{-1,-1},{-1,1},{1,-1},{1,1}},
//								{{-1,-1},{-1,1},{1,1},{1,-1}},
//								{{-1,-1},{1,1},{-1,1},{1,-1}}};
//	
//	int nnt, b,tn,i,j,k,l,m,i2,stt,stp,minx,miny,maxx,maxy,cperm,nperm,MS,nlp,nj,ns;
//	int *ab,**TempOrigin,**TempDest,**TempNum;
//	long long TempSZ,ntd,nc,*TempCount;
//	int *RadiationLookupX,*RadiationLookupY;
//	double *RadiationLookupR;
//	double t,s,f,R2max;
//	float *cpb;
//	char outname[1024];
//	FILE *dat;
//
//	fprintf(stderr,"Initialising mobility vectors...\n");
//	if(!(State.InvCDFArray=(int *) malloc(P.NCP*(RAD_INVCDF+1)*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
//	ab=State.InvCDFArray;
//	nnt=0;
//	for(i=0;i<P.NCP;i++)
//		{
//		b=(int) (CellLookup[i]-Cells);
//		if(Cells[b].n>0)
//			{
//			Cells[b].InvCDFLookupRad=ab+nnt;
//			nnt+=RAD_INVCDF+1;
//			}
//		}
//	if(P.DoRadiationMobility)
//		{
//		if(!(dat=fopen(RadiationFile,"r"))) ERR_CRITICAL("Unable to open radiation mobility file\n");
//		fscanf(dat,"%lg",&R2max);
//		t=P.KernelRadiusMax/900;
//		if(t>R2max) t=R2max;
//		R2max=t*t;
//		MS=(int) ((t+1)*(t+1));
//		
//		if(!(RadiationLookupX=(int *) malloc(MS*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage 1\n");
//		if(!(RadiationLookupY=(int *) malloc(MS*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage 2\n");
//		if(!(RadiationLookupR=(double *) malloc(MS*sizeof(double)))) ERR_CRITICAL("Unable to allocate cell storage 3\n");
//		i=0;
//		while(!feof(dat))
//			{
//			fscanf(dat,"%i %i",&j,&k);
//			l=j*j+k*k;
//			if(l<=R2max)
//				{
//				RadiationLookupX[i]=j;
//				RadiationLookupY[i]=k;
//				RadiationLookupR[i]=sqrt((double) l)*900; //approx
//				i++;
//				}
//			}
//		nlp=i;
//		fclose(dat);
//		ntd=((long long) P.NumThreads);
//		TempSZ=RAD_TRIPS_PP*((long) P.N)/ntd/RAD_TEMP_MEM_SCALE;
//		if(!(TempOrigin=(int **) malloc(P.NumThreads*sizeof(int *)))) ERR_CRITICAL("Unable to allocate cell storage 4\n");
//		if(!(TempDest=(int **) malloc(P.NumThreads*sizeof(int *)))) ERR_CRITICAL("Unable to allocate cell storage 5\n");
//		if(!(TempNum=(int **) malloc(P.NumThreads*sizeof(int *)))) ERR_CRITICAL("Unable to allocate cell storage 6\n");
//		if(!(TempCount=(long long*) malloc(P.NumThreads*CACHE_LINE_SIZE*sizeof(long long)))) ERR_CRITICAL("Unable to allocate cell storage 7\n");
//		if(!(TempOrigin[0]=(int *) malloc(ntd*TempSZ*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage 8\n");
//		if(!(TempDest[0]=(int *) malloc(ntd*TempSZ*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage 9\n");
//		if(!(TempNum[0]=(int *) malloc(ntd*TempSZ*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage 10\n");
//		for(i=1;i<P.NumThreads;i++)
//			{
//			TempOrigin[i]=TempOrigin[i-1]+TempSZ;
//			TempDest[i]=TempDest[i-1]+TempSZ;
//			TempNum[i]=TempNum[i-1]+TempSZ;
//			}
//#pragma omp parallel for schedule(static,4096)
//		for(i=0;i<P.NC;i++) 
//			{
//			Cells[i].nn=Cells[i].tot_treat=0;
//			Cells[i].TotProbRad=0;
//			}
//#pragma omp parallel for private(tn,i,i2,b,j,k,l,m,stt,cperm,nperm,minx,miny,maxx,maxy,ns,nc,nj,s,t,f) schedule(static,1)
//		for(tn=0;tn<P.NumThreads;tn++)
//			{
//			nc=0;
//			for(i2=tn;i2<P.NCP;i2+=P.NumThreads)
//				{
//				if(i2%10000==0) fprintf(stderr,"Cell %i        \r",i2);
//				b=(int) (CellLookup[i2]-Cells);
//				minx=b/P.nch;
//				miny=b%P.nch;
//				stt=1;
//				i=0;
//				t=(double) Cells[b].n*P.MoveKernelP3+P.MoveKernelP4;  //John's radiation model ver 3 & 4
//				nj=RAD_TRIPS_PP*Cells[b].n;
//				while(stt)
//					{
//					j=RadiationLookupX[i];
//					k=RadiationLookupY[i];
//					if((j==0)||(k==0)) 
//						{nperm=2;cperm=((int) (2*ranf_mt(tn)));}
//					else
//						{nperm=4;cperm=(int) (24*ranf_mt(tn));}
//					for(l=0;(l<nperm)&&(stt);l++)
//						{
//						maxx=minx+ord[cperm][l][0]*j;
//						maxy=miny+ord[cperm][l][1]*k;
//						if((maxx>=0)&&(maxx<P.ncw)&&(maxy>=0)&&(maxy<P.nch))
//							{
//							m=maxy+P.nch*maxx;
//							if(Cells[m].n>0)
//								{
//								s=t+((double) Cells[m].n);
////								s=t+((double) Cells[m].n)*(P.MoveKernelP3+(1-P.MoveKernelP3)*(1.0-exp(-pow(RadiationLookupR[i]/P.MoveKernelScale,P.MoveKernelShape))));
//								f=1-t/s;
//								t=s;
////								s=f*((double) nj);
////								ns=(s<1e-4)?((ranf_mt(tn)<s)?1:0):((int) ignbin_mt((long) nj,f,tn));
//								ns=((int) ignbin_mt((long) nj,f,tn));
//								if(ns>0)
//									{
//									nj-=ns;
//									if(nj==0) stt=0;
//									TempOrigin[tn][nc]=b;
//									TempDest[tn][nc]=m;
//									TempNum[tn][nc]=ns;
//									nc++;
//									if(nc==TempSZ) ERR_CRITICAL("Storage for radiation model setup insufficient\n");
//									}
//								}
//							}
//						}
//					i++;
//					if(i==nlp) stt=0;
//					}
//				}
//			TempCount[tn*CACHE_LINE_SIZE]=nc;
//			}
//		nc=0;
//		for(tn=0;tn<P.NumThreads;tn++) nc+=TempCount[tn*CACHE_LINE_SIZE];
//		fprintf(stderr,"%li cell to cell connections formed.\n",nc);
//		nc+=(long long) P.NCP;
//		if(!(State.MobilityLinks=(int *) malloc(nc*sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage 11\n");
//		if(!(State.MobilityProbs=(float *) malloc(nc*sizeof(float)))) ERR_CRITICAL("Unable to allocate cell storage 12\n");
//		ab=State.MobilityLinks;
//		cpb=State.MobilityProbs;
//		State.NumMobilityLinks=nc;
//		for(tn=0;tn<P.NumThreads;tn++)
//			for(nc=0;nc<TempCount[tn*CACHE_LINE_SIZE];nc++)
//				// change this from TempDest to TempOrig, so that we can cdf of journeys from the cell instead of to the cell: ggilani 09/02/15
//				Cells[TempOrigin[tn][nc]].nn++;
//		i=0;
//		for(i2=0;i2<P.NCP;i2++)
//			{
//			b=(int) (CellLookup[i2]-Cells);
//			Cells[b].NeighboursIndexRad=ab+i;
//			Cells[b].CumProbRad=cpb+i;
//			Cells[b].NeighboursIndexRad[0]=b;
//			Cells[b].CumProbRad[0]=0;  //In our case, CumProb[0] stores the probability of within cell transmission// CumProb[0] stores prop of resident population staying away on any night
//			i+=(Cells[b].nn+1);
//			}
//		for(i=0;i<MAX_DIST;i++) PlaceDistDistrib[0][i]=0;
//		for(tn=0;tn<P.NumThreads;tn++)
//			{
//			for(nc=0;nc<TempCount[tn*CACHE_LINE_SIZE];nc++)
//				{
//				j=TempDest[tn][nc];
//				l=TempOrigin[tn][nc];
//				k=1+Cells[j].tot_treat;
//				//switched round j and l so that we are counting journies made from origin cell rather than to the cell
//				Cells[l].NeighboursIndexRad[k]=j;
//				if(k==1)
//					Cells[l].CumProbRad[k]=TempNum[tn][nc];
//				else
//					Cells[l].CumProbRad[k]=Cells[l].CumProbRad[k-1]+TempNum[tn][nc];
//				Cells[l].tot_treat++;
//				//Cells[l].CumProbRad[0]+=TempNum[tn][nc];
//				m=(int) (sqrt(dist2_cc(Cells+l,Cells+j))/OUTPUT_DIST_SCALE);
//				if(m<MAX_DIST) PlaceDistDistrib[0][m]+=TempNum[tn][nc];
//				}
//			}
//#pragma omp parallel for private(tn,i2,b,j,m) schedule(static,1)
//		for(tn=0;tn<P.NumThreads;tn++)
//			{
//			for(i2=tn;i2<P.NCP;i2+=P.NumThreads)
//				{
//				if(i2%10000==0) fprintf(stderr,"Cell %i        \r",i2);
//				b=(int) (CellLookup[i2]-Cells);
//				Cells[b].tot_treat=0;
//				Cells[b].CumProbRad[0]=1-P.PropNightsAway*Cells[b].CumProb[0]/((float) Cells[b].n)/RAD_TRIPS_PP;
//				Cells[b].TotProbRad=Cells[b].CumProbRad[Cells[b].nn];
//				for(j=1;j<=Cells[b].nn;j++)
//					Cells[b].CumProbRad[j]/=Cells[b].TotProbRad;
//				Cells[b].TotProbRad*=(P.PropNightsAway/RAD_TRIPS_PP);
//				m=(Cells[b].nn==0)?0:1;
//				for(j=0;j<=RAD_INVCDF;j++)
//					{
//					while((Cells[b].CumProbRad[m]*RAD_INVCDF<((float) j))&&(m<Cells[b].nn)) m++;
//					Cells[b].InvCDFLookupRad[j]=m;
//					}
//				}
//			}
//		free(TempOrigin[0]);
//		free(TempDest[0]);
//		free(TempNum[0]);
//		free(TempOrigin);
//		free(TempDest);
//		free(TempNum);
//		free(TempCount);
//		free(RadiationLookupX);
//		free(RadiationLookupY);
//		free(RadiationLookupR);
//		sprintf(outname,"%s.traveldist.xls",OutFile);
//		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
//		fprintf(dat,"dist\tfreq\n");
//		for(i=0;i<MAX_DIST;i++) fprintf(dat,"%i\t%i\n",i+1,PlaceDistDistrib[0][i]);
//		fclose(dat);
//		fprintf(stderr,"Mobility vectors initialised.\n");
//		}
//	
//}


#ifdef COUNTRY_WA
#define MEAN_CHILD_AGE_GAP 2.8 /* Mean number of years separating adjacent children */
#define MIN_ADULT_AGE 17
#define MAX_MF_PARTNER_AGE_GAP 7  /*  Max number of years older male partner is than female */
#define MAX_FM_PARTNER_AGE_GAP 2  /* Max number of years younger male partner is than female */
#define MIN_PARENT_AGE_GAP 16 /* Min number of years older a mother is than their child */
#define MAX_PARENT_AGE_GAP 33   /* Max number of years older a mother is than their child */
#define MAX_CHILD_AGE 29   //28   //30
#define ONE_CHILD_TWO_PERS_PROB 0.03
#define TWO_CHILD_THREE_PERS_PROB 0.0
#define TWO_PERS_HOUSE_PROB_OLD 0.5
#define ONE_PERS_HOUSE_PROB_OLD 0.5
#define ONE_PERS_HOUSE_PROB_YOUNG 0.0
#define TWO_PERS_HOUSE_PROB_YOUNG 0.0
#define ONE_CHILD_PROB_YOUNGEST_CHILD_UNDER_FIVE 0.9
#define TWO_CHILDREN_PROB_YOUNGEST_UNDER_FIVE 0.5
#define PROB_YOUNGEST_CHILD_UNDER_FIVE 0.45
#define ZERO_CHILD_THREE_PERS_PROB 0.35
#define ONE_CHILD_FOUR_PERS_PROB 0.25
#define YOUNG_AND_SINGLE 30
#define YOUNG_AND_SINGLE_SLOPE 0
#define NOCHILD_PERS_AGE 40 //30
#define OLD_PERS_AGE 55
#define THREE_CHILD_FIVE_PERS_PROB 0.5
#define FRAC_CHILDREN_BIG_HOUSEHOLDS 0.78
#define OLDER_GEN_GAP 8
#endif

#define PROP_OTHER_PARENT_AWAY 0.0

/* Complex household age distribution model
	- picks number of children (nc)
	- tries to space them reasonably
	- picks parental ages to be consistent with childrens' and each other
	- other adults in large housholds are assumed to be grandparents
	- for Thailand, 2 person households are 95% couples without children, 5% 1 parent families
*/
void AssignHouseholdAges(int n, int pers,int tn)
{
	int i,j,k,l,nc,ad;
	int a[MAX_HOUSEHOLD_SIZE+2];
	double f;

	ad=((P.DoAdunitDemog)&&(P.DoAdUnits))? Mcells[Hosts[pers].mcell].adunit:0;
	if(!P.DoHouseholds)
		{
		for(i=0;i<n;i++)
			a[i]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];
		}
	else
		{
		if(n==1)
			{
			if(ranf_mt(tn)<ONE_PERS_HOUSE_PROB_OLD)
				{
				do
					{a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];}
#ifdef COUNTRY_THAILAND
				while(a[0]<NOCHILD_PERS_AGE);
#else
				while((a[0]<NOCHILD_PERS_AGE)
					||(ranf_mt(tn)>(((double) a[0])-NOCHILD_PERS_AGE+1)/(OLD_PERS_AGE-NOCHILD_PERS_AGE+1)));
#endif
				}
			else if((ONE_PERS_HOUSE_PROB_YOUNG>0)&&(ranf_mt(tn)<ONE_PERS_HOUSE_PROB_YOUNG/(1-ONE_PERS_HOUSE_PROB_OLD)))
				{
				do
					{a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];}
				while((a[0]>YOUNG_AND_SINGLE)||(a[0]<MIN_ADULT_AGE)
					||(ranf_mt(tn)>1-YOUNG_AND_SINGLE_SLOPE*(((double) a[0])-MIN_ADULT_AGE)/(YOUNG_AND_SINGLE-MIN_ADULT_AGE)));
				}
			else
				while((a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))])<MIN_ADULT_AGE);
			}
		else if(n==2)
			{
			if(ranf_mt(tn)<TWO_PERS_HOUSE_PROB_OLD)
				{
				do
					{a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];}
#ifdef COUNTRY_THAILAND
				while(a[0]<NOCHILD_PERS_AGE);
#else
				while((a[0]<NOCHILD_PERS_AGE)
					||(ranf_mt(tn)>(((double) a[0])-NOCHILD_PERS_AGE+1)/(OLD_PERS_AGE-NOCHILD_PERS_AGE+1)));
#endif
				do
					{a[1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];}
#ifdef COUNTRY_THAILAND
				while((a[1]>=a[0]+MAX_MF_PARTNER_AGE_GAP)||(a[1]<a[0]-MAX_FM_PARTNER_AGE_GAP)||(a[1]<MIN_ADULT_AGE));
#else
				while((a[1]>a[0]+MAX_MF_PARTNER_AGE_GAP)||(a[1]<a[0]-MAX_FM_PARTNER_AGE_GAP)||(a[1]<NOCHILD_PERS_AGE)
					||(ranf_mt(tn)>(((double) a[1])-NOCHILD_PERS_AGE+1)/(OLD_PERS_AGE-NOCHILD_PERS_AGE+1)));
#endif
				}
			else if(ranf_mt(tn)<ONE_CHILD_TWO_PERS_PROB/(1-TWO_PERS_HOUSE_PROB_OLD))
				{
				while((a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))])>MAX_CHILD_AGE);
				do
					{a[1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];}
#ifdef COUNTRY_THAILAND
				while((a[1]>=a[0]+MAX_PARENT_AGE_GAP)||(a[1]<a[0]+MIN_PARENT_AGE_GAP));
#else
				while((a[1]>a[0]+MAX_PARENT_AGE_GAP)||(a[1]<a[0]+MIN_PARENT_AGE_GAP)||(a[1]<MIN_ADULT_AGE));
#endif
				}
			else if((TWO_PERS_HOUSE_PROB_YOUNG>0) &&(ranf_mt(tn)<TWO_PERS_HOUSE_PROB_YOUNG/(1-TWO_PERS_HOUSE_PROB_OLD-ONE_CHILD_TWO_PERS_PROB)))
				{
				do
					{a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];}
				while((a[0]<MIN_ADULT_AGE)||(a[0]>YOUNG_AND_SINGLE)
					||(ranf_mt(tn)>1-YOUNG_AND_SINGLE_SLOPE*(((double) a[0])-MIN_ADULT_AGE)/(YOUNG_AND_SINGLE-MIN_ADULT_AGE)));
				do
					{a[1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];}
#ifdef COUNTRY_THAILAND
				while((a[1]>=a[0]+MAX_MF_PARTNER_AGE_GAP)||(a[1]<a[0]-MAX_FM_PARTNER_AGE_GAP)||(a[1]<MIN_ADULT_AGE));
#else
				while((a[1]>a[0]+MAX_MF_PARTNER_AGE_GAP)||(a[1]<a[0]-MAX_FM_PARTNER_AGE_GAP)||(a[1]<MIN_ADULT_AGE));
#endif
				}
			else
				{
				do
					{a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];}
				while(a[0]<MIN_ADULT_AGE);
				do
					{a[1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];}
#ifdef COUNTRY_THAILAND
				while((a[1]>a[0]+MAX_MF_PARTNER_AGE_GAP)||(a[1]<a[0]-MAX_FM_PARTNER_AGE_GAP)||(a[1]<MIN_ADULT_AGE));
#else
				while((a[1]>a[0]+MAX_MF_PARTNER_AGE_GAP)||(a[1]<a[0]-MAX_FM_PARTNER_AGE_GAP)||(a[1]<MIN_ADULT_AGE));
#endif
				}

			}
		else
			{
			if(n==3)
				{
				if((ZERO_CHILD_THREE_PERS_PROB>0)||(TWO_CHILD_THREE_PERS_PROB>0))
					nc=(ranf_mt(tn)<ZERO_CHILD_THREE_PERS_PROB)?0:((ranf_mt(tn)<TWO_CHILD_THREE_PERS_PROB)?2:1);
				else
					nc=1;
				}
			else if(n==4)
				nc=(ranf_mt(tn)<ONE_CHILD_FOUR_PERS_PROB)?1:2;
			else if(n==5)
				nc=(ranf_mt(tn)<THREE_CHILD_FIVE_PERS_PROB)?3:2;
			else
#ifdef COUNTRY_INDONESIA
				do
					{
					nc=ignpoi_mt(FRAC_CHILDREN_BIG_HOUSEHOLDS*((double) n),tn);
					}
				while((nc>n-2)||(nc<2));
#else
				nc=n-2-(int) (3*ranf_mt(tn));
#endif
			if(nc==0)
				{
				do
					{
					a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];
					a[1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];
					}
#ifdef COUNTRY_THAILAND
				while((a[0]<MIN_ADULT_AGE)||(a[1]>=a[0]+MAX_PARENT_AGE_GAP)
					||(a[1]<a[0]+MIN_PARENT_AGE_GAP)||(a[1]>NOCHILD_PERS_AGE));
#else
				while((a[1]<MIN_ADULT_AGE)||(a[0]<MIN_ADULT_AGE));
#endif
				do
					{
					a[2]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];
					}
#ifdef COUNTRY_THAILAND
				while((a[2]>=a[1]+MAX_MF_PARTNER_AGE_GAP)||(a[2]<a[1]-MAX_FM_PARTNER_AGE_GAP)); 
#else
				while((a[2]>=a[1]+MAX_MF_PARTNER_AGE_GAP)||(a[2]<a[1]-MAX_FM_PARTNER_AGE_GAP)); 
#endif
				}
			else
				{
#ifdef COUNTRY_THAILAND
				a[0]=0;
				for(i=1;i<nc;i++)
					a[i]=a[i-1]+1+((int) ignpoi_mt(MEAN_CHILD_AGE_GAP-1,tn));
				j=a[nc-1]-(MAX_PARENT_AGE_GAP-MIN_PARENT_AGE_GAP);
				if(j>0)
					j+=MAX_PARENT_AGE_GAP;
				else
					j=MAX_PARENT_AGE_GAP;
				do
					{
					a[nc]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];
					a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];
					k=((nc>1)?a[nc-1]:0)+a[0];
					l=k-MAX_CHILD_AGE;
					if(ranf_mt(tn)<ONE_CHILD_PROB_YOUNGEST_CHILD_UNDER_FIVE) l=k-5;
					}
				while((l>0)||(a[nc]>a[0]+j)||(a[nc]<k+MIN_PARENT_AGE_GAP));
				for(i=1;i<nc;i++) a[i]+=a[0];
				if((n>nc+1)&&(ranf_mt(tn)>PROP_OTHER_PARENT_AWAY))
					{
					do
						{a[nc+1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];}
					while((a[nc+1]>=a[nc]+MAX_MF_PARTNER_AGE_GAP)
						||(a[nc+1]<a[nc]-MAX_FM_PARTNER_AGE_GAP));
					l=nc+2;
					}
				else
					l=nc+1;
#else
				do
					{
					a[0]=0;
					for(i=1;i<nc;i++)
						a[i]=a[i-1]+1+((int) ignpoi_mt(MEAN_CHILD_AGE_GAP-1,tn));
					a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))]-a[(int) (ranf_mt(tn)*((double) nc))];
					for(i=1;i<nc;i++) a[i]+=a[0];
					k=(((nc==1)&&(ranf_mt(tn)<ONE_CHILD_PROB_YOUNGEST_CHILD_UNDER_FIVE))||((nc==2)&&(ranf_mt(tn)<TWO_CHILDREN_PROB_YOUNGEST_UNDER_FIVE))
						||((nc>2)&&(ranf_mt(tn)<PROB_YOUNGEST_CHILD_UNDER_FIVE)))?5:MAX_CHILD_AGE;
					}
				while((a[0]<0)||(a[0]>k)||(a[nc-1]>MAX_CHILD_AGE));
				j=a[nc-1]-a[0]-(MAX_PARENT_AGE_GAP-MIN_PARENT_AGE_GAP);
				if(j>0)
					j+=MAX_PARENT_AGE_GAP;
				else
					j=MAX_PARENT_AGE_GAP;
				do
					{
					a[nc]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];
					k=a[nc-1];
					l=k-MAX_CHILD_AGE;
					}
				while((a[nc]>a[0]+j)||(a[nc]<k+MIN_PARENT_AGE_GAP)||(a[nc]<MIN_ADULT_AGE));
				if((n>nc+1)&&(ranf_mt(tn)>PROP_OTHER_PARENT_AWAY))
					{
					do
						{a[nc+1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))];}
					while((a[nc+1]>a[nc]+MAX_MF_PARTNER_AGE_GAP)||(a[nc+1]<a[nc]-MAX_FM_PARTNER_AGE_GAP)
						||(a[nc+1]>a[0]+j)||(a[nc+1]<k+MIN_PARENT_AGE_GAP)||(a[nc+1]<MIN_ADULT_AGE));
					l=nc+2;
					}
				else
					l=nc+1;
#endif
				if(n>nc+2)
					{
#ifdef COUNTRY_THAILAND
					j=a[nc]+MIN_PARENT_AGE_GAP;
#else
					j=((a[nc+1]>a[nc])?a[nc+1]:a[nc])+OLDER_GEN_GAP;
#endif
					if(j>=NUM_AGE_GROUPS*AGE_GROUP_WIDTH) j=NUM_AGE_GROUPS*AGE_GROUP_WIDTH-1;
					if(j<NOCHILD_PERS_AGE) j=NOCHILD_PERS_AGE;
					for(i=nc+2;i<n;i++)
						while((a[i]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(tn))])<j);
					}
				}
			}
		}
	for(i=0;i<n;i++)
		{
		Hosts[pers+i].age=(unsigned short int) a[i];
#pragma omp atomic
		AgeDist[HOST_AGE_GROUP(pers+i)]++;
		}
}

#define NO_EMP_WEIGHT 100

void AssignPeopleToPlaces(void)
{
	int i, i2, j, j2, k, k2, l, m, m2, tp, f, f2, f3, f4, ic, mx, my, a, cnt, tn, ca, nt, nn;
	int *PeopleArray;
	int *NearestPlaces[MAX_NUM_THREADS];
	double s, t, s2, *NearestPlacesProb[MAX_NUM_THREADS];
	cell *ct;
	int dbg;
	int npt;
	int nh_assigned, i3; //nh_assigned keeps track of the number of people in each household assigned to the same place, i3 is just another counting variable - ggilani 10/02/17
	int g, g1, g2, maxph; //added this variables to keep track of place variables when assigning members of the same household to the same place - ggilani 13/02/17

	npt = NUM_PLACE_TYPES;

	if (P.DoPlaces)
	{
		fprintf(stderr, "Assigning people to places....\n");
		for (i = 0; i < P.NC; i++)
		{
			Cells[i].infected = Cells[i].susceptible;
			if (!(Cells[i].susceptible = (int *)calloc(Cells[i].n, sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			Cells[i].cumTC = Cells[i].n;
		}

		//PropPlaces initialisation is only valid for non-overlapping places.

		for (i = 0; i < P.N; i++)
		{
			for (tp = 0; tp < npt; tp++) //Changed from 'for(tp=0;tp<P.PlaceTypeNum;tp++)' to try and assign -1 early and avoid problems when using less than the default number of placetypes later 
			{
				Hosts[i].PlaceLinks[tp] = -1;
				Hosts[i].PlacePresence[tp] = 0; //this is binary variable, so set to zero initially
			}
		}

		for (tp = 0; tp < P.PlaceTypeNum; tp++)
		{
			if (tp != HOTEL_PLACE_TYPE)
			{
				cnt = 0;
				for (a = 0; a < P.NCP; a++)
				{
					i = (int)(CellLookup[a] - Cells);
					Cells[i].n = 0;
					for (j = 0; j < Cells[i].cumTC; j++)
					{
						k = HOST_AGE_YEAR(Cells[i].members[j]);
						f = ((PropPlaces[k][tp] > 0) && (ranf() < PropPlaces[k][tp]));
						if (f)
						{
							for (k = 0; (k < tp) && (f); k++)
							{
								if ((Hosts[Cells[i].members[j]].PlaceLinks[k] >= 0) && (tp != P.HospPlaceTypeNum)) f = 0; //(ranf()<P.PlaceExclusivityMatrix[tp][k]); 
							}
							// Am assuming people can only belong to 1 place (and a hotel) at present
						}
						if (f)
						{
							Cells[i].susceptible[Cells[i].n] = Cells[i].members[j];
							Cells[i].n++;
							cnt++;
							Hosts[Cells[i].members[j]].PlacePresence[tp] = 1; //marking that this person does actually belong to this place type - will use this later to ensure we can assign some household members to similar places
						}
					}
					Cells[i].S = Cells[i].n;
					Cells[i].I = 0;
				}
				if (!(PeopleArray = (int*)calloc(cnt, sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
				j2 = 0;
				for (a = 0; a < P.NCP; a++)
				{
					i = (int)(CellLookup[a] - Cells);
					for (j = 0; j < Cells[i].n; j++)
					{
						PeopleArray[j2] = Cells[i].susceptible[j];
						j2++;
					}
				}
				for (i2 = 0; i2 < 2; i2++)
					for (k = 0; k < cnt; k++)
					{
						while ((l = (int)(((double)cnt) * ranf())) == k);
						j2 = PeopleArray[k];
						PeopleArray[k] = PeopleArray[l];
						PeopleArray[l] = j2;
					}
				m = 0;
				if (tp < P.nsp)
				{
					for (i = 0; i < P.Nplace[tp]; i++)
					{
						m += (int)(Places[tp][i].treat_end_time = (float)Places[tp][i].n);
						Places[tp][i].n = 0;
					}
				}
				else if (P.PlaceTypeSizePower[tp] == 0 && P.PlaceTypeSizeSD[tp] == 0)
				{
					for (i = 0; i < P.Nplace[tp]; i++)
					{
						Places[tp][i].n = 0;
						j = 1 + ((int)ignpoi(P.PlaceTypeMeanSize[tp] - 1));
						//if (j > USHRT_MAX - 1) j = USHRT_MAX - 1;
						m += (int)(Places[tp][i].nhcws = (float)j);
					}
				}
				//added this code to allow a place size to be specified according to a lognormal distribution - ggilani 09/02/17
				else if (P.PlaceTypeSizePower[tp] == 0 && P.PlaceTypeSizeSD[tp] > 0)
				{
					for (i = 0; i < P.Nplace[tp]; i++)
					{
						Places[tp][i].n = 0;
						j = (int)gen_lognormal(P.PlaceTypeMeanSize[tp], P.PlaceTypeSizeSD[tp]);
						//if (j > USHRT_MAX - 1) j = USHRT_MAX - 1;
						m += (int)(Places[tp][i].nhcws = (float)j); //changed this to nhcws (which hasn't been set yet) to overcome problems with ushrt and max size of places
					}
				}
				else
				{
					s = pow(P.PlaceTypeSizeOffset[tp] / (P.PlaceTypeSizeOffset[tp] + P.PlaceTypeSizeMax[tp] - 1), P.PlaceTypeSizePower[tp]);
					for (i = 0; i < P.Nplace[tp]; i++)
					{
						j = (int)floor(P.PlaceTypeSizeOffset[tp] * pow((1 - s) * ranf() + s, -1 / P.PlaceTypeSizePower[tp]) + 1 - P.PlaceTypeSizeOffset[tp]);
						if (j > USHRT_MAX - 1) j = USHRT_MAX - 1;
						m += (int)(Places[tp][i].treat_end_time = (float)j);
						Places[tp][i].n = 0;
					}
				}
				if (tp < P.nsp)
				{
					t = ((double)m) / ((double)P.Nplace[tp]);
					fprintf(stderr, "Adjusting place weights by cell (Capacity=%i Demand=%i  Av place size=%lg)\n", m, cnt, t);
					for (i = 0; i < P.Nplace[tp]; i++)
						if (Places[tp][i].treat_end_time > 0)
						{
							j = (int)(Places[tp][i].loc_x / P.cwidth);
							k = j * P.nch + ((int)(Places[tp][i].loc_y / P.cheight));
							Cells[k].I += (int)Places[tp][i].treat_end_time;
						}
					for (k = 0; k < P.NC; k++)
					{
						i = k % P.nch;
						j = k / P.nch;
						f2 = Cells[k].I; f3 = Cells[k].S;
						if ((i > 0) && (j > 0))
						{
							f2 += Cells[(j - 1) * P.nch + (i - 1)].I; f3 += Cells[(j - 1) * P.nch + (i - 1)].S;
						}
						if (i > 0)
						{
							f2 += Cells[j * P.nch + (i - 1)].I; f3 += Cells[j * P.nch + (i - 1)].S;
						}
						if ((i > 0) && (j < P.ncw - 1))
						{
							f2 += Cells[(j + 1) * P.nch + (i - 1)].I; f3 += Cells[(j + 1) * P.nch + (i - 1)].S;
						}
						if (j > 0)
						{
							f2 += Cells[(j - 1) * P.nch + i].I; f3 += Cells[(j - 1) * P.nch + i].S;
						}
						if (j < P.ncw - 1)
						{
							f2 += Cells[(j + 1) * P.nch + i].I; f3 += Cells[(j + 1) * P.nch + i].S;
						}
						if ((i < P.nch - 1) && (j > 0))
						{
							f2 += Cells[(j - 1) * P.nch + (i + 1)].I; f3 += Cells[(j - 1) * P.nch + (i + 1)].S;
						}
						if (i < P.nch - 1)
						{
							f2 += Cells[j * P.nch + (i + 1)].I; f3 += Cells[j * P.nch + (i + 1)].S;
						}
						if ((i < P.nch - 1) && (j < P.ncw - 1))
						{
							f2 += Cells[(j + 1) * P.nch + (i + 1)].I; f3 += Cells[(j + 1) * P.nch + (i + 1)].S;
						}
						Cells[k].L = f3; Cells[k].R = f2;
					}
					m = f2 = f3 = f4 = 0;
					for (k = 0; k < P.NC; k++)
						if ((Cells[k].S > 0) && (Cells[k].I == 0))
						{
							f2 += Cells[k].S; f3++;
							if (Cells[k].R == 0) f4 += Cells[k].S;
						}
					fprintf(stderr, "Demand in cells with no places=%i in %i cells\nDemand in cells with no places <=1 cell away=%i\n", f2, f3, f4);
					for (i = 0; i < P.Nplace[tp]; i++)
						if (Places[tp][i].treat_end_time > 0)
						{
							j = (int)(Places[tp][i].loc_x / P.cwidth);
							k = j * P.nch + ((int)(Places[tp][i].loc_y / P.cheight));
							if ((Cells[k].L > 0) && (Cells[k].R > 0))
							{
								s = ((double)Cells[k].L) / ((double)Cells[k].R);
								Places[tp][i].treat_end_time = ceil(Places[tp][i].treat_end_time * s);
							}
							m += ((int)Places[tp][i].treat_end_time);
						}
					for (i = 0; i < P.NC; i++) Cells[i].L = Cells[i].I = Cells[i].R = 0;
				}
				t = ((double)m) / ((double)P.Nplace[tp]);
				fprintf(stderr, "Adjusting place weights (Capacity=%i Demand=%i  Av place size=%lg)\n", m, cnt, t);
				for (i = m = 0; i < P.Nplace[tp]; i++)
				{
					s = ((double)Places[tp][i].treat_end_time) * 43 / 40 - 1;
					m += (int)(Places[tp][i].treat_end_time = 1 + ((float)ignpoi(s)));
				}
				if (tp < P.nsp)
					s = ((double)cnt) * 1.075;
				else
					s = ((double)cnt) * 1.125;
				j2 = ((int)s) - m;
				for (i = 0; i < j2; i++)
				{
					Places[tp][(int)(((double)P.Nplace[tp]) * ranf())].treat_end_time++;
				}
				j2 = -j2;
				for (i = 0; i < j2; i++)
				{
					while (Places[tp][j = (int)(((double)P.Nplace[tp]) * ranf())].treat_end_time < 2);
					Places[tp][j].treat_end_time--;
				}
				if (P.PlaceTypeNearestNeighb[tp] == 0)
				{
					for (i = 0; i < P.NC; i++) Cells[i].S = 0;
					for (j = 0; j < P.Nplace[tp]; j++)
					{
						i = P.nch * ((int)(Places[tp][j].loc_x / P.cwidth)) + ((int)(Places[tp][j].loc_y / P.cheight));
						Cells[i].S += (int)Places[tp][j].treat_end_time;
					}
					for (i = 0; i < P.NC; i++)
					{
						if (Cells[i].S > Cells[i].cumTC)
						{
							free(Cells[i].susceptible);
							if (!(Cells[i].susceptible = (int*)calloc(Cells[i].S, sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
						}
						Cells[i].S = 0;
					}
					for (j = 0; j < P.Nplace[tp]; j++)
					{
						i = P.nch * ((int)(Places[tp][j].loc_x / P.cwidth)) + ((int)(Places[tp][j].loc_y / P.cheight));
						k = (int)Places[tp][j].treat_end_time;
						for (j2 = 0; j2 < k; j2++)
						{
							Cells[i].susceptible[Cells[i].S] = j;
							Cells[i].S++;
						}
					}
				}
				for (i = 0; i < P.NumThreads; i++)
				{
					if (!(NearestPlaces[i] = (int*)calloc(P.PlaceTypeNearestNeighb[tp] + CACHE_LINE_SIZE, sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
					if (!(NearestPlacesProb[i] = (double*)calloc(P.PlaceTypeNearestNeighb[tp] + CACHE_LINE_SIZE, sizeof(double)))) ERR_CRITICAL("Unable to allocate cell storage\n");
				}
				P.KernelType = P.PlaceTypeKernelType[tp];
				P.KernelScale = P.PlaceTypeKernelScale[tp];
				P.KernelShape = P.PlaceTypeKernelShape[tp];
				P.KernelP3 = P.PlaceTypeKernelP3[tp];
				P.KernelP4 = P.PlaceTypeKernelP4[tp];
				InitKernel(1, 1.0);
				UpdateProbs(1);
				ca = 0;
				fprintf(stderr, "Allocating people to place type %i\n", tp);
				a = cnt;
				nt = P.NumThreads;
				nn = P.PlaceTypeNearestNeighb[tp];
				if (P.PlaceTypeNearestNeighb[tp] > 0)
				{
					//#pragma omp parallel for private(i,i2,j,j2,k,k2,l,m,m2,f,f2,ic,cnt,tn,s,t,mx,my,i3,nh_assigned) firstprivate(a,nt,nn) reduction(+:ca) schedule(static,1) //added i3, nh_assigned to private
					for (tn = 0; tn < P.NumThreads; tn++)
					{
						for (j = tn; j < a; j += nt)
						{
							if (j % 1000 == 0) fprintf(stderr, "(%i) %i      \r", tp, j);
							for (i2 = 0; i2 < nn; i2++)
								NearestPlacesProb[tn][i2] = 0;
							l = 1; k = m = m2 = f2 = 0;
							i = PeopleArray[j];
							ic = Hosts[i].mcell;
							mx = ic / P.nmch;
							my = ic % P.nmch;
							if (Hosts[i].PlaceLinks[tp] < 0) //added this so that if any hosts have already be assigned due to their household membership, they will not be reassigned
							{
								while (((k < nn) || (l < 4)) && (l < P.nmcw))
								{
									if ((mx >= 0) && (my >= 0) && (mx < P.nmcw) && (my < P.nmch))
									{
										ic = mx * P.nmch + my;
										if (Mcells[ic].country == Mcells[Hosts[i].mcell].country)
										{
											for (cnt = 0; cnt < Mcells[ic].np[tp]; cnt++)
											{
												if (Mcells[ic].places[tp][cnt] >= P.Nplace[tp]) fprintf(stderr, "#%i %i %i  ", tp, ic, cnt);
												t = dist2_raw(Households[Hosts[i].hh].loc_x, Households[Hosts[i].hh].loc_y,
													Places[tp][Mcells[ic].places[tp][cnt]].loc_x, Places[tp][Mcells[ic].places[tp][cnt]].loc_y);
												s = numKernel(t);
												if (tp < P.nsp)
												{
													t = ((double)Places[tp][Mcells[ic].places[tp][cnt]].treat_end_time);
													if (HOST_AGE_YEAR(i) < P.PlaceTypeMaxAgeRead[tp])
													{
														if ((t > 0) && (Places[tp][Mcells[ic].places[tp][cnt]].AvailByAge[HOST_AGE_YEAR(i)] > 0))
															s *= t;
														else
															s = 0;
													}
													else if (t > 0)
														s *= t;
												}
												k2 = 0; j2 = 0; t = 1e10;
												if (s > 0)
												{
													if (k < nn)
													{
														NearestPlaces[tn][k] = Mcells[ic].places[tp][cnt];
														NearestPlacesProb[tn][k] = s;
														k++;
													}
													else
													{
														for (i2 = 0; i2 < nn; i2++)
														{
															if (NearestPlacesProb[tn][i2] < t)
															{
																t = NearestPlacesProb[tn][i2]; j2 = i2;
															}
														}
														if (s > t)
														{
															NearestPlacesProb[tn][j2] = s;
															NearestPlaces[tn][j2] = Mcells[ic].places[tp][cnt];
														}
													}
												}
											}
										}
									}
									if (m2 == 0)
										mx = mx + 1;
									else if (m2 == 1)
										my = my - 1;
									else if (m2 == 2)
										mx = mx - 1;
									else if (m2 == 3)
										my = my + 1;
									f2 = (f2 + 1) % l;
									if (f2 == 0)
									{
										m2 = (m2 + 1) % 4;
										m = (m + 1) % 2;
										if (m == 0) l++;
									}
								}

								s = 0;
								if (k > nn) fprintf(stderr, "*** k>P.PlaceTypeNearestNeighb[tp] ***\n");
								if (k == 0)
								{
									fprintf(stderr, "# %i %i     \r", i, j);
									Hosts[i].PlaceLinks[tp] = -1;
								}
								else
								{
									for (i2 = 1; i2 < k; i2++)
										NearestPlacesProb[tn][i2] += NearestPlacesProb[tn][i2 - 1];
									s = NearestPlacesProb[tn][k - 1];
									t = ranf_mt(tn);
									f = 0;
									for (i2 = 0; (i2 < k) && (!f); i2++)
									{
										if ((f = (t < NearestPlacesProb[tn][i2] / s)))
										{
											Hosts[i].PlaceLinks[tp] = NearestPlaces[tn][i2];
											ca++;
											//add something about assigning family members here - ggilani 10/02/17
											if (!(P.HospPlaceTypeNum) && (P.PlaceHouseholdOverlap == 1))
											{
												nh_assigned = 1;
												i3 = 0;

												maxph = (int)ceil(Households[Hosts[i].hh].nh / P.PlaceHouseholdDivisor); //changed from 2.0 to 3.0

												while ((nh_assigned < maxph) && (i3 < Households[Hosts[i].hh].nh)) //so that households are roughly divided between places
												{
													if ((Hosts[Households[Hosts[i].hh].FirstPerson + i3].PlaceLinks[tp] < 0) && (Hosts[Households[Hosts[i].hh].FirstPerson + i3].PlacePresence[tp] == 1))
													{
														//then assign this person to the same place
														Hosts[Households[Hosts[i].hh].FirstPerson + i3].PlaceLinks[tp] = NearestPlaces[tn][i2];
														nh_assigned++;
														ca++; //think this should be here - need to check! ggilani 10/02/2017
													}
													i3++;
												}
											}
											if (tp < P.nsp)
											{
												//#pragma omp critical (places_treat_time)
												Places[tp][Hosts[i].PlaceLinks[tp]].treat_end_time--;
											}
										}
										if (!f) Hosts[i].PlaceLinks[tp] = -1;
										if (NearestPlaces[tn][i2] >= P.Nplace[tp]) fprintf(stderr, "@%i %i %i  ", tp, i, j);
									}
								}
							}
						}

					}
				}
				else
				{
					k2 = cnt - ca;
					m2 = cnt;
					s2 = 0;
					a = k2 / 1000;
					f = k2;
					for (ic = 0; ic <= 30; ic++)
					{
						UpdateProbs(1);
						m2 = f - 1;
						if (ic < 9)
							f = 100 * (9 - ic) * a;
						else if (ic < 18)
							f = 10 * (18 - ic) * a;
						else if (ic < 27)
							f = (27 - ic) * a;
						else
						{
							m2 = k2 - 1; f = 0;
						}
						//#pragma omp parallel for private(i,i2,j,k,l,m,f2,f3,t,ct,s,g,g1,g2,i3,nh_assigned) reduction(+:ca) /* schedule(dynamic,500)*/ //add s to private variables, added g,g1,g2,i3 and nh_assigned to private variables
						for (i2 = m2; i2 >= f; i2--)
						{
							if (i2 % 1000 == 0)
								fprintf(stderr, "(%i) %i            \r", tp, i2);
							k = PeopleArray[i2];
							i = Hosts[k].pcell;
							f2 = 1;
							f3 = (HOST_AGE_YEAR(k) >= P.PlaceTypeMaxAgeRead[tp]);
							if (Hosts[k].PlaceLinks[tp] < 0)
								while ((f2 > 0) && (f2 < 1000))
								{
									do
									{
										s = ranf();
										l = Cells[i].InvCDF[(int)floor(s * 1024)];
										while (Cells[i].cum_trans[l] < s) l++;
										ct = CellLookup[l];
										m = (int)(ranf() * ((double)ct->S));
										j = -1;
										//#pragma omp critical
										{
											if (ct->susceptible[m] >= 0)
												if ((f3) || (Places[tp][ct->susceptible[m]].AvailByAge[HOST_AGE_YEAR(k)] > 0))
												{
													j = ct->susceptible[m];
													ct->susceptible[m] = -1;
												}
										}
									} while (j < 0);
									if (j >= P.Nplace[tp])
									{
										fprintf(stderr, "*%i %i: %i %i\n", k, tp, j, P.Nplace[tp]);
									}
									t = dist2_raw(Households[Hosts[k].hh].loc_x, Households[Hosts[k].hh].loc_y, Places[tp][j].loc_x, Places[tp][j].loc_y);
									s = ((double)ct->S) / ((double)ct->S0) * numKernel(t) / Cells[i].max_trans[l];
									if ((P.DoAdUnits) && (P.InhibitInterAdunitPlaceAssignment[tp] > 0))
									{
										if (Mcells[Hosts[k].mcell].adunit != Mcells[Places[tp][j].mcell].adunit) s *= (1 - P.InhibitInterAdunitPlaceAssignment[tp]);
									}
									if (ranf() < s)
									{
										//#pragma omp critical
										l = (--ct->S);
										if (m < l) ct->susceptible[m] = ct->susceptible[l];
										//#pragma omp critical (places_treat_time)
										Places[tp][j].treat_end_time--;
										ca++;
										Hosts[k].PlaceLinks[tp] = j;
										//add something about assigning family members here - ggilani 10/02/17
										if (!(P.HospPlaceTypeNum) && (P.PlaceHouseholdOverlap == 1))
										{
											nh_assigned = 1;
											i3 = 0;

											while ((nh_assigned <= (int)ceil(Households[Hosts[k].hh].nh / P.PlaceHouseholdDivisor)) && (i3 < Households[Hosts[k].hh].nh)) //so that households are roughly divided between places - changed from 2.0 to 3.0

											{
												//also need to check for space in this scenario???
												if ((Hosts[Households[Hosts[k].hh].FirstPerson + i3].PlaceLinks[tp] < 0) && (Hosts[Households[Hosts[k].hh].FirstPerson + i3].PlacePresence[tp] == 1) && (Places[tp][j].treat_end_time > 0))
												{
													//need to find a value of m for which place is equal to j
													g = -1; //this is flag variable
													g1 = m; //g1, g2 are counting variables that allow us to move away from m and find the next closest 
													g2 = m;
													while (g < 0)
													{
														g1++;
														g2--;
														if (g1 < (ct->S))
														{
															if (ct->susceptible[g1] == j)
															{
																m = g1;
																g = 0;
															}
														}
														else if (g2 >= 0)
														{
															if (ct->susceptible[g2] == j)
															{
																m = g2;
																g = 0;
															}
														}
													}
													//then assign this person to the same place
//#pragma omp critical
													l = (--ct->S);
													if (m < l) ct->susceptible[m] = ct->susceptible[l];
													Hosts[Households[Hosts[k].hh].FirstPerson + i3].PlaceLinks[tp] = j;
													nh_assigned++;
													//#pragma omp critical (places_treat_time)
													Places[tp][j].treat_end_time--;
													ca++;
												}
												i3++;
											}
										}
										f2 = 0;
									}
									else
									{
										ct->susceptible[m] = j;
										f2++;
									}
								}
						}
					}
				}
				fprintf(stderr, "%i hosts assigned to placetype %i\n", ca, tp);
				free(PeopleArray);
				for (i = 0; i < P.Nplace[tp]; i++)
				{
					Places[tp][i].treat_end_time = 0;
					Places[tp][i].n = 0;
				}
				for (i = 0; i < P.NumThreads; i++)
				{
					free(NearestPlacesProb[i]);
					free(NearestPlaces[i]);
				}
			}
		}

			for (i = 0; i < P.NC; i++)
			{
				Cells[i].n = Cells[i].cumTC;
				Cells[i].cumTC = 0;
				Cells[i].S = Cells[i].I = Cells[i].L = Cells[i].R = 0;
				free(Cells[i].susceptible);
				Cells[i].susceptible = Cells[i].infected;
			}
			P.KernelScale = P.MoveKernelScale;
			P.KernelShape = P.MoveKernelShape;
			P.KernelType = P.MoveKernelType;
			P.KernelP3 = P.MoveKernelP3;
			P.KernelP4 = P.MoveKernelP4;
		

	}
}


void StratifyPlaces(void)
{
	int i,j,k,l,m,n,tn;
	double t;
	int nhcws, nflws;

	if(P.DoPlaces)
		{
		fprintf(stderr,"Initialising groups in places\n");
#pragma omp parallel for private(i,j) schedule(static,500)
		for(i=0;i<P.N;i++)
			for (j = 0; j < NUM_PLACE_TYPES; j++)
			{
				Hosts[i].PlaceGroupLinks[j] = 0;
			}
		for(j=0;j<P.PlaceTypeNum;j++)
			for (i = 0; i < P.Nplace[j]; i++)
			{
				Places[j][i].n = 0;
				Places[j][i].nhcws = 0;
				Places[j][i].nflws = 0;
				Places[j][i].n_current = 0;
			}

#pragma omp parallel for private(i,j,k,l,m,n,t,tn,nhcws,nflws) schedule(static,1)
		for(tn=0;tn<P.NumThreads;tn++)
		for(j=tn;j<P.PlaceTypeNum;j+=P.NumThreads)
			{
			if(j==HOTEL_PLACE_TYPE)
				{
				l=2*((int) P.PlaceTypeMeanSize[j]);
				for(i=0;i<P.Nplace[j];i++)
					{
					if(!(Places[j][i].members=(int *) calloc(l,sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
					Places[j][i].n=0;
					}
				}
			else
				{
				for(i=0;i<P.N;i++)
					{
					if(Hosts[i].PlaceLinks[j]>=0)
						Places[j][Hosts[i].PlaceLinks[j]].n++;
					}
				for(i=0;i<P.Nplace[j];i++)
					{
					if(Places[j][i].n>0)
						{
						if(!(Places[j][i].members=(int *) calloc(Places[j][i].n,sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
						if (!(Places[j][i].flwmembers = (int*)calloc(Places[j][i].n, sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
						}
					Places[j][i].n=0;
					}
				for(i=0;i<P.N;i++)
					{
					k=Hosts[i].PlaceLinks[j];
					if(k>=0)
						{
						Places[j][k].members[Places[j][k].n]=i;
						Places[j][k].n++;
						}
					}
				if(!P.PlaceTypeDoHousehold[j])
					for(k=0;k<P.Nplace[j];k++)
						for(i=0;i<Places[j][k].n;i++)
							if(!(Hosts[Places[j][k].members[i]].nc_plus_hh_disabled & HH_DISABLED))
								Hosts[Places[j][k].members[i]].nc_plus_hh_disabled=HH_DISABLED;
				for(i=0;i<P.Nplace[j];i++)
					if(Places[j][i].n>0)
						{
						t=((double) Places[j][i].n)/P.PlaceTypeGroupSizeParam1[j]-1.0;
						if((t<0)||(j==P.HospPlaceTypeNum))
							Places[j][i].ng=1;
						else
							Places[j][i].ng=1+(int) ignpoi_mt(t,tn);
						if(!(Places[j][i].group_start=(int *) calloc(Places[j][i].ng,sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
						if(!(Places[j][i].group_size=(int *) calloc(Places[j][i].ng,sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
						m=Places[j][i].n-Places[j][i].ng;
						for(k=l=0;k<Places[j][i].ng;k++)
							{
							t=1/((double) (Places[j][i].ng-k));
							Places[j][i].group_start[k]=l;
							Places[j][i].group_size[k]=1+ignbin_mt((long) m,t,tn);
							m-=(Places[j][i].group_size[k]-1);
							l+=Places[j][i].group_size[k];
							}
						for(k=0;k<Places[j][i].n;k++)
							{
							l=(int) (((double) Places[j][i].n)*ranf_mt(tn));
							n=Places[j][i].members[l];
							Places[j][i].members[l]=Places[j][i].members[k];
							Places[j][i].members[k]=n;
							}
						for(k=l=0;k<Places[j][i].ng;k++)
							for(m=0;m<Places[j][i].group_size[k];m++)
								{
								Hosts[Places[j][i].members[l]].PlaceGroupLinks[j]=k;
								l++;
								}
						}
				}
			}

		//vaccinate hcws and assign flws
	//assign frontline workers and adjust susceptibilities
		if (P.IncludeHospitalPlaceType)
		{
			for (i = 0; i < P.Nplace[P.HospPlaceTypeNum]; i++)
			{
				//first, calculate expected number of hcws and flws for the place, based on place size and hcws/flws ratios

				nhcws = (int)min((((double)Places[P.HospPlaceTypeNum][i].n / (double)1000) * P.HCWPerThousand),1); //added min to ensure we don't have zero hcws/flws
				nflws = (int)min((((double)Places[P.HospPlaceTypeNum][i].n / (double)1000) * P.FLWPerThousand),1);

				//assign healthcare workers and rearrange
				for (j = 0; j < nhcws;)
				{
					//select random place member
					l = (int)(((double)Places[P.HospPlaceTypeNum][i].n) * ranf());
					k = Places[P.HospPlaceTypeNum][i].members[l];
					//check to see if person is in correct age range to be a hcw and is not a hcw already
					if ((Hosts[k].age > P.MinAgeHCWFLW) && (Hosts[k].age < P.MaxAgeHCWFLW) && (Hosts[k].hcw == 0))
					{
						Hosts[k].hcw = Hosts[k].keyworker = 1;
						Places[P.HospPlaceTypeNum][i].members[l] = Places[P.HospPlaceTypeNum][i].members[j];
						Places[P.HospPlaceTypeNum][i].members[j] = k;
						Places[P.HospPlaceTypeNum][i].nhcws++;
						j++;
					}
				}
				Places[P.HospPlaceTypeNum][i].n_current = Places[P.HospPlaceTypeNum][i].nhcws;

				//set max capacity
				Places[P.HospPlaceTypeNum][i].maxcapacity = Places[P.HospPlaceTypeNum][i].nhcws + P.HospCaseCapacity;

				if (P.IncludeFLWs)
				{
					for (j = 0; j < nflws;)
					{
						//select random place member
						l = (int)(((double)Places[P.HospPlaceTypeNum][i].n) * ranf());
						k = Places[P.HospPlaceTypeNum][i].members[l];
						//check to see if person is in correct age range to be a hcw and is not a hcw already
						if ((Hosts[k].age > P.MinAgeHCWFLW) && (Hosts[k].age < P.MaxAgeHCWFLW) && (Hosts[k].hcw == 0) && (Hosts[k].flw == 0))
						{
							Hosts[k].flw = Hosts[k].keyworker = 1;
							Places[P.HospPlaceTypeNum][i].flwmembers[j] = k;
							Hosts[k].susc *= P.RelSuscFLW;
							Places[P.HospPlaceTypeNum][i].nflws++;
							j++;
						}
					}
				}
			}
		}
		for(i=0;i<P.N;i++)
			if(Hosts[i].nc_plus_hh_disabled & HH_DISABLED) Households[Hosts[i].hh].nhr--;
#pragma omp parallel for private(i,j,k,l) schedule(static,1)
		for(i=0;i<P.NumThreads;i++)
			{
			for(k=0;k<P.PlaceTypeNum;k++)
				{
				if(P.DoPlaceGroupTreat)
					{
					l=0;
					for(j=0;j<P.Nplace[k];j++)
						l+=(int) Places[k][j].ng;
					m=2*m/P.NumThreads;
					if(!(StateT[i].p_queue[k]=(int *) calloc(l,sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
					if(!(StateT[i].pg_queue[k]=(int *) calloc(l,sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
					}
				else
					{
					if(!(StateT[i].p_queue[k]=(int *) calloc(P.Nplace[k],sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
					if(!(StateT[i].pg_queue[k]=(int *) calloc(P.Nplace[k],sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
					}
				}
			}
		fprintf(stderr,"Groups initialised\n");
/*		s2=t2=0;
		for(j=0;j<P.PlaceTypeNum;j++)
			{
			t=s=0;
			for(i=0;i<P.Nplace[j];i++)
				if(Places[j][i].ng>0)
					{
					for(k=0;k<Places[j][i].ng;k++)
						t+=(double) Places[j][i].group_size[k];
					s+=(double) Places[j][i].ng;
					}
			s2+=s;
			t2+=t;
			fprintf(stderr,"Mean group size for place type %i = %lg\n",j,t/s);
			}
		t=0;
		for(i=0;i<P.N;i++)
			for(j=0;j<P.PlaceTypeNum;j++)
				if(Hosts[i].PlaceLinks[j]>=0)
					t+=(double) Places[j][Hosts[i].PlaceLinks[j]].group_size[Hosts[i].PlaceGroupLinks[j]];
		fprintf(stderr,"Overall mean group size = %lg (%lg)\n",t/((double) P.N),t2/s2);
*/		}
}


void LoadPeopleToPlaces(char *NetworkFile)
{
	int i,j,k,l,m,n,npt,i2;
	long s1,s2;
	FILE *dat;

	if(!(dat=fopen(NetworkFile,"rb"))) ERR_CRITICAL("Unable to open network file\n");
	npt=NUM_PLACE_TYPES_NOAIR;
	fread_big(&i,sizeof(int),1,dat);
	fread_big(&s1,sizeof(long),1,dat);
	fread_big(&s2,sizeof(long),1,dat);
	if(i!=npt) ERR_CRITICAL("Number of place types does not match saved value\n");
	if((s1!=P.seed1)||(s2!=P.seed2)) ERR_CRITICAL("Random number seeds do not match saved values\n");
	k=(P.N+999999)/1000000;
	for(i=0;i<P.N;i++)
		for(j=0;j<P.PlaceTypeNum;j++)
			Hosts[i].PlaceLinks[j]=-1;
	for(i=i2=0;i<k;i++)
		{
		l=(i<k-1)?1000000:(P.N-1000000*(k-1));
		fread_big(&netbuf,sizeof(int),npt*l,dat);
		for(j=0;j<l;j++)
			{
			n=j*npt;
			for(m=0;m<npt;m++)
				{
				Hosts[i2].PlaceLinks[m]=netbuf[n+m];
				if(Hosts[i2].PlaceLinks[m]>=P.Nplace[m])
					fprintf(stderr,"*%i %i: %i %i",i2,m,Hosts[i2].PlaceLinks[m],P.Nplace[m]);
				}
			i2++;
			}	
		fprintf(stderr,"%i loaded            \r",i*1000000+l);
		}

/*	for(i=0;i<P.N;i++)
		{
		if((i+1)%100000==0) fprintf(stderr,"%i loaded            \r",i+1);
		fread_big(&(Hosts[i].PlaceLinks[0]),sizeof(int),P.PlaceTypeNum,dat);
		}
*/	fprintf(stderr,"\n");
	fclose(dat);
	P.KernelScale=P.MoveKernelScale;
	P.KernelShape=P.MoveKernelShape;
	P.KernelP3=P.MoveKernelP3;
	P.KernelP4=P.MoveKernelP4;
}

void SavePeopleToPlaces(char *NetworkFile)
{
	int i,j,npt;
	FILE *dat;

	npt=NUM_PLACE_TYPES_NOAIR;
	if(!(dat=fopen(NetworkFile,"wb"))) ERR_CRITICAL("Unable to open network file\n");
	if(P.PlaceTypeNum>0)
		{
		fwrite_big(&npt,sizeof(int),1,dat);
		fwrite_big(&P.seed1,sizeof(long),1,dat);
		fwrite_big(&P.seed2,sizeof(long),1,dat);
		for(i=0;i<P.N;i++)
			{
			if((i+1)%100000==0) fprintf(stderr,"%i saved            \r",i+1);
/*			fwrite_big(&(Hosts[i].spatial_norm),sizeof(float),1,dat);
*/			fwrite_big(&(Hosts[i].PlaceLinks[0]),sizeof(int),npt,dat);
			for(j=0;j<npt;j++)
				if(Hosts[i].PlaceLinks[j]>=P.Nplace[j])
					fprintf(stderr,"*%i %i: %i %i",i,j,Hosts[i].PlaceLinks[j],P.Nplace[j]);
			}
		}
	fprintf(stderr,"\n");
	fflush(dat);
	fclose(dat);
}

