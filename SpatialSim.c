/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "SpatialSim.h"
#include "binio.h"
#include "ranf.h"

void ReadParams(char *,char *);
void ReadInterventions(char *);
int GetXMLNode(FILE *,char *,char *,char *,int);
void ReadAirTravel(char *);
void SetupModel(char *,char *,char *,char *);
void SetupPopulation(char *,char *,char *);
void SetupAirports(void);
void SetupRoads(void); //added new function to take care of roads: ggilani - 12/02/15
void InitKernel(int,double);
void AssignHouseholdAges(int, int,int);
void AssignPeopleToPlaces(void);
void StratifyPlaces(void);
void LoadPeopleToPlaces(char *);
void SavePeopleToPlaces(char *);
void InitModel(int); //adding run number as a parameter for event log: ggilani - 15/10/2014
void SeedInfection(double, int *,int, int); //adding run number as a parameter for event log: ggilani - 15/10/2014
int CalcSeedResist(void);
void RunModel(int); //adding run number as a parameter for event log: ggilani - 15/10/2014
void TravelReturnSweep(double);
void TravelDepartSweep(double);
double CalcHouseInf(int , unsigned short int );
double CalcPlaceInf(int , int , unsigned short int );
double CalcSpatialInf(int, unsigned short int );
double CalcPersonInf(int , unsigned short int );
double CalcHouseSusc(int, unsigned short int,int,int);
void EquilibPersonDemog(int,int);
void UpdatePersonDemog(int ,int );
void UpdatePersonDemogSIR(int ,int );
void UpdateSuscPersonDemogSIR(int ,int );
double CalcPlaceSusc(int ,int, unsigned short int,int,int);
double CalcSpatialSusc(int, unsigned short int,int,int);
double CalcPersonSusc(int, unsigned short int,int,int);
void InfectSweep(double, int); //added int as argument to InfectSweep to record run number: ggilani - 15/10/14
void IncubRecoverySweep(double,int); //added int as argument to record run number: ggilani - 15/10/14
void DemogSweep(double);
void SIASweep(unsigned short int);
int TreatSweep(double);
//void HospitalSweep(double); //added hospital sweep function: ggilani - 10/11/14
void HospitalSweepAdunits(double); //added hospitalisation by adunit sweep function: ggilani - 24/11/14
void ContactTracingSweep(double); // added function to update contact tracing number
void VaccSweep(double); //added function to process ring vaccination queue: ggilani - 21/08/19
void UpdateHospitals(double); //added function to update hospital parameters on each time step: ggilani - 11/03/2017
void UpdateContactTracing(double); //added function to update contact tracing capacity
void UpdateSDB(double); //added function to update safe burial capacity
void UpdateVaccination(double,int); //added function to update vaccination parameters at each time step: ggilani - 29/05/2019
void UpdateCaseDetection(double); //added function to update vaccination parameters at each time step: ggilani - 29/05/2019
void SaveAgeDistrib(void);
void SaveAgeDistrib2(void);
void SaveDistribs(void);
void SaveOriginDestMatrix(void); //added function to save origin destination matrix so it can be done separately to the main results: ggilani - 13/02/15
void SaveResults(void);
void SaveSummaryResults(void);
void SaveRandomSeeds(void); //added this function to save random seeds for each run: ggilani - 09/03/17
void SaveEvents(void); //added this function to save infection events from all realisations: ggilani - 15/10/14
void LoadSnapshot(void);
void SaveSnapshot(void);
void UpdateProbs(int);
void RecordInfTypes(void);
void RecordSample(double,int);
//adding function to record an event: ggilani - 10/10/2014
void RecordEvent(double, int, int, int, int); //added int as argument to InfectSweep to record run number: ggilani - 15/10/14
void DoImmune(int);
void DoReborn(int);
void DoInfect(int,double,int,int); //added int as argument to InfectSweep to record run number: ggilani - 15/10/14
void DoIncub(int,unsigned short int,int,int); //added int as argument to record run number: ggilani - 23/10/14
void DoDetectedCase(int,double,unsigned short int,int);
void DoCase(int,double,unsigned short int,int);
void DoFalseCase(int ,double,unsigned short int,int );
void DoRecover(int,int,int); //added int as argument to record run number: ggilani - 23/10/14
void DoDeath(int,int,int); //added int as argument to record run number: ggilani - 23/10/14
void DoPlaceClose(int ,int ,unsigned short int,int,int);
void DoTreatCase(int,unsigned short int,int);
void DoProph(int,unsigned short int,int);
void DoPrivateTreatCase(int,unsigned short int,int);
void DoPrivateProph(int,unsigned short int,int);
void DoProphNoDelay(int,unsigned short int,int,int);
int DoVacc(int,int,int);
void DoVaccNoDelay(int,int);
void CalcOriginDestMatrix_adunit(void); //added function to calculate origin destination matrix: ggilani 28/01/15
int BedsAvailablePerAdUnit(double,int); //added function to check current number of beds available in an admin unit
void DetermineCellsWithCapitalCities(void); //added function to quickly check capital cities
double ExpKernel(double);
double PowerKernel(double);
double PowerKernelB(double);
double PowerKernelUS(double);
double PowerExpKernel(double);
double GaussianKernel(double);
double StepKernel(double);
double numKernel(double);
double dist2UTM(double ,double ,double ,double );
double dist2(person *, person *);
double dist2_cc(cell *, cell *);
double dist2_cc_min(cell *, cell *);
double dist2_mm(microcell *, microcell *);
double dist2_mm_min(microcell *, microcell *); //added min distance from microcell to microcell for origin-destination matrix: ggilani 28/01/15
double dist2_raw(double,double,double,double);
int GetInputParameter(FILE *,FILE *,char *,char *,void *, int,int,int);
int GetInputParameter2(FILE *,FILE *,char *,char *,void *, int,int,int);
int GetInputParameter3(FILE*, const char*, const char*, void*, int, int, int);
//int GetInputParameter3(FILE *,char *,char *,void *, int,int,int);
void CaptureBitmap(int,int);
void CaptureMeanBitmap(int);
void OutputBitmap(double,int);
void InitBMHead();
void HSB2RGB(double, double, double, int *, int *, int *);
void HandleBreak(int);
double PrevalenceDepTransmission(int);

#ifdef FRESSCA
void FYShuffle_mt(int*, int, int);
void SetupDistrNet(void);
void DistrVaccSweep(double );
int DoDistribVacc(int,unsigned short int);
void UpdateVaccStatus(int ,int );
#endif

param P;
person *Hosts;
household *Households;
popvar State,StateT[MAX_NUM_THREADS];
cell *Cells, **CellLookup;
int *RevCellLookup;
microcell *Mcells,**McellLookup;
place **Places;
adminunit AdUnits[MAX_ADUNITS];
results *TimeSeries,*TSMean,*TSVar,*TSMeanNE,*TSVarNE,*TSMeanE,*TSVarE;
results *TimeSeries_G,*TSMean_G,*TSVar_G,*TSMeanNE_G,*TSVarNE_G,*TSMeanE_G,*TSVarE_G;
results *TimeSeries_L,*TSMean_L,*TSVar_L,*TSMeanNE_L,*TSVarNE_L,*TSMeanE_L,*TSVarE_L;

//added declaration of pointer to events log: ggilani - 10/10/2014
events *InfEventLog;
int *nEvents;

airport *Airports;
double (*Kernel)(double );
double *nKernel,*nKernelHR,**PopDensity,*mcell_dens;
int *mcell_adunits,*mcell_num,*mcell_country;
double inftype[INFECT_TYPE_MASK],inftype_av[INFECT_TYPE_MASK],infcountry[MAX_COUNTRIES],infcountry_av[MAX_COUNTRIES],infcountry_num[MAX_COUNTRIES];
double indivR0[MAX_SEC_REC][MAX_GEN_REC],indivR0_av[MAX_SEC_REC][MAX_GEN_REC];
double inf_household[MAX_HOUSEHOLD_SIZE+1][MAX_HOUSEHOLD_SIZE+1],denom_household[MAX_HOUSEHOLD_SIZE+1];
double inf_household_av[MAX_HOUSEHOLD_SIZE+1][MAX_HOUSEHOLD_SIZE+1],AgeDist[NUM_AGE_GROUPS],AgeDist2[NUM_AGE_GROUPS];
double vaccdose_dist[NUM_VACCDOSE_GROUPS], vaccdosering_dist[NUM_VACCDOSE_GROUPS],vaccdosecell_dist[NUM_VACCDOSECELL_GROUPS],vaccdistance_dist[NUM_VACCDIST_GROUPS], vaccpop_dist[NUM_POP_GROUPS];
double case_household[MAX_HOUSEHOLD_SIZE+1][MAX_HOUSEHOLD_SIZE+1],case_household_av[MAX_HOUSEHOLD_SIZE+1][MAX_HOUSEHOLD_SIZE+1];
double PropPlaces[NUM_AGE_GROUPS*AGE_GROUP_WIDTH][NUM_PLACE_TYPES];
double PropPlacesC[NUM_AGE_GROUPS*AGE_GROUP_WIDTH][NUM_PLACE_TYPES],AirTravelDist[MAX_DIST];
double sinx[361],cosx[361],asin2sqx[1001];
double PeakHeightSum,PeakHeightSS,PeakTimeSum,PeakTimeSS;
bitmap_header *bmh,*bmh2,*bmh3;
void *BinFileBuf;
bin_file *BF;
unsigned char *bmf,*bm,*bmp,*bmf2,*bm2,*bmp2,*bmf3,*bm3,*bmp3;
float *bmi,*bmi2,*bmi3,*bmi4,*bmi5,*bmi6,*bmi2min,*bmi2max,*bmi2mean,*bmi3min,*bmi3max,*bmi3mean,*bmi7,*bmi7min,*bmi7max,*bmi7mean;

#ifdef WIN32_BM
HAVI avi,avi2,avi3;
HBITMAP bmpdib,bmpdib2,bmpdib3;
ULONG_PTR m_gdiplusToken;
CLSID  encoderClsid;
#endif
char OutFile[1024],OutFileBase[1024],OutDensFile[1024],SnapshotLoadFile[1024],SnapshotSaveFile[1024],RadiationFile[1024],RoadNetworkFile[1024]; //added radiation file for radiation model: ggilani 09/02/15

int ns,DoInitUpdateProbs,InterruptRun=0;
int netbuf[NUM_PLACE_TYPES_NOAIR*1000000];
unsigned int cntr=0;
int PlaceDistDistrib[NUM_PLACE_TYPES][MAX_DIST],PlaceSizeDistrib[NUM_PLACE_TYPES][MAX_PLACE_SIZE];


#ifdef FRESSCA
GEO_DISTR_NETWORK *VaccineDistributionNetwork;
int **MCellDistrIndex,noCenterMcellNum;
int GotDN;
char DistribNetworkFile[1024]; //Making this global so as to not change the syntax of functions
#endif
FILE *KMLFile,*KMLFile2;

/* int NumPC,NumPCD; */
#define MAXINTFILE 10

int main(int argc,char *argv[])
{
	char ParamFile[1024],DensityFile[1024],NetworkFile[1024],AirTravelFile[1024],SchoolFile[1024],RegDemogFile[1024],InterventionFile[MAXINTFILE][1024],PreParamFile[1024],buf[2048],*sep;
	int i,j,k,GotP,GotPP,GotO,GotD,GotL,GotS,GotA,GotScF,GotIF,Perr,cl;
	double t,t2;

	Perr=0;
	fprintf(stderr,"sizeof(int)=%i sizeof(long)=%i sizeof(float)=%i sizeof(double)=%i sizeof(unsigned short int)=%i sizeof(int *)=%i\n",sizeof(int),sizeof(long),sizeof(float),sizeof(double),sizeof(unsigned short int),sizeof(int *));
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
#ifndef FRESSCA
				fprintf(stderr,"This Version of SpatialSim doesn't have FRESSCA compiled into it\n");
#else
				GotDN=1;
				P.DoDistributionVaccination = 1;
				sscanf(&argv[i][3],"%s",DistribNetworkFile);
				fprintf(stderr,"Going to do the distribution Network\n");
#endif
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
			ERR_CRITICAL("Number of runs %i is greater than the size of the seed array %i!\n",P.NR,MAX_FIXED_SEEDS);
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
	//SaveAgeDistrib2();
	
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

#ifdef FRESSCA
	if(P.DoDistributionVaccination)
		{
		VaccineDistributionNetwork->EndGE(KMLFile);
		fclose(KMLFile);
		}
#endif
#ifdef WIN32_BM
	Gdiplus::GdiplusShutdown(m_gdiplusToken);
#endif
	fprintf(stderr,"Extinction in %i out of %i runs\n",P.NRactE,P.NRactNE+P.NRactE);
	fprintf(stderr,"Model ran in %lg seconds\n",((double) (clock()-cl))/CLOCKS_PER_SEC);
	fprintf(stderr,"Model finished\n");
}


void ReadParams(char *ParamFile,char *PreParamFile)
{
	FILE *dat,*dat2,*dat3;
	double s,t,CumAgeDist[NUM_AGE_GROUPS+1],AgeSuscScale;
	int i,j,k,nc,na;
	char buf[1024],*CountryNames[MAX_COUNTRIES],CountryNameBuf[128*MAX_COUNTRIES];
	char **AdunitNames, *AdunitNamesBuf;

	if(!(dat=fopen(ParamFile,"r"))) ERR_CRITICAL("Unable to open parameter file\n");
	dat2=fopen(PreParamFile,"r");

	AgeSuscScale=1.0;
	GetInputParameter(dat,dat2,"Update timestep","%lf",(void *) &(P.TimeStep),1,1,0);
	GetInputParameter(dat,dat2,"Sampling timestep","%lf",(void *) &(P.SampleStep),1,1,0);
	if(P.TimeStep>P.SampleStep) ERR_CRITICAL("Update step must be smaller than sampling step\n");
	t=ceil(P.SampleStep/P.TimeStep-1e-6);
	P.UpdatesPerSample = (int)t;
	P.TimeStep=P.SampleStep/t;
	P.TimeStepsPerDay=ceil(1.0/P.TimeStep-1e-6);
	P.TimeStepsPerDayInt=(int) P.TimeStepsPerDay;
	P.TimeStepsPerYear=P.TimeStepsPerDay*DAYS_PER_YEAR;
	fprintf(stderr,"Update step = %lf\nSampling step = %lf\nUpdates per sample=%i\nTimeStepsPerDay=%lf\n",P.TimeStep,P.SampleStep,P.UpdatesPerSample,P.TimeStepsPerDay);
	GetInputParameter(dat,dat2,"Sampling time","%lf",(void *) &(P.SampleTime),1,1,0);
	P.NumSamples=1+(int) ceil(P.SampleTime/P.SampleStep);
	GetInputParameter(dat,dat2,"Population size","%i",(void *) &(P.N),1,1,0);
	GetInputParameter(dat,dat2,"Number of realisations","%i",(void *) &(P.NR),1,1,0);
	if(!GetInputParameter2(dat,dat2,"Number of non-extinct realisations","%i",(void *) &(P.NRN),1,1,0)) P.NRN=P.NR;
	if(!GetInputParameter2(dat,dat2,"Maximum number of cases defining small outbreak","%i",(void *) &(P.SmallEpidemicCases),1,1,0)) P.SmallEpidemicCases=-1;
	GetInputParameter(dat,dat2,"Number of spatial cells","%i",(void *) &(P.NC),1,1,0);
	GetInputParameter(dat,dat2,"Number of micro-cells per spatial cell width","%i",(void *) &(P.NMCL),1,1,0);
	//added parameter to reset seeds after every run
	if(!GetInputParameter2(dat,dat2,"Reset seeds for every run","%i",(void *) &(P.ResetSeeds),1,1,0)) P.ResetSeeds=0;
	if (P.ResetSeeds)
	{
		if (!GetInputParameter2(dat, dat2, "Keep same seeds for every run", "%i", (void*) &(P.KeepSameSeeds), 1, 1, 0)) P.KeepSameSeeds = 0; //added this to control which seeds are used: ggilani 27/11/19
		if (!GetInputParameter2(dat, dat2, "Use fixed input seeds", "%i", (void*) &(P.DoFixedSeeds), 1, 1, 0)) P.DoFixedSeeds = 0; //added this to allow a fixed list of seeds for Janetta: ggilani 08/03/2023
	}
	if (!GetInputParameter2(dat, dat2, "Reset seeds after intervention", "%i", (void*) &(P.ResetSeedsPostIntervention), 1, 1, 0)) P.ResetSeedsPostIntervention = 0;
	if (P.ResetSeedsPostIntervention)
	{
		if (!GetInputParameter2(dat, dat2, "Time to reset seeds after intervention", "%i", (void*) & (P.TimeToResetSeeds), 1, 1, 0)) P.TimeToResetSeeds = 1e6;
	}
	if(!GetInputParameter2(dat,dat2,"Include households","%i",(void *) &(P.DoHouseholds),1,1,0)) P.DoHouseholds=0;
	if(P.DoHouseholds)
		{
		GetInputParameter(dat,dat2,"Household size distribution","%lf",(void *) P.HouseholdSizeDistrib[0],MAX_HOUSEHOLD_SIZE,1,0);
		GetInputParameter(dat,dat2,"Household attack rate","%lf",(void *) &(P.HouseholdTrans),1,1,0);
		GetInputParameter(dat,dat2,"Household transmission denominator power","%lf",(void *) &(P.HouseholdTransPow),1,1,0);
		if(!GetInputParameter2(dat,dat2,"Median household income","%lf",(void *) &(P.MedianIncome),1,1,0)) P.MedianIncome=10000;
		if(!GetInputParameter2(dat,dat2,"Household income weibull power","%lf",(void *) &(P.IncomeWeibullPower),1,1,0)) P.IncomeWeibullPower=1.5;
		if (!GetInputParameter2(dat, dat2, "Output household file", "%i", (void*)&(P.DoHouseholdOutput), 1, 1, 0)) P.DoHouseholdOutput = 0;
		}
	else
		{
		P.HouseholdTrans=0.0;
		P.HouseholdTransPow=1.0;
		P.HouseholdSizeDistrib[0][0]=1.0;
		for(i=1;i<MAX_HOUSEHOLD_SIZE;i++)
			P.HouseholdSizeDistrib[0][i]=0;
		P.MedianIncome=10000;
		P.IncomeWeibullPower=1.5;
		}
	for(i=1;i<MAX_HOUSEHOLD_SIZE;i++)
		P.HouseholdSizeDistrib[0][i]=P.HouseholdSizeDistrib[0][i]+P.HouseholdSizeDistrib[0][i-1];
	for(i=0;i<MAX_HOUSEHOLD_SIZE;i++)
		P.HouseholdDenomLookup[i]=1/pow(((double)(i+1)),P.HouseholdTransPow);
	if(!GetInputParameter2(dat,dat2,"Include administrative units within countries","%i",(void *) &(P.DoAdUnits),1,1,0)) P.DoAdUnits=0;
	if(!GetInputParameter2(dat,dat2,"Divisor for countries","%i",(void *) &(P.CountryDivisor),1,1,0)) P.CountryDivisor=1;
	if (!GetInputParameter2(dat, dat2, "Output country file", "%i", (void*)&(P.DoCountryOutput), 1, 1, 0)) P.DoCountryOutput = 0;
	if(P.DoAdUnits)
		{
		if(!(AdunitNames=(char **) malloc(3*ADUNIT_LOOKUP_SIZE*sizeof(char *)))) ERR_CRITICAL("Unable to allocate temp storage\n");
		if(!(AdunitNamesBuf=(char *) malloc(3*ADUNIT_LOOKUP_SIZE*360*sizeof(char)))) ERR_CRITICAL("Unable to allocate temp storage\n");

		for(i=0;i<ADUNIT_LOOKUP_SIZE;i++) 
			{
			P.AdunitLevel1Lookup[i]=-1;
			AdunitNames[3*i]=AdunitNamesBuf+3*i*360;
			AdunitNames[3*i+1]=AdunitNamesBuf+3*i*360+60;
			AdunitNames[3*i+2]=AdunitNamesBuf+3*i*360+160;
			}
		for(i=0;i<MAX_COUNTRIES;i++) {CountryNames[i]=CountryNameBuf+128*i;CountryNames[i][0]=0;}
		na=(GetInputParameter2(dat,dat2,"Codes and country/province names for admin units","%s",(void *) AdunitNames,3*ADUNIT_LOOKUP_SIZE,1,0))/3;
		if((na>0)&&(GetInputParameter2(dat,dat2,"Number of countries to include","%i",(void *) &nc,1,1,0)))
			{
			P.DoAdunitBoundaries=(nc>0);
			nc=abs(nc);
			P.AdunitLevel1Divisor=1;
			P.AdunitLevel1Mask=1000000000;
			GetInputParameter(dat,dat2,"List of names of countries to include","%s",(nc>1)?((void *) CountryNames):((void *)CountryNames[0]),nc,1,0);
			P.NumAdunits=0;
			for(i=0;i<na;i++)
				for(j=0;j<nc;j++)
					{
					if((AdunitNames[3*i+1][0])&&(!strcmp(AdunitNames[3*i+1],CountryNames[j]))&&(atoi(AdunitNames[3*i])!=0))
						{
						AdUnits[P.NumAdunits].id=atoi(AdunitNames[3*i]);
						P.AdunitLevel1Lookup[AdUnits[P.NumAdunits].id]=P.NumAdunits;
						if(strlen(AdunitNames[3*i+1])<100) strcpy(AdUnits[P.NumAdunits].cnt_name,AdunitNames[3*i+1]);
						if(strlen(AdunitNames[3*i+2])<200) strcpy(AdUnits[P.NumAdunits].ad_name,AdunitNames[3*i+2]);
//						fprintf(stderr,"%i %s %s ## ",AdUnits[P.NumAdunits].id,AdUnits[P.NumAdunits].cnt_name,AdUnits[P.NumAdunits].ad_name);
						P.NumAdunits++;
						}
					}
			}
		else
			{
			if(!GetInputParameter2(dat,dat2,"Number of level 1 administrative units to include","%i",(void *) &(P.NumAdunits),1,1,0)) P.NumAdunits=0;
			if(!GetInputParameter2(dat,dat2,"Divisor for level 1 administrative units","%i",(void *) &(P.AdunitLevel1Divisor),1,1,0)) P.AdunitLevel1Divisor=1;
			if(!GetInputParameter2(dat,dat2,"Mask for level 1 administrative units","%i",(void *) &(P.AdunitLevel1Mask),1,1,0)) P.AdunitLevel1Mask=1000000000;
			if(P.NumAdunits>0)
				{
				P.DoAdunitBoundaries=1;
				int AdunitList[MAX_ADUNITS];
				if(P.DoAdunitBoundaries>MAX_ADUNITS) ERR_CRITICAL("MAX_ADUNITS too small.\n");
				GetInputParameter(dat,dat2,"List of level 1 administrative units to include","%i",(void *) AdunitList,P.NumAdunits,1,0);
				for(i=0;i<P.NumAdunits;i++)
					{
					AdUnits[i].id=AdunitList[i];
					P.AdunitLevel1Lookup[(AdUnits[i].id%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor]=i;
					for(j=0;j<na;j++)
						if(atoi(AdunitNames[3*j])==AdunitList[i])
							{
							if(strlen(AdunitNames[3*j+1])<100) strcpy(AdUnits[i].cnt_name,AdunitNames[3*j+1]);
							if(strlen(AdunitNames[3*j+2])<200) strcpy(AdUnits[i].ad_name,AdunitNames[3*j+2]);
							j=na;
							}
					}
				if(!GetInputParameter2(dat,dat2,"Median income of level 1 administrative units","%lf",(void *) &(P.MedianIncomeByAdunit[0]),P.DoAdunitBoundaries,1,0))
					{
					for(i=0;i<P.NumAdunits;i++)
						P.MedianIncomeByAdunit[i]=P.MedianIncome;
					}
				}
			else
				P.DoAdunitBoundaries=0;
			free(AdunitNames);
			free(AdunitNamesBuf);
			}
		if(!GetInputParameter2(dat,dat2,"Output incidence by administrative unit","%i",(void *) &(P.DoAdunitOutput),1,1,0)) P.DoAdunitOutput=0;
		if(!GetInputParameter2(dat,dat2,"Draw administrative unit boundaries on maps","%i",(void *) &(P.DoAdunitBoundaryOutput),1,1,0)) P.DoAdunitBoundaryOutput=0;
		if(!GetInputParameter2(dat,dat2,"Correct administrative unit populations","%i",(void *) &(P.DoCorrectAdunitPop),1,1,0)) P.DoCorrectAdunitPop=0;
		if(!GetInputParameter2(dat,dat2,"Fix population size at specified value","%i",(void *) &(P.DoSpecifyPop),1,1,0)) P.DoSpecifyPop=0;
		fprintf(stderr,"Using %i administrative units\n",P.NumAdunits);
		if(!GetInputParameter2(dat,dat2,"Divisor for administrative unit codes for boundary plotting on bitmaps","%i",(void *) &(P.AdunitBitmapDivisor),1,1,0)) P.AdunitBitmapDivisor=1;
		if(!GetInputParameter2(dat,dat2,"Only output household to place distance distribution for one administrative unit","%i",(void *) &(P.DoOutputPlaceDistForOneAdunit),1,1,0)) P.DoOutputPlaceDistForOneAdunit=0;
		if(P.DoOutputPlaceDistForOneAdunit)
			{
			if(!GetInputParameter2(dat,dat2,"Administrative unit for which household to place distance distribution to be output","%i",(void *) &(P.OutputPlaceDistAdunit),1,1,0)) P.DoOutputPlaceDistForOneAdunit=0;
			}
		}
	else
		{P.DoAdunitBoundaries=P.DoAdunitBoundaryOutput=P.DoAdunitOutput=P.DoCorrectAdunitPop=P.DoSpecifyPop=0;P.AdunitLevel1Divisor=1;P.AdunitLevel1Mask=1000000000;P.AdunitBitmapDivisor=P.AdunitLevel1Divisor;}
	//added flag for outputting summary results
	if (!GetInputParameter2(dat, dat2, "Output summary results", "%i", (void*)&(P.DoSummaryOutput), 1, 1, 0)) P.DoSummaryOutput = 0;

	if(!GetInputParameter2(dat,dat2,"Include age","%i",(void *) &(P.DoAge),1,1,0)) P.DoAge=1;
	if(!P.DoAge)
		{
		P.DoAge=0;
		for(i=0;i<NUM_AGE_GROUPS;i++)
			P.PropAgeGroup[0][i]=1.0/NUM_AGE_GROUPS;
		for(i=0;i<NUM_AGE_GROUPS;i++)
			{
			P.InitialImmunity[i]=0;
			P.AgeInfectiousness[i]=P.AgeSusceptibility[i]=1;
			P.RelativeSpatialContact[i]=P.RelativeTravelRate[i]=1.0;
			}
		}
	else
		{

		if (!GetInputParameter2(dat, dat2, "Output age file", "%i", (void*)&(P.DoAgeOutput), 1, 1, 0)) P.DoAgeOutput = 0;
		if(P.DoHouseholds)
			{if(!GetInputParameter2(dat,dat2,"Initial immunity applied to all household members","%i",(void *) &(P.DoWholeHouseholdImmunity),1,1,0)) P.DoWholeHouseholdImmunity=0;}
		else
			P.DoWholeHouseholdImmunity=0;
		if(!GetInputParameter2(dat,dat2,"Initial immunity profile by age","%lf",(void *) P.InitialImmunity,NUM_AGE_GROUPS,1,0))
			for(i=0;i<NUM_AGE_GROUPS;i++)
				P.InitialImmunity[i]=0;
		if(!GetInputParameter2(dat,dat2,"Relative spatial contact rates by age","%lf",(void *) P.RelativeSpatialContact,NUM_AGE_GROUPS,1,0))
			for(i=0;i<NUM_AGE_GROUPS;i++)
				P.RelativeSpatialContact[i]=1;
		if(!GetInputParameter2(dat,dat2,"Age-dependent infectiousness","%lf",(void *) P.AgeInfectiousness,NUM_AGE_GROUPS,1,0))
			for(i=0;i<NUM_AGE_GROUPS;i++)
				P.AgeInfectiousness[i]=1.0;
		if(!GetInputParameter2(dat,dat2,"Age-dependent susceptibility","%lf",(void *) P.AgeSusceptibility,NUM_AGE_GROUPS,1,0))
			for(i=0;i<NUM_AGE_GROUPS;i++)
				P.AgeSusceptibility[i]=1.0;
		GetInputParameter(dat,dat2,"Age distribution of population","%lf",(void *) P.PropAgeGroup[0],NUM_AGE_GROUPS,1,0);
		t=0;
		for(i=0;i<NUM_AGE_GROUPS;i++)
			t+=P.PropAgeGroup[0][i];
		for(i=0;i<NUM_AGE_GROUPS;i++)
			P.PropAgeGroup[0][i]/=t;
		t=0;
		for(i=0;i<NUM_AGE_GROUPS;i++)
			if(P.AgeSusceptibility[i]>t) t=P.AgeSusceptibility[i];  //peak susc has to be 1
		for(i=0;i<NUM_AGE_GROUPS;i++)
			P.AgeSusceptibility[i]/=t;
		AgeSuscScale=t;
		if(P.DoHouseholds) P.HouseholdTrans*=AgeSuscScale;
		if(!GetInputParameter2(dat,dat2,"Relative travel rates by age","%lf",(void *) P.RelativeTravelRate,NUM_AGE_GROUPS,1,0))
			for(i=0;i<NUM_AGE_GROUPS;i++)
				P.RelativeTravelRate[i]=1;
		if(!GetInputParameter2(dat,dat2,"WAIFW matrix","%lf",(void *) P.WAIFW_Matrix,NUM_AGE_GROUPS,NUM_AGE_GROUPS,0))
			{
			for(i=0;i<NUM_AGE_GROUPS;i++)
				for(j=0;j<NUM_AGE_GROUPS;j++)
					P.WAIFW_Matrix[i][j]=1.0;
			}
		else
			{
/* WAIFW matrix needs to be scaled to have max value of 1.
1st index of matrix specifies host being infected, second the infector.
Overall age variation in infectiousness/contact rates/susceptibility should be factored
out of WAIFW_matrix and put in Age dep infectiousness/susceptibility for efficiency. */
			t=0;
			for(i=0;i<NUM_AGE_GROUPS;i++)
				for(j=0;j<NUM_AGE_GROUPS;j++)
					if(P.WAIFW_Matrix[i][j]>t) t=P.WAIFW_Matrix[i][j];
			if(t>0)
				{
				for(i=0;i<NUM_AGE_GROUPS;i++)
					for(j=0;j<NUM_AGE_GROUPS;j++)
						P.WAIFW_Matrix[i][j]/=t;
				}
			else
				{
				for(i=0;i<NUM_AGE_GROUPS;i++)
					for(j=0;j<NUM_AGE_GROUPS;j++)
						P.WAIFW_Matrix[i][j]=1.0;
				}
			}
#ifdef NEW_AGE_MODEL
		if(!GetInputParameter2(dat,dat2,"Time steps per demography update (SIR only)","%i",(void *) &(P.UpdatesPerDemogUpdate),1,1,0)) P.UpdatesPerDemogUpdate=P.UpdatesPerSample;
		if(!GetInputParameter2(dat,dat2,"Probability of death by year of age","%lf",(void *) &(P.MortalityByAge[0][0]),NUM_AGE_GROUPS*AGE_GROUP_WIDTH,1,0))
			{
			for(i=0;i<NUM_AGE_GROUPS*AGE_GROUP_WIDTH;i++)
				{
				P.MortalityByAge[0][i]=0;
				P.CumulPropDead[0][i]=0;
				}
			for(i=0;i<=1000;i++) P.InvLifeExpecDist[0][i]=NUM_AGE_GROUPS*AGE_GROUP_WIDTH-1e-6;
			P.DoDeath=0;
			}
		else
			{
			t=0;
			for(i=0;i<NUM_AGE_GROUPS*AGE_GROUP_WIDTH;i++) t+=P.MortalityByAge[0][i];
			if(t>0)
				{
				for(i=0;i<NUM_AGE_GROUPS*AGE_GROUP_WIDTH;i++) P.MortalityByAge[0][i]/=t;
				P.CumulPropDead[0][0]=0;
				for(i=1;i<=NUM_AGE_GROUPS*AGE_GROUP_WIDTH;i++)
					P.CumulPropDead[0][i]=P.CumulPropDead[0][i-1]+P.MortalityByAge[0][i-1];
				P.DoDeath=1;
				}
			else
				{
				for(i=1;i<=NUM_AGE_GROUPS*AGE_GROUP_WIDTH;i++)
					P.CumulPropDead[0][i]=0;
				P.DoDeath=0;
				}		
			P.MeanAnnualDeathRate=0;
			for(i=0;i<NUM_AGE_GROUPS*AGE_GROUP_WIDTH;i++) P.MeanAnnualDeathRate+=P.MortalityByAge[0][i]*P.PropAgeGroup[0][i/AGE_GROUP_WIDTH]/AGE_GROUP_WIDTH;
			for(i=j=0;i<1000;i++)
				{
				t=((double) i)/1000;
				while ((t>=P.CumulPropDead[0][j+1])&&(j+1<NUM_AGE_GROUPS*AGE_GROUP_WIDTH)) j++;
				if(t<P.CumulPropDead[0][j+1])
					{
					t=(((double) j)+(t-P.CumulPropDead[0][j])/(P.CumulPropDead[0][j+1]-P.CumulPropDead[0][j]));
					P.InvLifeExpecDist[0][i]=t;
					}
				else
					P.InvLifeExpecDist[0][i]=NUM_AGE_GROUPS*AGE_GROUP_WIDTH-1e-9;
				}
			P.InvLifeExpecDist[0][1000]=NUM_AGE_GROUPS*AGE_GROUP_WIDTH-1e-9;
/*			for(i=0;i<NUM_AGE_GROUPS*AGE_GROUP_WIDTH;i++)
				fprintf(stderr,"(%i:%lf:%lf) ",i,P.CumulPropDead[0][i],P.InvLifeExpecDist[0][(int) ceil(1000*P.CumulPropDead[0][i])]);
*/
			if(!GetInputParameter2(dat,dat2,"Start time for routine immunisation","%lf",(void *) &(P.RoutineImmunisationStartTime),1,1,0)) P.RoutineImmunisationStartTime=0;
			if(!GetInputParameter2(dat,dat2,"Routine immunisation effective coverage","%lf",(void *) &(P.RoutineImmunisationEffectiveCoverage),1,1,0)) P.RoutineImmunisationEffectiveCoverage=0;
			if(!GetInputParameter2(dat,dat2,"Minimum age in years for routine immunisation","%lf",(void *) &(P.RoutineImmunisationMinAge),1,1,0)) P.RoutineImmunisationMinAge=1;
			if(!GetInputParameter2(dat,dat2,"Maximum age in years for routine immunisation","%lf",(void *) &(P.RoutineImmunisationMaxAge),1,1,0)) P.RoutineImmunisationMaxAge=3;
			P.usRoutineImmunisationStartTime=(unsigned short int) (P.RoutineImmunisationStartTime*P.TimeStepsPerDay);
			P.usRoutineImmunisationMinAge=(unsigned short int) (P.RoutineImmunisationMinAge*P.TimeStepsPerYear);
			P.usRoutineImmunisationMaxAge=(unsigned short int) (P.RoutineImmunisationMaxAge*P.TimeStepsPerYear);
			if(!GetInputParameter2(dat,dat2,"Start time for SIA","%lf",(void *) &(P.SIAStartTime),1,1,0)) P.SIAStartTime=0;
			if(!GetInputParameter2(dat,dat2,"SIA effective coverage","%lf",(void *) &(P.SIAEffectiveCoverage),1,1,0)) P.SIAEffectiveCoverage=0;
			if(!GetInputParameter2(dat,dat2,"SIA campaign duration in days","%lf",(void *) &(P.SIADuration),1,1,0)) P.SIADuration=7;
			if(!GetInputParameter2(dat,dat2,"Minimum age in years for SIA","%lf",(void *) &(P.SIAMinAge),1,1,0)) P.SIAMinAge=2;
			if(!GetInputParameter2(dat,dat2,"Maximum age in years for SIA","%lf",(void *) &(P.SIAMaxAge),1,1,0)) P.SIAMaxAge=10;
			if(!GetInputParameter2(dat,dat2,"Repeat interval for SIA in years","%lf",(void *) &(P.SIARepeatInterval),1,1,0)) P.SIARepeatInterval=1e10;
			P.usSIAStartTime=(unsigned short int) (P.SIAStartTime*P.TimeStepsPerDay);
			P.usSIAMinAge=(unsigned short int) (P.SIAMinAge*P.TimeStepsPerYear);
			P.usSIAMaxAge=(unsigned short int) (P.SIAMaxAge*P.TimeStepsPerYear);
			P.usSIADuration=(unsigned short int) (P.SIADuration*P.TimeStepsPerDay);
			if(P.SampleTime/DAYS_PER_YEAR<P.SIARepeatInterval)
				P.SIARepeatInterval=P.SampleTime+1;
			else
				P.SIARepeatInterval*=DAYS_PER_YEAR;
			P.usSIARepeatInterval=(unsigned short int) (P.SIARepeatInterval*P.TimeStepsPerDay);
			if(P.DoDistributionVaccination)
				{
				if(!GetInputParameter2(dat,dat2,"Country code for vaccination distribution network (-1 for all)","%i",(void *) &(P.DistribNetCountry),1,1,0)) P.DistribNetCountry=-1;
				if(!GetInputParameter2(dat,dat2,"Routine immunisation effective coverage for countries outside modelled supply chain","%lf",(void *) &(P.RoutineCoverageNonDistrCountry),1,1,0)) P.RoutineCoverageNonDistrCountry=P.RoutineImmunisationEffectiveCoverage;
				if(!GetInputParameter2(dat,dat2,"Apply SIAs to countries outside modelled supply chain","%i",(void *) &(P.SIADoAllCountries),1,1,0)) P.SIADoAllCountries=0;
				if(!GetInputParameter2(dat,dat2,"Vaccine doses per phial (to calculate wastage)","%i",(void *) &(P.VaccDosesPerPhial),1,1,0)) P.VaccDosesPerPhial=1;
				if(!GetInputParameter2(dat,dat2,"Vaccine phial lifetime once opened in days","%lf",(void *) &(P.VaccPhialLifetime),1,1,0)) P.VaccPhialLifetime=1;
				}
			else
				{
				P.SIADoAllCountries=0;
				}
			}
		if(P.DoDeath)
			{
			fprintf(stderr,"Implementing simple host birth/death process\n");
			if(!GetInputParameter2(dat,dat2,"Start at SIR equilibrium","%i",(void *) &(P.DoInitEquilib),1,1,0)) P.DoInitEquilib=0;
			}
#else
		P.DoDeath=0;
#endif
		t=0;
		for(i=0;i<NUM_AGE_GROUPS;i++)
			t+=P.AgeInfectiousness[i]*P.PropAgeGroup[0][i];
		for(i=0;i<NUM_AGE_GROUPS;i++)
			P.AgeInfectiousness[i]/=t;
		}
	if(!GetInputParameter2(dat,dat2,"Include spatial transmission","%i",(void *) &(P.DoSpatial),1,1,0)) P.DoSpatial=1;
	GetInputParameter(dat,dat2,"Kernel type","%i",(void *) &(P.MoveKernelType),1,1,0);
	GetInputParameter(dat,dat2,"Kernel scale","%lf",(void *) &(P.MoveKernelScale),1,1,0);
	if(P.KernelOffsetScale!=1)
	{
		P.MoveKernelScale*=P.KernelOffsetScale;
	}
	if(!GetInputParameter2(dat,dat2,"Kernel 3rd param","%lf",(void *) &(P.MoveKernelP3),1,1,0)) P.MoveKernelP3=0;
	if(!GetInputParameter2(dat,dat2,"Kernel 4th param","%lf",(void *) &(P.MoveKernelP4),1,1,0)) P.MoveKernelP4=0;
	if(!GetInputParameter2(dat,dat2,"Kernel Shape","%lf",(void *) &(P.MoveKernelShape),1,1,0)) P.MoveKernelShape=1.0;
	if(P.KernelPowerScale!=1)
	{
		P.MoveKernelShape*=P.KernelPowerScale;
	}
	if(!GetInputParameter2(dat,dat2,"Airport Kernel Type","%i",(void *) &(P.AirportKernelType),1,1,0)) P.AirportKernelType=P.MoveKernelType;
	if(!GetInputParameter2(dat,dat2,"Airport Kernel Scale","%lf",(void *) &(P.AirportKernelScale),1,1,0)) P.AirportKernelScale=P.MoveKernelScale;
	if(!GetInputParameter2(dat,dat2,"Airport Kernel Shape","%lf",(void *) &(P.AirportKernelShape),1,1,0)) P.AirportKernelShape=P.MoveKernelShape;
	if(!GetInputParameter2(dat,dat2,"Airport Kernel 3rd param","%lf",(void *) &(P.AirportKernelP3),1,1,0)) P.AirportKernelP3=P.MoveKernelP3;
	if(!GetInputParameter2(dat,dat2,"Airport Kernel 4th param","%lf",(void *) &(P.AirportKernelP4),1,1,0)) P.AirportKernelP4=P.MoveKernelP4;
	if(!GetInputParameter2(dat,dat2,"Include places","%i",(void *) &(P.DoPlaces),1,1,0)) P.DoPlaces=P.PlaceTypeNum=P.DoAirports=0;
	if(P.DoPlaces)
		{
		GetInputParameter(dat,dat2,"Number of types of places","%i",(void *) &(P.PlaceTypeNum),1,1,0);
		if(P.PlaceTypeNum==0) P.DoPlaces=P.DoAirports=0;
		}
	if(P.DoPlaces)
		{
		if(P.PlaceTypeNum>NUM_PLACE_TYPES) ERR_CRITICAL("Too many place types\n");
		GetInputParameter(dat,dat2,"Minimum age for age group 1 in place types","%i",(void *) P.PlaceTypeAgeMin,P.PlaceTypeNum,1,0);
		GetInputParameter(dat,dat2,"Maximum age for age group 1 in place types","%i",(void *) P.PlaceTypeAgeMax,P.PlaceTypeNum,1,0);
		GetInputParameter(dat,dat2,"Proportion of age group 1 in place types","%lf",(void *) &(P.PlaceTypePropAgeGroup),P.PlaceTypeNum,1,0);
		if(!GetInputParameter2(dat,dat2,"Proportion of age group 2 in place types","%lf",(void *) &(P.PlaceTypePropAgeGroup2),P.PlaceTypeNum,1,0))
			{
			for(i=0;i<NUM_PLACE_TYPES;i++)
				{
				P.PlaceTypePropAgeGroup2[i]=0;
				P.PlaceTypeAgeMin2[i]=0;
				P.PlaceTypeAgeMax2[i]=1000;
				}
			}
		else
			{
			GetInputParameter(dat,dat2,"Minimum age for age group 2 in place types","%i",(void *) P.PlaceTypeAgeMin2,P.PlaceTypeNum,1,0);
			GetInputParameter(dat,dat2,"Maximum age for age group 2 in place types","%i",(void *) P.PlaceTypeAgeMax2,P.PlaceTypeNum,1,0);
			}
		if(!GetInputParameter2(dat,dat2,"Proportion of age group 3 in place types","%lf",(void *) &(P.PlaceTypePropAgeGroup3),P.PlaceTypeNum,1,0))
			{
			for(i=0;i<NUM_PLACE_TYPES;i++)
				{
				P.PlaceTypePropAgeGroup3[i]=0;
				P.PlaceTypeAgeMin3[i]=0;
				P.PlaceTypeAgeMax3[i]=1000;
				}
			}
		else
			{
			GetInputParameter(dat,dat2,"Minimum age for age group 3 in place types","%i",(void *) P.PlaceTypeAgeMin3,P.PlaceTypeNum,1,0);
			GetInputParameter(dat,dat2,"Maximum age for age group 3 in place types","%i",(void *) P.PlaceTypeAgeMax3,P.PlaceTypeNum,1,0);
			}		
		if(!GetInputParameter2(dat,dat2,"Kernel shape params for place types","%lf",(void *) &(P.PlaceTypeKernelShape),P.PlaceTypeNum,1,0))
			{
			for(i=0;i<NUM_PLACE_TYPES;i++)
				{
				P.PlaceTypeKernelShape[i]=P.MoveKernelShape;
				P.PlaceTypeKernelScale[i]=P.MoveKernelScale;
				}
			}
		else
			GetInputParameter(dat,dat2,"Kernel scale params for place types","%lf",(void *) &(P.PlaceTypeKernelScale),P.PlaceTypeNum,1,0);
		if(!GetInputParameter2(dat,dat2,"Kernel 3rd param for place types","%lf",(void *) &(P.PlaceTypeKernelP3),P.PlaceTypeNum,1,0))
			{
			for(i=0;i<NUM_PLACE_TYPES;i++)
				{
				P.PlaceTypeKernelP3[i]=P.MoveKernelP3;
				P.PlaceTypeKernelP4[i]=P.MoveKernelP4;
				}
			}
		else
			GetInputParameter(dat,dat2,"Kernel 4th param for place types","%lf",(void *) &(P.PlaceTypeKernelP4),P.PlaceTypeNum,1,0);
		if(!GetInputParameter2(dat,dat2,"Number of closest places people pick from (0=all) for place types","%i",(void *) &(P.PlaceTypeNearestNeighb),P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++)
				P.PlaceTypeNearestNeighb[i]=0;
		if(P.DoAdUnits)
			{
			if(!GetInputParameter2(dat,dat2,"Degree to which crossing administrative unit boundaries to go to places is inhibited","%lf",(void *) &(P.InhibitInterAdunitPlaceAssignment),P.PlaceTypeNum,1,0))
				for(i=0;i<NUM_PLACE_TYPES;i++)
					P.InhibitInterAdunitPlaceAssignment[i]=0;
			}
		if(P.DoHouseholds)
			{
			if(!GetInputParameter2(dat,dat2,"Do place type members have household contacts","%i",(void *) &(P.PlaceTypeDoHousehold),P.PlaceTypeNum,1,0))
				for(i=0;i<NUM_PLACE_TYPES;i++)
					P.PlaceTypeDoHousehold[i]=1;
			}
		//add a condition about household and place membership overlapping - ggilani 13/02/17
		if(P.DoHouseholds)
		{
			if(!GetInputParameter2(dat,dat2,"Do place and household membership overlap","%i",(void *) &(P.PlaceHouseholdOverlap),1,1,0))
			{
				P.PlaceHouseholdOverlap=0;
			}
			if (P.PlaceHouseholdOverlap)
			{
				if (!GetInputParameter2(dat, dat2, "Divisor for households into places", "%lf", (void*)&(P.PlaceHouseholdDivisor), 1, 1, 0)) P.PlaceHouseholdDivisor = 2.0;
			}
		}
		if(!GetInputParameter2(dat,dat2,"Include air travel","%i",(void *) &(P.DoAirports),1,1,0)) P.DoAirports=0;
#ifdef FAST_US
		P.DoAirports=0;
#endif
		if(P.DoAirports)
			{
			if(!GetInputParameter2(dat,dat2,"Scaling factor for input file to convert to daily traffic","%lf",(void *) &(P.AirportTrafficScale),1,1,0)) P.AirportTrafficScale=1.0;
			if(!GetInputParameter2(dat,dat2,"Proportion of hotel attendees who are local","%lf",(void *) &(P.HotelPropLocal),1,1,0)) P.HotelPropLocal=0;
			if(!GetInputParameter2(dat,dat2,"Distribution of duration of air journeys","%lf",(void *) &(P.JourneyDurationDistrib),MAX_TRAVEL_TIME,1,0))
				{
				P.JourneyDurationDistrib[0]=1;
				for(i=0;i<MAX_TRAVEL_TIME;i++)
					P.JourneyDurationDistrib[i]=0;
				}
			if(!GetInputParameter2(dat,dat2,"Distribution of duration of local journeys","%lf",(void *) &(P.LocalJourneyDurationDistrib),MAX_TRAVEL_TIME,1,0))
				{
				P.LocalJourneyDurationDistrib[0]=1;
				for(i=0;i<MAX_TRAVEL_TIME;i++)
					P.LocalJourneyDurationDistrib[i]=0;
				}
			P.MeanJourneyTime=P.MeanLocalJourneyTime=0;
			for(i=0;i<MAX_TRAVEL_TIME;i++)
				{
				P.MeanJourneyTime+=((double) (i))*P.JourneyDurationDistrib[i];
				P.MeanLocalJourneyTime+=((double) (i))*P.LocalJourneyDurationDistrib[i];
				}
			fprintf(stderr,"Mean duration of local journeys = %lf days\n",P.MeanLocalJourneyTime);
			for(i=1;i<MAX_TRAVEL_TIME;i++)
				{
				P.JourneyDurationDistrib[i]+=P.JourneyDurationDistrib[i-1];
				P.LocalJourneyDurationDistrib[i]+=P.LocalJourneyDurationDistrib[i-1];
				}
			for(i=j=0;i<=1024;i++)
				{
				s=((double) i)/1024;
				while(P.JourneyDurationDistrib[j]<s)j++;
				P.InvJourneyDurationDistrib[i]=j;
				}
			for(i=j=0;i<=1024;i++)
				{
				s=((double) i)/1024;
				while(P.LocalJourneyDurationDistrib[j]<s)j++;
				P.InvLocalJourneyDurationDistrib[i]=j;
				}
			}
		GetInputParameter(dat,dat2,"Mean size of place types","%lf",(void *) P.PlaceTypeMeanSize,P.PlaceTypeNum,1,0);
		GetInputParameter(dat,dat2,"Param 1 of place group size distribution","%lf",(void *) P.PlaceTypeGroupSizeParam1,P.PlaceTypeNum,1,0);
		if(!GetInputParameter2(dat,dat2,"Power of place size distribution","%lf",(void *) P.PlaceTypeSizePower,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++)
				P.PlaceTypeSizePower[i]=0;
		//added to enable lognormal distribution - ggilani 09/02/17
		if(!GetInputParameter2(dat,dat2,"Standard deviation of place size distribution","%lf",(void *) P.PlaceTypeSizeSD,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++)
				P.PlaceTypeSizeSD[i]=0;
		if(!GetInputParameter2(dat,dat2,"Offset of place size distribution","%lf",(void *) P.PlaceTypeSizeOffset,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++)
				P.PlaceTypeSizeOffset[i]=0;
		if(!GetInputParameter2(dat,dat2,"Maximum of place size distribution","%lf",(void *) P.PlaceTypeSizeMax,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++)
				P.PlaceTypeSizeMax[i]=1e20;
		GetInputParameter(dat,dat2,"Proportion of between group place links","%lf",(void *) P.PlaceTypePropBetweenGroupLinks,P.PlaceTypeNum,1,0);
		if(!GetInputParameter2(dat,dat2,"Kernel type for place types","%i",(void *) P.PlaceTypeKernelType,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++)
				P.PlaceTypeKernelType[i]=P.MoveKernelType;
		GetInputParameter(dat,dat2,"Place overlap matrix","%lf",(void *) P.PlaceExclusivityMatrix,P.PlaceTypeNum*P.PlaceTypeNum,1,0); //changed from P.PlaceTypeNum,P.PlaceTypeNum,0);
/* Note P.PlaceExclusivityMatrix not used at present - places assumed exclusive (each person belongs to 0 or 1 place) */
		GetInputParameter(dat,dat2,"Relative transmission rates for place types","%lf",(void *) P.PlaceTypeTrans,P.PlaceTypeNum,1,0);
		for(i=0;i<P.PlaceTypeNum;i++) P.PlaceTypeTrans[i]*=AgeSuscScale;
		}
    if(!GetInputParameter2(dat,dat2,"Daily seasonality coefficients","%lf",(void *) P.Seasonality,DAYS_PER_YEAR,1,0)) 
		{
		P.DoSeasonality=0;
		for(i=0;i<DAYS_PER_YEAR;i++)
		  P.Seasonality[i]=1;
		}
	  else
		{
		P.DoSeasonality=1;
		s=0;
		for(i=0;i<DAYS_PER_YEAR;i++)
		  s+=P.Seasonality[i];
		s+=1e-20;
		s/=DAYS_PER_YEAR;
		for(i=0;i<DAYS_PER_YEAR;i++)
		  P.Seasonality[i]/=s;
		}
	if(!GetInputParameter2(dat,dat2,"Mutation probability per infection of step increase in infectiousness","%lf",(void *) &(P.EvolInfectMutationRate),1,1,0)) P.EvolInfectMutationRate=0;
	if(!GetInputParameter2(dat,dat2,"Step size of increase in infectiousness due to mutation","%lf",(void *) &(P.EvolInfectStep),1,1,0)) P.EvolInfectStep=0;
	if(!GetInputParameter2(dat,dat2,"Maximum increase in infectiousness due to mutation","%lf",(void *) &(P.EvolInfectMax),1,1,0)) P.EvolInfectMax=0;

	if(!GetInputParameter2(dat,dat2,"Number of seed locations","%i",(void *) &(P.NumSeedLocations),1,1,0)) P.NumSeedLocations=1;
	if(P.NumSeedLocations>MAX_NUM_SEED_LOCATIONS)
		{
		fprintf(stderr,"Too many seed locations\n");
		P.NumSeedLocations=MAX_NUM_SEED_LOCATIONS;
		}
	GetInputParameter(dat,dat2,"Initial number of infecteds","%i",(void *) P.NumInitialInfections,P.NumSeedLocations,1,0);
	GetInputParameter(dat,dat2,"Location of initial infecteds","%lf",(void *) &(P.LocationInitialInfection[0][0]),P.NumSeedLocations*2,1,0);
	if(!GetInputParameter2(dat,dat2,"Minimum population in microcell of initial infection","%i",(void *) &(P.MinPopDensForInitialInfection),1,1,0)) P.MinPopDensForInitialInfection=0;
	GetInputParameter(dat,dat2,"Maximum population in microcell of initial infection","%i",(void *) &(P.MaxPopDensForInitialInfection),1,1,0);
	GetInputParameter(dat,dat2,"Randomise initial infection location","%i",(void *) &(P.DoRandomInitialInfectionLoc),1,1,0);
	GetInputParameter(dat,dat2,"All initial infections located in same microcell","%i",(void *) &(P.DoAllInitialInfectioninSameLoc),1,1,0);
	if(P.DoAdUnits)
		{
		if(!GetInputParameter2(dat,dat2,"Administrative unit to seed initial infection into","%i",(void *) &(P.InitialInfectionsAdminUnit[0]),P.NumSeedLocations,1,0)) 
			for(i=0;i<P.NumSeedLocations;i++) P.InitialInfectionsAdminUnit[i]=0;
		}
	else
		{for(i=0;i<P.NumSeedLocations;i++) P.InitialInfectionsAdminUnit[i]=0;}
	GetInputParameter(dat,dat2,"Reproduction number","%lf",(void *) &(P.R0),1,1,0);
	GetInputParameter(dat,dat2,"Infectious period","%lf",(void *) &(P.InfectiousPeriod),1,1,0);
	if(!GetInputParameter2(dat,dat2,"Assume SIS model","%i",(void *) &(P.DoSIS),1,1,0)) P.DoSIS=0;
	if(!GetInputParameter2(dat,dat2,"SD of individual variation in infectiousness","%lf",(void *) &(P.InfectiousnessSD),1,1,0)) P.InfectiousnessSD=0;
	if(P.InfectiousnessSD>0)
	{
		P.InfectiousnessGamA=P.InfectiousnessGamR=1/(P.InfectiousnessSD*P.InfectiousnessSD);
		//GetInputParameter(dat,dat2,"Infectiousness Beta A","%lf",(void *) &(P.InfectiousnessBetaA),1,1,0);
		//GetInputParameter(dat,dat2,"Infectiousness Beta B","%lf",(void *) &(P.InfectiousnessBetaB),1,1,0);
	}
	if(P.DoSIS)
		{
		if(!GetInputParameter2(dat,dat2,"Factor by which susceptibility is multiplied per infection","%lf",(void *) &(P.SuscReductionFactorPerInfection),1,1,0)) P.SuscReductionFactorPerInfection=1;
		}
	else
		P.SuscReductionFactorPerInfection=0;

	if (!GetInputParameter2(dat, dat2, "Output R0 file", "%i", (void*)&(P.DoROutput), 1, 1, 0)) P.DoROutput = 0;
	if (!GetInputParameter2(dat, dat2, "Output infection type file", "%i", (void*)&(P.DoInftypeOutput), 1, 1, 0)) P.DoInftypeOutput = 0;
	
	if(!GetInputParameter2(dat,dat2,"Model time varying infectiousness","%i",(void *) &(P.DoInfectiousnessProfile),1,1,0)) P.DoInfectiousnessProfile=0;
	if(!GetInputParameter2(dat,dat2,"Power of scaling of spatial R0 with density","%lf",(void *) &(P.R0DensityScalePower),1,1,0)) P.R0DensityScalePower=0;
	if(P.DoInfectiousnessProfile)
		{
		if(!GetInputParameter2(dat,dat2,"Infectiousness profile","%lf",(void *) P.infectious_prof,INFPROF_RES,1,0))
			{
			for(i=0;i<INFPROF_RES;i++)
				P.infectious_prof[i]=1;
			}
		k=(int) ceil(P.InfectiousPeriod/P.TimeStep);
		if(k>=MAX_INFECTIOUS_STEPS) ERR_CRITICAL("MAX_INFECTIOUS_STEPS not big enough\n");
		s=0;
		P.infectious_prof[INFPROF_RES]=0;
		for(i=0;i<k;i++)
			{
			t=(((double) i)*P.TimeStep/P.InfectiousPeriod*INFPROF_RES);
			j=(int) t;
			t-=(double) j;
			if(j<INFPROF_RES)
				s+=(P.infectiousness[i]=P.infectious_prof[j]*(1-t)+P.infectious_prof[j+1]*t);
			else
				s+=(P.infectiousness[i]=P.infectious_prof[INFPROF_RES]);
			}
		s/=((double) k);
		for(i=0;i<=k;i++) P.infectiousness[i]/=s;
		for(i=0;i<=CDF_RES;i++) P.infectious_icdf[i]=exp(-1.0);
		}
	else
		{
		if(!GetInputParameter2(dat,dat2,"Infectious period inverse CDF","%lf",(void *) P.infectious_icdf,CDF_RES+1,1,0))
			{
			P.infectious_icdf[CDF_RES]=100;
			for(i=0;i<CDF_RES;i++)
				P.infectious_icdf[i]=-log(1-((double)i)/CDF_RES);
			}
		k=(int) ceil(P.InfectiousPeriod*P.infectious_icdf[CDF_RES]/P.TimeStep);
		if(k>=MAX_INFECTIOUS_STEPS) ERR_CRITICAL("MAX_INFECTIOUS_STEPS not big enough\n");
		for(i=0;i<k;i++) P.infectiousness[i]=1.0;
		P.infectiousness[k]=0;
		for(i=0;i<=CDF_RES;i++) P.infectious_icdf[i]=exp(-P.infectious_icdf[i]);
		}
	//if(!GetInputParameter2(dat,dat2,"Proportion of cases hospitalised","%lf",(void *) &(P.ProportionHospitalised),1,1,0)) P.ProportionHospitalised=0;
	if(!GetInputParameter2(dat,dat2,"Include mortality","%i",(void *) &(P.DoMortality),1,1,0)) P.DoMortality=0;
	if(P.DoMortality) //changes to mortality to allow for probability of death/recovery to be decided when leaving the infectious class: ggilani - 23/10/14
	{
		if (!GetInputParameter2(dat, dat2, "Include event-dependent mortality", "%i", (void*)&(P.DoEventMortality), 1, 1, 0)) P.DoEventMortality = 0;
		if(P.DoEventMortality)
		{
			if (!GetInputParameter2(dat, dat2, "Amplitude of recovery function", "%lf", (void*)&(P.RecoveryAmp), 1, 1, 0)) P.RecoveryAmp = 0;
			if (!GetInputParameter2(dat, dat2, "Weibull shape of recovery function", "%lf", (void*)&(P.RecoveryShape), 1, 1, 0)) P.RecoveryShape = 0;
			if (!GetInputParameter2(dat, dat2, "Weibull scale of recovery function", "%lf", (void*)&(P.RecoveryScale), 1, 1, 0)) P.RecoveryScale = 0;
			if (P.RecoveryAmp > 0)
			{
				for (i = 0; i < RECOVERY_RES; i++)	P.RecoveryProb[i] = P.RecoveryAmp * (WeibullCDF(i, P.RecoveryShape, P.RecoveryScale) - WeibullCDF(i - 1, P.RecoveryShape, P.RecoveryScale));
			}
			else
			{
				for (i = 0; i < RECOVERY_RES; i++)	P.RecoveryProb[i] = 1.0; //i.e. if we don't provide any information about the shape of the recovery function, we'll assume everyone recovers
			}
		}
		else
		{
			//These first parameters are scalings that are used for both age-dependent and age-independent mortality
			if (!GetInputParameter2(dat, dat2, "Relative mortality of vaccinated cases dying", "%lf", (void*)&(P.DiseaseMortalityVacc), 1, 1, 0)) P.DiseaseMortalityVacc = 1;
			if (!GetInputParameter2(dat, dat2, "Relative duration of infectiousness for dying cases", "%lf", (void*)&(P.LethalInfectiousPeriod), 1, 1, 0)) P.LethalInfectiousPeriod = 1; // P.LethalInfectiousPeriod = P.InfectiousPeriod;
			
			if (!GetInputParameter2(dat, dat2, "Include age-dependent mortality", "%i", (void*)&(P.DoAgeMortality), 1, 1, 0)) P.DoAgeMortality = 0;
			if (P.DoAgeMortality)
			{
				if (!GetInputParameter2(dat, dat2, "Mortality by age", "%lf", (void*)P.AgeMortality, NUM_AGE_GROUPS, 1, 0));
			}
			else
			{
				if (!GetInputParameter2(dat, dat2, "Proportion of cases dying", "%lf", (void*)&(P.DiseaseMortality), 1, 1, 0)) P.DiseaseMortality = 0;
			}
		}
	}
	if(!GetInputParameter2(dat,dat2,"Include funeral transmission","%i",(void *) &(P.DoFuneralTransmission),1,1,0)) P.DoFuneralTransmission=0;
	if(P.DoFuneralTransmission) //Funeral transmission parameters: ggilani - 26/10/14
	{
		if(!GetInputParameter2(dat,dat2,"Duration of funeral transmission","%lf",(void *) &(P.FuneralTransmissionDuration),1,1,0)) P.FuneralTransmissionDuration=0;
		if(!GetInputParameter2(dat,dat2,"Relative infectiousness of a funeral transmission","%lf",(void *) &(P.RelativeInfectiousnessFuneral),1,1,0)) P.RelativeInfectiousnessFuneral=1;
		if (!GetInputParameter2(dat, dat2, "Safe burial trigger incidence per cell", "%lf", (void*)&(P.FuneralControlCellIncThresh), 1, 1, 0)) P.FuneralControlCellIncThresh = 1000000000;
		if (!GetInputParameter2(dat, dat2, "Safe burial start time", "%lf", (void*)&(P.FuneralControlTimeStartBase), 1, 1, 0)) P.FuneralControlTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
		if (!GetInputParameter2(dat, dat2, "Relative infectiousness of a safe burial", "%lf", (void*)&(P.RelInfSafeFuneral), 1, 1, 0)) P.RelInfSafeFuneral = 1;
		if (!GetInputParameter2(dat, dat2, "Proportion of burials conducted safely", "%lf", (void*)&(P.ProportionSafeFuneral), 1, 1, 0)) P.ProportionSafeFuneral = 1;
		if (!GetInputParameter2(dat, dat2, "Daily burials per admin unit", "%i", (void*)&(P.AdunitSDBCapacity), 1, 1, 0)) P.AdunitSDBCapacity = 0;
		if (!GetInputParameter2(dat, dat2, "Delay to increase burial capacity", "%lf", (void*)&(P.DelayToSDB), 1, 1, 0)) P.DelayToSDB = 0;
		if (!GetInputParameter2(dat, dat2, "Capacity when burial capacity increases", "%lf", (void*)&(P.CapacityToMoreSDB), 1, 1, 0)) P.CapacityToMoreSDB = 1;
		if (!GetInputParameter2(dat, dat2, "Increase in burial capacity", "%i", (void*)&(P.incCapacitySDB), 1, 1, 0)) P.incCapacitySDB = 1;
		
		//if(!GetInputParameter2(dat,dat2,"Funeral controls by admin unit","%i",(void *) &(P.DoFuneralByAdUnit),1,1,0)) P.DoFuneralByAdUnit=0;
		//if((P.DoFuneralByAdUnit)&&(P.DoAdUnits))
		//{
		//	//Set up arrays to temporarily store funeral parameters per admin unit
		//	int AdunitCaseDetectPreFuneralControl[MAX_ADUNITS];
		//	double AdunitDelayToFuneralControl[MAX_ADUNITS];
		//	double AdunitInitPropSafeFunerals[MAX_ADUNITS];
		//	double AdunitSecPropSafeFunerals[MAX_ADUNITS];
		//	double AdunitInitRelInfSafeFuneral[MAX_ADUNITS];
		//	double AdunitSecRelInfSafeFuneral[MAX_ADUNITS];
		//	double AdunitFuneralTimeToControl[MAX_ADUNITS];
		//	double AdunitFuneralControlStart[MAX_ADUNITS];
		//	double AdunitFuneralControlEnd[MAX_ADUNITS];
		//
		//	if(!GetInputParameter2(dat,dat2,"Number of detected cases before funeral controls implementation per admin unit","%i",(void *) AdunitCaseDetectPreFuneralControl,P.NumAdunits,1,0)) for(i=0;i<P.NumAdunits;i++) AdunitCaseDetectPreFuneralControl[i]=0;
		//	if(!GetInputParameter2(dat,dat2,"Delay between detection and start of funeral controls per admin unit","%lf",(void *) AdunitDelayToFuneralControl,P.NumAdunits,1,0)) for(i=0;i<P.NumAdunits;i++) AdunitDelayToFuneralControl[i]=0;
		//	if(!GetInputParameter2(dat,dat2,"Initial proportion of safe funerals per admin unit","%lf",(void *) AdunitInitPropSafeFunerals,P.NumAdunits,1,0)) for(i=0;i<P.NumAdunits;i++) AdunitInitPropSafeFunerals[i]=0;
		//	if(!GetInputParameter2(dat,dat2,"Secondary proportion of safe funerals per admin unit","%lf",(void *) AdunitSecPropSafeFunerals,P.NumAdunits,1,0)) for(i=0;i<P.NumAdunits;i++) AdunitSecPropSafeFunerals[i]=0;
		//	if(!GetInputParameter2(dat,dat2,"Initial relative infectiousness of a safe funeral per admin unit","%lf",(void *) AdunitInitRelInfSafeFuneral,P.NumAdunits,1,0)) for(i=0;i<P.NumAdunits;i++) AdunitInitRelInfSafeFuneral[i]=1;
		//	if(!GetInputParameter2(dat,dat2,"Secondary relative infectiousness of a safe funeral per admin unit","%lf",(void *) AdunitSecRelInfSafeFuneral,P.NumAdunits,1,0)) for(i=0;i<P.NumAdunits;i++) AdunitSecRelInfSafeFuneral[i]=1;
		//	if(!GetInputParameter2(dat,dat2,"Time taken to fully implement funeral controls per admin unit","%lf",(void *) AdunitFuneralTimeToControl,P.NumAdunits,1,0)) for(i=0;i<P.NumAdunits;i++) AdunitFuneralTimeToControl[i]=0;
		//	if(!GetInputParameter2(dat,dat2,"Start time of funeral controls per admin unit","%lf",(void *) AdunitFuneralControlStart,P.NumAdunits,1,0)) for(i=0;i<P.NumAdunits;i++) AdunitFuneralControlStart[i]=USHRT_MAX-1;
		//	if(!GetInputParameter2(dat,dat2,"End time of funeral controls per admin unit","%lf",(void *) AdunitFuneralControlEnd,P.NumAdunits,1,0)) for(i=0;i<P.NumAdunits;i++) AdunitFuneralControlEnd[i]=USHRT_MAX-1;
		//	//Transfer funeral parameter per admin unit to admin unit structs
		//	for(i=0;i<P.NumAdunits;i++)
		//	{
		//		AdUnits[i].caseDetPreFuneralControl=AdunitCaseDetectPreFuneralControl[i];
		//		AdUnits[i].delayDetFuneralControl=AdunitDelayToFuneralControl[i];
		//		AdUnits[i].initPropSafeFunerals=AdunitInitPropSafeFunerals[i];
		//		AdUnits[i].secPropSafeFunerals=AdunitSecPropSafeFunerals[i];
		//		AdUnits[i].initRelInfSafeFuneral=AdunitInitRelInfSafeFuneral[i];
		//		AdUnits[i].secRelInfSafeFuneral=AdunitSecRelInfSafeFuneral[i];
		//		AdUnits[i].timeToSafeFuneral=AdunitFuneralTimeToControl[i];
		//		AdUnits[i].startFuneralControl=AdunitFuneralControlStart[i];
		//		AdUnits[i].endFuneralControl=AdunitFuneralControlEnd[i];
		//	}
		//}
		//else
		//{
		//if(!GetInputParameter2(dat,dat2,"Number of detected cases before funeral controls implementation","%i",(void *) &(P.CasesDetPreFuneralControl),1,1,0)) P.CasesDetPreFuneralControl=0;
		//if(!GetInputParameter2(dat,dat2,"Delay between detection and start of funeral controls","%lf",(void *) &(P.DelayDetectFuneralControl),1,1,0)) P.DelayDetectFuneralControl=0;
		//if(!GetInputParameter2(dat,dat2,"Initial proportion of funerals considered safe","%lf",(void *) &(P.InitProportionSafeFuneral),1,1,0)) P.InitProportionSafeFuneral=0;
		//if(!GetInputParameter2(dat,dat2,"Secondary proportion of funerals considered safe","%lf",(void *) &(P.SecProportionSafeFuneral),1,1,0)) P.SecProportionSafeFuneral=0;
		//if(!GetInputParameter2(dat,dat2,"Initial relative infectiousness of a safe funeral compared to an unsafe funeral","%lf",(void *) &(P.InitRelInfSafeFuneral),1,1,0)) P.InitRelInfSafeFuneral=1;
		//if(!GetInputParameter2(dat,dat2,"Secondary relative infectiousness of a safe funeral compared to an unsafe funeral","%lf",(void *) &(P.SecRelInfSafeFuneral),1,1,0)) P.SecRelInfSafeFuneral=1;
		//if(!GetInputParameter2(dat,dat2,"Time taken to fully implement funeral controls","%lf",(void *) &(P.TimeToSafeFuneral),1,1,0)) P.TimeToSafeFuneral=0;
		//if(!GetInputParameter2(dat,dat2,"Start time of funeral controls","%lf",(void *) &(P.FuneralControlStartTime),1,1,0)) P.FuneralControlStartTime=1e10;
		//	if(!GetInputParameter2(dat,dat2,"End time of funeral controls","%lf",(void *) &(P.FuneralControlEndTime),1,1,0)) P.FuneralControlEndTime=1e10;
		//}
	}
	if(!GetInputParameter2(dat,dat2,"Include hospitalisation","%i",(void *) &(P.DoHospitalisation),1,1,0)) P.DoHospitalisation=0; //Hospitalisation parameters: ggilani 28/10/14
	if(P.DoHospitalisation)
	{
		
		//leave the terminology the same at the moment but consider making these more general. These relate to both hospitals and ETUs
		if(!GetInputParameter2(dat,dat2,"Mean time to hospitalisation","%lf",(void *) &(P.HospitalisationTime),1,1,0)) P.HospitalisationTime=0;
		if(!GetInputParameter2(dat,dat2,"Hospital waiting time","%lf",(void *) &(P.HospWaitingTime),1,1,0)) P.HospWaitingTime=0.25; //To ensure acceptance 
		if(!GetInputParameter2(dat,dat2,"Time to hospitalisation inverse CDF","%lf",(void *) P.hospital_icdf,CDF_RES+1,1,0))
			{
			P.hospital_icdf[CDF_RES]=1e10; 
			for(i=0;i<CDF_RES;i++)
				P.hospital_icdf[i]=-log(1-((double)i)/CDF_RES);
			}
		for(i=0;i<=CDF_RES;i++)
			P.hospital_icdf[i]=exp(-P.hospital_icdf[i]);
		if(!GetInputParameter2(dat,dat2,"Number of times to hospitalisation","%i",(void *) &(P.NMeanTimeToHosp),1,1,0)) P.NMeanTimeToHosp=0;
		if(P.NMeanTimeToHosp>0)
		{
			if(P.NMeanTimeToHosp>=MAX_CHANGE_POINTS) ERR_CRITICAL("MAX_CHANGE_POINTS too small\n");
			GetInputParameter(dat,dat2,"Change points times to hospitalisation","%lf",(void *) P.ChangePointMeanTimeToHosp,P.NMeanTimeToHosp,1,0);
			GetInputParameter(dat,dat2,"Times to hospitalisation","%lf",(void *) P.MeanTimeToHosp,P.NMeanTimeToHosp,1,0);
		}
		P.CurrIndMeanTimeToHosp=0;

		//
		GetInputParameter(dat, dat2, "Proportion of cases seeking care before outbreak declared", "%lf", (void*) &P.PropHospSeekPreOutbreak, 1, 1, 0);
		if (!GetInputParameter2(dat, dat2, "Relative change in care seeking behaviour after outbreak declared", "%lf", (void*)&P.RelChangeHospSeekPostOutbreak, 1, 1, 0)) P.RelChangeHospSeekPostOutbreak = 1;

		// Currently commented this out in order to make 
		//if(!GetInputParameter2(dat,dat2,"Number of hospital beds change points","%i",(void *) &(P.NHospBeds),1,1,0)) P.NHospBeds=0;
		//if(P.NHospBeds>0)
		//{
		//	if(P.NHospBeds>=MAX_CHANGE_POINTS) ERR_CRITICAL("MAX_CHANGE_POINTS too small\n");
		//	GetInputParameter(dat,dat2,"Change points hospital beds","%lf",(void *) P.ChangePointHospBeds,P.NHospBeds,1,0);
		//	GetInputParameter(dat,dat2,"Number of hospital beds","%i",(void *) P.HospBeds,P.NHospBeds,1,0);
		//}
		//P.CurrIndHospBeds=0;

		//if we want to model ETUs on an admin unit level
		if(P.DoAdUnits)
		{
			int AdunitETUCapacity[MAX_ADUNITS];
			double AdunitETUTime[MAX_ADUNITS]; //time when beds become available in each admin unit
			if(!GetInputParameter2(dat,dat2,"Do ETUs by admin unit","%i",(void *) &(P.DoETUByAdUnit),1,1,0)) P.DoETUByAdUnit=0;
			
			if(P.DoETUByAdUnit)
			{
				if(!GetInputParameter2(dat,dat2,"Initial bed capacity per ETU admin unit","%i",(void *) AdunitETUCapacity,P.NumAdunits,1,0))
				{
					for(i=0;i<P.NumAdunits;i++) AdunitETUCapacity[i]=0;
				}
				if(!GetInputParameter2(dat,dat2,"Time ETU beds available per admin unit","%lf",(void *) AdunitETUTime,P.NumAdunits,1,0))
				{
					for(i=0;i<P.NumAdunits;i++) AdunitETUTime[i]=1e10;
				}
				
				//save these to admin unit variables
				for(i=0;i<P.NumAdunits;i++)
				{
					//AdUnits[i].totalBeds=AdunitHospCapacity[i];
					AdUnits[i].totalETUBeds = 0;
					AdUnits[i].initialETUBeds=AdunitETUCapacity[i];
					AdUnits[i].timeETUBedsAvailable=AdunitETUTime[i];
					AdUnits[i].currentETUBeds=0;
					AdUnits[i].nextETUBeds=0;
					AdUnits[i].nextTimeToETUBeds=1e10;
					AdUnits[i].flagETUBeds=0;
				}
				
				if(!GetInputParameter2(dat, dat2, "Do reactive ETU bed allocation", "%i", (void *) &(P.DoReactETUBeds), 1, 1, 0)) P.DoReactETUBeds = 0;

				if(P.DoReactETUBeds)
				{
					// for reactive provisioning of beds
					if (!GetInputParameter2(dat, dat2, "Reactive ETU bed allocation trigger incidence per cell", "%lf", (void*)&(P.ETUCellIncThresh), 1, 1, 0)) P.ETUCellIncThresh = 1000000000;
					if (!GetInputParameter2(dat, dat2, "Reactive ETU bed allocation start time", "%lf", (void*)&(P.ETUTimeStartBase), 1, 1, 0)) P.ETUTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
					//if (!GetInputParameter2(dat, dat2, "Time to start reactive bed policy", "%i", (void *) &(P.StartTimeReactiveBeds), 1, 1, 0)) P.StartTimeReactiveBeds= 1e6;
					if (!GetInputParameter2(dat, dat2, "Initial number of detected cases to trigger ETU beds in an admin unit", "%i", (void *) &(P.InitCasesToETUBeds), 1, 1, 0)) P.InitCasesToETUBeds = 1e6;
					if (!GetInputParameter2(dat, dat2, "Initial time to install ETU beds", "%lf", (void *) &(P.InitDelayToETUBeds), 1, 1, 0)) P.InitDelayToETUBeds = 0;
					if (!GetInputParameter2(dat, dat2, "Initial number of ETU beds", "%i", (void *) &(P.InitNumETUBeds), 1, 1, 0)) P.InitNumETUBeds = 0;
					if (!GetInputParameter2(dat, dat2, "Hospital capacity reached before increasing ETU bed numbers", "%lf", (void *) &(P.CapacityToMoreETUBeds), 1, 1, 0)) P.CapacityToMoreETUBeds = 0;
					if (!GetInputParameter2(dat, dat2, "Subsequent time to install ETU beds", "%lf", (void *) &(P.SubDelayToETUBeds), 1, 1, 0)) P.SubDelayToETUBeds = 0;
					if (!GetInputParameter2(dat, dat2, "Subsequent number of additional ETU beds", "%i", (void *) &(P.SubNumETUBeds), 1, 1, 0)) P.SubNumETUBeds = 0;
				}
				if (!GetInputParameter2(dat, dat2, "Output ETU capacity", "%i", (void*)&(P.DoOutputETUCapacity), 1, 1, 0)) P.DoOutputETUCapacity = 0;
			}
		}
		else P.DoETUByAdUnit=0;

		//This section is going to be specific to hospitals rather than ETUs
		if (!GetInputParameter2(dat, dat2, "Include hospitals as places", "%i", (void*)&(P.IncludeHospitalPlaceType), 1, 1, 0)) P.IncludeHospitalPlaceType = 0;
		if (P.IncludeHospitalPlaceType)
		{
			if (!GetInputParameter2(dat, dat2, "Output keyworker file", "%i", (void*)&(P.DoKeyworkerOutput), 1, 1, 0)) P.DoKeyworkerOutput = 0;
			GetInputParameter(dat, dat2, "Hospital place type number", "%i", (void*)&(P.HospPlaceTypeNum), 1, 1, 0);
			GetInputParameter(dat, dat2, "Max number of beds for Ebola cases in hospital", "%i", (void*)&(P.HospCaseCapacity), 1, 1, 0);
			GetInputParameter(dat, dat2, "Number of HCWs per 1000 population", "%lf", (void*)&(P.HCWPerThousand), 1, 1, 0);
			GetInputParameter(dat, dat2, "Proportion of HCWs and FLWs vaccinated before outbreak", "%lf", (void*) &(P.PropHCWFLWVacc),1, 1, 0);
			if (!GetInputParameter2(dat, dat2, "Relative susceptibility of HCWs and FLWs vaccinated before outbreak", "%lf", (void*)&(P.HCWVaccSuscDrop), 1, 1, 0)) P.HCWVaccSuscDrop = -1;
			if (!GetInputParameter2(dat, dat2, "Revaccinate HCW and FLWs after outbreak alert", "%i", (void*)&(P.RevaccHCWs), 1, 1, 0)) P.RevaccHCWs = 0;
			GetInputParameter(dat, dat2, "Number of days between HCW and FLW vaccination and outbreak", "%i", (void*)&(P.DayHCWFLWVacc), 1, 1, 0);
			GetInputParameter(dat, dat2, "Min age for HCWs and FLWs", "%i", (void*)&(P.MinAgeHCWFLW), 1, 1, 0);
			GetInputParameter(dat, dat2, "Max age for HCWs and FLWs", "%i", (void*)&(P.MaxAgeHCWFLW), 1, 1, 0);
			if (!GetInputParameter2(dat, dat2, "Include FLWs", "%i", (void*)&(P.IncludeFLWs), 1, 1, 0)) P.IncludeFLWs = 0;
			if (P.IncludeFLWs)
			{
				GetInputParameter(dat, dat2, "Number of FLWs per 1000 population", "%lf", (void*) &(P.FLWPerThousand), 1, 1, 0);
				GetInputParameter(dat, dat2, "Relative suspectibility of FLWs", "%lf", (void*) &(P.RelSuscFLW), 1, 1, 0);
			}
			if(!GetInputParameter2(dat, dat2, "Relative susceptibility of HCWs and FLWs due to PPE post outbreak detection", "%lf", (void*&)P.RelSuscPPE, 1, 1, 0)) P.RelSuscPPE=1;
		}
		
		if(!GetInputParameter2(dat,dat2,"Relative infectiousness of an ETU case","%lf",(void *) &(P.RelativeInfectiousnessETU),1,1,0)) P.RelativeInfectiousnessETU=1;
	}
	
	if (!GetInputParameter2(dat, dat2, "Do interrupt interventions", "%i", (void*) & (P.DoInterruptIntervention), 1, 1, 0)) P.DoInterruptIntervention = 0; //Interruption of interventions parameters: ggilani 28/10/14
	if (P.DoInterruptIntervention)
	{
		if (!GetInputParameter2(dat, dat2, "Number of days of interruptions", "%i", (void*) & (P.NDaysInterrupt), 1, 1, 0)) P.NDaysInterrupt = 0;
		if (P.NDaysInterrupt > 0)
		{
			if (P.NDaysInterrupt >= MAX_CHANGE_POINTS) ERR_CRITICAL("MAX_CHANGE_POINTS too small\n");
			GetInputParameter(dat, dat2, "Days of interruptions", "%i", (void*) &(P.DaysInterruptIntervention), P.NDaysInterrupt, 1, 0);
		}
	}

	if(!GetInputParameter2(dat,dat2,"Do prevalence dependent transmission","%i",(void *) &(P.DoPrevDependTrans),1,1,0)) P.DoPrevDependTrans=0;
	//added prevalent dependent transmission - ggilani, 22/12/14
	if(P.DoPrevDependTrans)
	{
		// determine type of prevalence dependent transmission - either by total number of cases so far (which will never go down i.e. once we're over the threshold, people's behaviour is permanently changed) or by current number of infections (i.e. it can go down again)
		if(!GetInputParameter2(dat,dat2,"Do prevalence dependent transmission by total cases","%i",(void *) &(P.DoPrevDepTransTotalCases),1,1,0)) P.DoPrevDepTransTotalCases=0;
		if(!GetInputParameter2(dat,dat2,"Do prevalence dependent transmission by current infections","%i",(void *) &(P.DoPrevDepTransCurrInf),1,1,0)) P.DoPrevDepTransCurrInf=0;
		// we can only do one version though
		if((P.DoPrevDepTransTotalCases)&&(P.DoPrevDepTransCurrInf))
		{
			ERR_CRITICAL("Only one type of prevalent dependent transmission can be specified!\n");
		}
		if(!GetInputParameter2(dat,dat2,"Threshold for prevalence dependent transmission","%lf",(void *) &(P.PrevDependTransThresh),1,1,0)) P.PrevDependTransThresh=1.0;
		if(!GetInputParameter2(dat,dat2,"Threshold for prevalence dependent transmission","%lf",(void *) &(P.PrevDependRelativeSusc),1,1,0)) P.PrevDependRelativeSusc=1.0;

	}
	if(!GetInputParameter2(dat,dat2,"Include latent period","%i",(void *) &(P.DoLatent),1,1,0)) P.DoLatent=0;
	if(P.DoLatent)
		{
		GetInputParameter(dat,dat2,"Latent period","%lf",(void *) &(P.LatentPeriod),1,1,0);
		if(!GetInputParameter2(dat,dat2,"Latent period inverse CDF","%lf",(void *) P.latent_icdf,CDF_RES+1,1,0))
			{
			P.latent_icdf[CDF_RES]=1e10;
			for(i=0;i<CDF_RES;i++)
				P.latent_icdf[i]=-log(1-((double)i)/CDF_RES);
			}
		for(i=0;i<=CDF_RES;i++)
			P.latent_icdf[i]=exp(-P.latent_icdf[i]);
		}

	if(!GetInputParameter2(dat,dat2,"Include symptoms","%i",(void *) &(P.DoSymptoms),1,1,0)) P.DoSymptoms=0;
	if(!P.DoSymptoms)
		{
		for(i=0;i<NUM_AGE_GROUPS;i++)
			P.ProportionSymptomatic[i]=0;
		P.FalsePositiveRate=0;
		P.SymptInfectiousness=1.0;
		P.LatentToSymptDelay=0;
		}
	else
		{
		if(P.DoAge)
			GetInputParameter(dat,dat2,"Proportion symptomatic by age group","%lf",(void *) P.ProportionSymptomatic,NUM_AGE_GROUPS,1,0);
		else
			{
			GetInputParameter(dat,dat2,"Proportion symptomatic","%lf",(void *) P.ProportionSymptomatic,1,1,0);
			for(i=1;i<NUM_AGE_GROUPS;i++)
				P.ProportionSymptomatic[i]=P.ProportionSymptomatic[0];
			}
		GetInputParameter(dat,dat2,"Delay from end of latent period to start of symptoms","%lf",(void *) &(P.LatentToSymptDelay),1,1,0);
		GetInputParameter(dat,dat2,"Relative rate of random contacts if symptomatic","%lf",(void *) &(P.SymptSpatialContactRate),1,1,0);
		GetInputParameter(dat,dat2,"Symptomatic infectiousness relative to asymptomatic","%lf",(void *) &(P.SymptInfectiousness),1,1,0);
		if(!GetInputParameter2(dat,dat2,"Model symptomatic withdrawal to home as true absenteeism","%i",(void *) &P.DoRealSymptWithdrawal,1,1,0)) P.DoRealSymptWithdrawal=0;
		if(P.DoPlaces)
			{
			GetInputParameter(dat,dat2,"Relative level of place attendance if symptomatic","%lf",(void *) P.SymptPlaceTypeContactRate,P.PlaceTypeNum,1,0);
			if(P.DoRealSymptWithdrawal)
				{
				for(j=0;j<NUM_PLACE_TYPES;j++)
					{
					P.SymptPlaceTypeWithdrawalProp[j]=1.0-P.SymptPlaceTypeContactRate[j];
					P.SymptPlaceTypeContactRate[j]=1.0;
					}
				}
			else
				for(j=0;j<NUM_PLACE_TYPES;j++) P.SymptPlaceTypeWithdrawalProp[j]=0.0;
			}
		if(!GetInputParameter2(dat,dat2,"Maximum age of child at home for whom one adult also stays at home","%i",(void *) &P.CaseAbsentChildAgeCutoff,1,1,0)) P.CaseAbsentChildAgeCutoff=0;
		if(!GetInputParameter2(dat,dat2,"Proportion of children at home for whom one adult also stays at home","%lf",(void *) &P.CaseAbsentChildPropAdultCarers,1,1,0)) P.CaseAbsentChildPropAdultCarers=0;
#ifdef ABSENTEEISM_PLACE_CLOSURE
		P.CaseAbsenteeismDelay=0;  // Set to zero for tracking absenteeism
#else
		if(!GetInputParameter2(dat,dat2,"Delay in starting place absenteeism for cases who withdraw","%lf",(void *) &P.CaseAbsenteeismDelay,1,1,0)) P.CaseAbsenteeismDelay=0;
#endif
		if(!GetInputParameter2(dat,dat2,"Duration of place absenteeism for cases who withdraw","%lf",(void *) &P.CaseAbsenteeismDuration,1,1,0)) P.CaseAbsenteeismDuration=7;

		if(!GetInputParameter2(dat,dat2,"False positive rate","%lf",(void *) &(P.FalsePositiveRate),1,1,0)) P.FalsePositiveRate=0.0;
		if(!GetInputParameter2(dat,dat2,"False positive per capita incidence","%lf",(void *) &(P.FalsePositivePerCapitaIncidence),1,1,0)) P.FalsePositivePerCapitaIncidence=0.0;
		if(!GetInputParameter2(dat,dat2,"False positive relative incidence by age","%lf",(void *) P.FalsePositiveAgeRate,NUM_AGE_GROUPS,1,0))
			for(j=0;j<NUM_AGE_GROUPS;j++) P.FalsePositiveAgeRate[j]=1.0;
		}
	if(!GetInputParameter2(dat,dat2,"Bounding box for bitmap","%lf",(void *) &(P.BoundingBox[0]),4,1,0))
		{
		P.BoundingBox[0]=P.BoundingBox[1]=0.0;
		P.BoundingBox[2]=P.BoundingBox[3]=1.0;
		}
	if(!GetInputParameter2(dat,dat2,"Spatial domain for simulation","%lf",(void *) &(P.SpatialBoundingBox[0]),4,1,0))
		{
		P.SpatialBoundingBox[0]=P.SpatialBoundingBox[1]=0.0;
		P.SpatialBoundingBox[2]=P.SpatialBoundingBox[3]=1.0;
		}
	if(!GetInputParameter2(dat,dat2,"Grid size","%lf",(void *) &(P.cwidth),1,1,0)) P.cwidth=1.0/120.0;
	if(!GetInputParameter2(dat,dat2,"Use long/lat coord system","%i",(void *) &(P.DoUTM_coords),1,1,0)) P.DoUTM_coords=0;
	if(!GetInputParameter2(dat,dat2,"Bitmap scale","%lf",(void *) &(P.BitmapScale),1,1,0)) P.BitmapScale=1.0;
	if(!GetInputParameter2(dat,dat2,"Bitmap y:x aspect scaling","%lf",(void *) &(P.BitmapAspectScale),1,1,0)) P.BitmapAspectScale=1.0;
	if(!GetInputParameter2(dat,dat2,"Bitmap movie frame interval","%i",(void *) &(P.BitmapMovieFrame),1,1,0)) P.BitmapMovieFrame=250;
	if(!GetInputParameter2(dat,dat2,"Output bitmap","%i",(void *) &(P.OutputBitmap),1,1,0)) P.OutputBitmap=0;
	if(!GetInputParameter2(dat,dat2,"Output bitmap detected","%i",(void *) &(P.OutputBitmapDetected),1,1,0)) P.OutputBitmapDetected=0;
	if(!GetInputParameter2(dat,dat2,"Output immunity on bitmap","%i",(void *) &(P.DoImmuneBitmap),1,1,0)) P.DoImmuneBitmap=0;
	if(!GetInputParameter2(dat,dat2,"Output infection tree","%i",(void *) &(P.DoInfectionTree),1,1,0)) P.DoInfectionTree=0;
	if(!GetInputParameter2(dat,dat2,"Do one generation","%i",(void *) &(P.DoOneGen),1,1,0)) P.DoOneGen=0;
	if(!GetInputParameter2(dat,dat2,"Output every realisation","%i",(void *) &(P.OutputAll),1,1,0)) P.OutputAll=1;
	if(!GetInputParameter2(dat,dat2,"Maximum number to sample for correlations","%i",(void *) &(P.MaxCorrSample),1,1,0)) P.MaxCorrSample=1000000000;
	if(!GetInputParameter2(dat,dat2,"Assume SI model","%i",(void *) &(P.DoSI),1,1,0)) P.DoSI=0;
	if(!GetInputParameter2(dat,dat2,"Assume periodic boundary conditions","%i",(void *) &(P.DoPeriodicBoundaries),1,1,0)) P.DoPeriodicBoundaries=0;
	if(!GetInputParameter2(dat,dat2,"Only output non-extinct realisations","%i",(void *) &(P.OutputNonExtinct),1,1,0)) P.OutputNonExtinct=0;

	if (!GetInputParameter2(dat, dat2, "Output control file", "%i", (void*)&(P.DoControlOutput), 1, 1, 0)) P.DoControlOutput = 0;
	if(!GetInputParameter2(dat,dat2,"Use cases per thousand threshold for area controls","%i",(void *) &(P.DoPerCapitaTriggers),1,1,0)) P.DoPerCapitaTriggers=0;
	if(!GetInputParameter2(dat,dat2,"Use global triggers for interventions","%i",(void *) &(P.DoGlobalTriggers),1,1,0)) P.DoGlobalTriggers=0;
	if(!GetInputParameter2(dat,dat2,"Divisor for per-capita area threshold (default 1000)","%i",(void *) &(P.IncThreshPop),1,1,0)) P.IncThreshPop=1000;
	if(!GetInputParameter2(dat,dat2,"Divisor for per-capita global threshold (default 1000)","%i",(void *) &(P.GlobalIncThreshPop),1,1,0)) P.GlobalIncThreshPop=1000;

	
	if(!GetInputParameter2(dat,dat2,"Number of sampling intervals over which cumulative incidence measured for global trigger","%i",(void *) &(P.TriggersSamplingInterval),1,1,0)) P.TriggersSamplingInterval=10000000;
	//if (!GetInputParameter2(dat, dat2, "Number of undetected infections before first case detected", "%i", (void*)&(P.NumUndetectedInfPreOutbreakAlert), 1, 1, 0)) P.NumUndetectedInfPreOutbreakAlert = 0;
	if(!GetInputParameter2(dat,dat2,"Proportion of cases detected after surveillance alert","%lf",(void *) &(P.PostAlertControlPropCasesId),1,1,0)) P.PostAlertControlPropCasesId=1;
	//if(!GetInputParameter2(dat,dat2,"Proportion of cases detected before surveillance alert","%lf",(void *) &(P.PreAlertControlPropCasesId),1,1,0)) P.PreAlertControlPropCasesId=1;
	//if(P.PreControlClusterIdCaseThreshold==0)
	//	{
	if(!GetInputParameter2(dat,dat2,"Number of detected cases before outbreak alert triggered","%i",(void *) &(P.PreControlClusterIdCaseThreshold),1,1,0)) P.PreControlClusterIdCaseThreshold=0;
	//	}
	if(!GetInputParameter2(dat,dat2,"Only use confirmed cases to trigger alert","%i",(void *) &(P.DoEarlyCaseDiagnosis),1,1,0)) P.DoEarlyCaseDiagnosis=0;
	if(!GetInputParameter2(dat,dat2,"Target country","%i",(void *) &(P.TargetCountry),1,1,0)) P.TargetCountry=1;
	fprintf(stderr,"Target country=%i\n",P.TargetCountry);
	//added extra target countries to help keep track of country wise parameters: ggilani
	if(!GetInputParameter2(dat,dat2,"Target country 2","%i",(void *) &(P.TargetCountry2),1,1,0)) P.TargetCountry2=1;
	if(!GetInputParameter2(dat,dat2,"Target country 3","%i",(void *) &(P.TargetCountry3),1,1,0)) P.TargetCountry3=1;
	if(!GetInputParameter2(dat,dat2,"Restrict treatment to target country","%i",(void *) &(P.RestrictTreatToTarget),1,1,0)) P.RestrictTreatToTarget=0;
	if(!GetInputParameter2(dat,dat2,"Only treat mixing groups within places","%i",(void *) &(P.DoPlaceGroupTreat),1,1,0)) P.DoPlaceGroupTreat=0;

	if(!GetInputParameter2(dat,dat2,"Treatment trigger incidence per cell","%lf",(void *) &(P.TreatCellIncThresh),1,1,0)) P.TreatCellIncThresh=1000000000;
	if(!GetInputParameter2(dat,dat2,"Relative susceptibility of treated individual","%lf",(void *) &(P.TreatSuscDrop),1,1,0)) P.TreatSuscDrop=1;
	if(!GetInputParameter2(dat,dat2,"Relative infectiousness of treated individual","%lf",(void *) &(P.TreatInfDrop),1,1,0)) P.TreatInfDrop=1;
	if(!GetInputParameter2(dat,dat2,"Proportion of symptomatic cases resulting in death prevented by treatment","%lf",(void *) &(P.TreatDeathDrop),1,1,0)) P.TreatDeathDrop=0;
	if(!GetInputParameter2(dat,dat2,"Proportion of symptomatic cases prevented by treatment","%lf",(void *) &(P.TreatSympDrop),1,1,0)) P.TreatSympDrop=0;
	if(!GetInputParameter2(dat,dat2,"Delay to treat cell","%lf",(void *) &(P.TreatDelayMean),1,1,0)) P.TreatDelayMean=0;
	if(!GetInputParameter2(dat,dat2,"Duration of course of treatment","%lf",(void *) &(P.TreatCaseCourseLength),1,1,0)) P.TreatCaseCourseLength=5;
	if(!GetInputParameter2(dat,dat2,"Duration of course of prophylaxis","%lf",(void *) &(P.TreatProphCourseLength),1,1,0)) P.TreatProphCourseLength=10;
	if(!GetInputParameter2(dat,dat2,"Proportion of detected cases treated","%lf",(void *) &(P.TreatPropCases),1,1,0)) P.TreatPropCases=1;
	if(!GetInputParameter2(dat,dat2,"Proportion of detected cases with private stockpile treated","%lf",(void *) &(P.PrivateTreatPropCases),1,1,0)) P.PrivateTreatPropCases=0;
	if(P.DoHouseholds)
		{
		if(!GetInputParameter2(dat,dat2,"Delay to treat households with private stockpiles","%lf",(void *) &(P.PrivateTreatDelayMean),1,1,0)) P.PrivateTreatDelayMean=0;
		if(!GetInputParameter2(dat,dat2,"Proportion of households of cases treated","%lf",(void *) &(P.TreatPropCaseHouseholds),1,1,0)) P.TreatPropCaseHouseholds=0;
		if(!GetInputParameter2(dat,dat2,"Proportion of households with private stockpiles treated","%lf",(void *) &(P.PrivateTreatPropCaseHouseholds),1,1,0)) P.PrivateTreatPropCaseHouseholds=0;
		if(!GetInputParameter2(dat,dat2,"Duration of household prophylaxis policy","%lf",(void *) &(P.TreatHouseholdsDuration),1,1,0)) P.TreatHouseholdsDuration=USHRT_MAX/P.TimeStepsPerDay;
		}
	if(!GetInputParameter2(dat,dat2,"Proportion treated","%lf",(void *) &(P.TreatPropRadial),1,1,0)) P.TreatPropRadial=1.0;
	if(!GetInputParameter2(dat,dat2,"Proportion treated in radial prophylaxis","%lf",(void *) &(P.TreatPropRadial),1,1,0)) P.TreatPropRadial=1.0;
	if(!GetInputParameter2(dat,dat2,"Treatment radius","%lf",(void *) &(P.TreatRadius),1,1,0)) P.TreatRadius=0;
	if(!GetInputParameter2(dat,dat2,"Duration of place/geographic prophylaxis policy","%lf",(void *) &(P.TreatPlaceGeogDuration),1,1,0)) P.TreatPlaceGeogDuration=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat,dat2,"Treatment start time","%lf",(void *) &(P.TreatTimeStartBase),1,1,0)) P.TreatTimeStartBase=USHRT_MAX/P.TimeStepsPerDay;
	if(P.DoPlaces)
		{
		if(!GetInputParameter2(dat,dat2,"Proportion of places treated after case detected","%lf",(void *) P.TreatPlaceProbCaseId,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++) P.TreatPlaceProbCaseId[i]=0;
		if(!GetInputParameter2(dat,dat2,"Proportion of people treated in targeted places","%lf",(void *) P.TreatPlaceTotalProp,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++) P.TreatPlaceTotalProp[i]=0;
		}
	if(!GetInputParameter2(dat,dat2,"Maximum number of doses available","%lf",(void *) &(P.TreatMaxCoursesBase),1,1,0)) P.TreatMaxCoursesBase=1e20;
	if(!GetInputParameter2(dat,dat2,"Start time of additional treatment production","%lf",(void *) &(P.TreatNewCoursesStartTime),1,1,0)) P.TreatNewCoursesStartTime=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat,dat2,"Rate of additional treatment production (courses per day)","%lf",(void *) &(P.TreatNewCoursesRate),1,1,0)) P.TreatNewCoursesRate=0;
	if(!GetInputParameter2(dat,dat2,"Maximum number of people targetted with radial prophylaxis per case","%i",(void *) &(P.TreatMaxCoursesPerCase),1,1,0)) P.TreatMaxCoursesPerCase=1000000000;
	if(!GetInputParameter2(dat,dat2,"Proportion of households buying private stockpile","%lf",(void *) &(P.PropPrivateStockpile),1,1,0)) P.PropPrivateStockpile=0;
	if(!GetInputParameter2(dat,dat2,"Income ordering of private stockpile (-1=low to high, 0=random, 1=high to low)","%i",(void *) &(P.PrivateStockpileOrderByIncome),1,1,0)) P.PrivateStockpileOrderByIncome=0;

	if(!GetInputParameter2(dat,dat2,"Number of resistance levels","%i",(void *) &(P.EvolResistNumTypes),1,1,0)) P.EvolResistNumTypes=0;
	P.EvolResistNumTypes++;
	if(P.EvolResistNumTypes>MAX_NUM_RESIST_TYPES) ERR_CRITICAL("Too many resistance levels\n");
	if(!GetInputParameter2(dat,dat2,"Mutation probability per treated infection of step increase in resistance","%lf",(void *) &(P.EvolResistTreatMutationRate),1,1,0)) P.EvolResistTreatMutationRate=0;
	if(!GetInputParameter2(dat,dat2,"Mutation probability per prophylaxed infection of step increase in resistance","%lf",(void *) &(P.EvolResistProphMutationRate),1,1,0)) P.EvolResistProphMutationRate=0;
	P.EvolResistRelInf[0]=1.0;
	P.EvolResistRelTreatInfDrop[0]=0.0;
	P.EvolResistRelProphSusc[0]=0.0;
	P.EvolResistRelTreatSympDrop[0]=0.0;
	P.EvolResistRelTreatDeathDrop[0]=0.0;
	P.EvolResistSeedProp[0]=1.0;
	if(!GetInputParameter2(dat,dat2,"Proportion of seed infections with different resistance levels","%lf",(void *) (P.EvolResistSeedProp+1),P.EvolResistNumTypes-1,1,0))
		{
		for(i=1;i<P.EvolResistNumTypes;i++) P.EvolResistSeedProp[i]=0.0;
		}
	else
		{
		for(i=1;i<P.EvolResistNumTypes;i++) P.EvolResistSeedProp[0]-=P.EvolResistSeedProp[i];
		}
	if(!GetInputParameter2(dat,dat2,"Relative infectiousness of resistance levels","%lf",(void *) (P.EvolResistRelInf+1),P.EvolResistNumTypes-1,1,0)) 
		{
		for(i=1;i<P.EvolResistNumTypes;i++) P.EvolResistRelInf[i]=1.0;
		}
	if(!GetInputParameter2(dat,dat2,"Reduction in drop in infectiousness from treatment caused by resistance levels","%lf",(void *) (P.EvolResistRelTreatInfDrop+1),P.EvolResistNumTypes-1,1,0)) 
		{
		for(i=1;i<P.EvolResistNumTypes;i++) P.EvolResistRelTreatInfDrop[i]=0.0;
		}
	if(!GetInputParameter2(dat,dat2,"Reduction in drop in susceptiblity from prophylaxis caused by resistance levels","%lf",(void *) (P.EvolResistRelProphSusc+1),P.EvolResistNumTypes-1,1,0)) 
		{
		for(i=1;i<P.EvolResistNumTypes;i++) P.EvolResistRelProphSusc[i]=0.0;
		}
	if(!GetInputParameter2(dat,dat2,"Reduction in drop in symptoms from treatment caused by resistance levels","%lf",(void *) (P.EvolResistRelTreatSympDrop+1),P.EvolResistNumTypes-1,1,0)) 
		{
		for(i=1;i<P.EvolResistNumTypes;i++) P.EvolResistRelTreatSympDrop[i]=0.0;
		}
	if(!GetInputParameter2(dat,dat2,"Reduction in drop in deaths from treatment caused by resistance levels","%lf",(void *) (P.EvolResistRelTreatDeathDrop+1),P.EvolResistNumTypes-1,1,0)) 
		{
		for(i=1;i<P.EvolResistNumTypes;i++) P.EvolResistRelTreatDeathDrop[i]=0.0;
		}
	for(i=0;i<P.EvolResistNumTypes;i++)
		{
		P.EvolResistRelTreatInfDrop[i]=1-(1-P.TreatInfDrop)*(1-P.EvolResistRelTreatInfDrop[i]);
		P.EvolResistRelProphSusc[i]=1-(1-P.TreatSuscDrop)*(1-P.EvolResistRelProphSusc[i]);
		P.EvolResistRelTreatSympDrop[i]=1-P.TreatSympDrop*(1-P.EvolResistRelTreatSympDrop[i]);
		P.EvolResistRelTreatDeathDrop[i]=1-P.TreatDeathDrop*(1-P.EvolResistRelTreatDeathDrop[i]);
		}
	if(P.DoAdUnits)
		{
		if(!GetInputParameter2(dat,dat2,"Treat administrative units rather than rings","%i",(void *) &(P.TreatByAdminUnit),1,1,0)) P.TreatByAdminUnit=0;
		if(!GetInputParameter2(dat,dat2,"Administrative unit divisor for treatment","%i",(void *) &(P.TreatAdminUnitDivisor),1,1,0)) P.TreatAdminUnitDivisor=1;
		if((P.TreatAdminUnitDivisor==0)||(P.TreatByAdminUnit==0)) {P.TreatByAdminUnit=0;P.TreatAdminUnitDivisor=1;}
		}
	else
		{P.TreatAdminUnitDivisor=1;P.TreatByAdminUnit=0;}

	if(!GetInputParameter2(dat,dat2,"Vaccination trigger incidence per cell","%lf",(void *) &(P.VaccCellIncThresh),1,1,0)) P.VaccCellIncThresh=1000000000;
	if(P.VaccCaseScale>1&&P.VaccCellIncThresh==1)
	{
		P.VaccCellIncThresh=P.VaccCellIncThresh*P.VaccCaseScale;
	}
	if(!GetInputParameter2(dat,dat2,"Relative susceptibility of vaccinated individual","%lf",(void *) &(P.VaccSuscDrop),1,1,0)) P.VaccSuscDrop=1;
	if (P.HCWVaccSuscDrop <0 )
	{
		P.HCWVaccSuscDrop = P.VaccSuscDrop;
	}
	if(!GetInputParameter2(dat,dat2,"Relative susceptibility of individual vaccinated after switch time","%lf",(void *) &(P.VaccSuscDrop2),1,1,0)) P.VaccSuscDrop2=1;
	if(!GetInputParameter2(dat,dat2,"Switch time at which vaccine efficacy increases","%lf",(void *) &(P.VaccTimeEfficacySwitch),1,1,0)) P.VaccTimeEfficacySwitch=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat,dat2,"Decay rate of vaccine efficacy (per year)","%lf",(void *) &(P.VaccEfficacyDecay),1,1,0)) P.VaccEfficacyDecay=0;
	P.VaccEfficacyDecay/=DAYS_PER_YEAR;
	if(!GetInputParameter2(dat,dat2,"Relative infectiousness of vaccinated individual","%lf",(void *) &(P.VaccInfDrop),1,1,0)) P.VaccInfDrop=1;
	if(!GetInputParameter2(dat,dat2,"Proportion of symptomatic cases resulting in death prevented by vaccination","%lf",(void *) &(P.VaccMortDrop),1,1,0)) P.VaccMortDrop=0;
	if(!GetInputParameter2(dat,dat2,"Proportion of symptomatic cases prevented by vaccination","%lf",(void *) &(P.VaccSympDrop),1,1,0)) P.VaccSympDrop=0;
	if(!GetInputParameter2(dat,dat2,"Delay to vaccinate","%lf",(void *) &(P.VaccDelayMean),1,1,0)) P.VaccDelayMean=0;
	if((P.VaccDelayMean==1)&(P.VaccDelayScale!=1))
	{
		P.VaccDelayMean=P.VaccDelayMean*P.VaccDelayScale; //added this to also us to scale delay to vaccination time from command line - ggilani 28/02/17
	}
	if(!GetInputParameter2(dat,dat2,"Delay from vaccination to full protection","%lf",(void *) &(P.VaccTimeToEfficacy),1,1,0)) P.VaccTimeToEfficacy=0;
	if((P.VaccTimeToEfficacy==1)&(P.VaccEffTimeScale!=1))
	{
		P.VaccTimeToEfficacy=P.VaccTimeToEfficacy*P.VaccEffTimeScale; //added this to also us to scale time to vaccination protection from command line - ggilani 28/02/17
	}
	if(!GetInputParameter2(dat,dat2,"Years between rounds of vaccination","%lf",(void *) &(P.VaccCampaignInterval),1,1,0)) P.VaccCampaignInterval=1e10;
	if(!GetInputParameter2(dat, dat2, "Base vaccine doses per day", "%i", (void*)&(P.BaseVaccDosePerDay), 1, 1, 0)) P.BaseVaccDosePerDay = -1;
	if(!GetInputParameter2(dat,dat2,"Max vaccine doses per day","%i",(void *) &(P.MaxVaccDosePerDay),1,1,0)) P.MaxVaccDosePerDay=-1;
	if (!GetInputParameter2(dat, dat2, "Base geo vaccine doses per day", "%i", (void*)&(P.BaseVaccGeoDosePerDay), 1, 1, 0)) P.BaseVaccGeoDosePerDay = -1;
	if (!GetInputParameter2(dat, dat2, "Max geo vaccine doses per day", "%i", (void*)&(P.MaxVaccGeoDosePerDay), 1, 1, 0)) P.MaxVaccGeoDosePerDay = -1;
	if (!GetInputParameter2(dat, dat2, "Reset vaccination queue each day", "%i", (void*)&(P.ResetVaccQueue), 1, 1, 0)) P.ResetVaccQueue = 0;
	P.VaccCampaignInterval*=DAYS_PER_YEAR;
	if(!GetInputParameter2(dat,dat2,"Maximum number of rounds of vaccination","%i",(void *) &(P.VaccMaxRounds),1,1,0)) P.VaccMaxRounds=1;
	if(P.DoHouseholds)
		{
		if(!GetInputParameter2(dat,dat2,"Proportion of households of cases vaccinated","%lf",(void *) &(P.VaccPropCaseHouseholds),1,1,0)) P.VaccPropCaseHouseholds=0;
		if(!GetInputParameter2(dat,dat2,"Duration of household vaccination policy","%lf",(void *) &(P.VaccHouseholdsDuration),1,1,0)) P.VaccHouseholdsDuration=USHRT_MAX/P.TimeStepsPerDay;
		}
	if(P.DoHouseholds&&P.DoPlaces)
	{
		if(!GetInputParameter2(dat,dat2,"Do ring vaccination","%i",(void *) &(P.DoRingVaccination),1,1,0)) P.DoRingVaccination=0;
		if(P.DoRingVaccination)
		{
			//if(!GetInputParameter2(dat,dat2,"Number of rings to vaccinate","%i",(void *) &(P.NVaccRings),1,1,0)) P.NVaccRings=1;
			//if(P.VaccRingScale>1&&P.NVaccRings==1)
			//{
			//	P.NVaccRings=P.NVaccRings*P.VaccRingScale;
			//}
			
			if (!GetInputParameter2(dat, dat2, "Ring vaccination trigger incidence per cell", "%lf", (void*)&(P.RingVaccCellIncThresh), 1, 1, 0)) P.RingVaccCellIncThresh = 1000000000;
			if (!GetInputParameter2(dat, dat2, "Ring vaccination start time", "%lf", (void*)&(P.RingVaccTimeStartBase), 1, 1, 0)) P.RingVaccTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;

			if (!GetInputParameter2(dat, dat2, "Number of proportion of ring to vaccinate", "%i", (void *) &(P.NPropRingVacc), 1, 1, 0)) P.NPropRingVacc = 0;
			if (P.NPropRingVacc > 0)
			{
				if (P.NPropRingVacc >= MAX_CHANGE_POINTS) ERR_CRITICAL("MAX_CHANGE_POINTS too small\n");
				GetInputParameter(dat, dat2, "Change points proportion of ring to vaccinate", "%lf", (void *)P.ChangePointPropRingVacc, P.NPropRingVacc, 1, 0);
				GetInputParameter(dat, dat2, "List of proportions of ring to vaccinate", "%lf", (void *)P.ListPropRingVacc, P.NPropRingVacc, 1, 0);
			}
			P.CurrIndPropRingVacc = 0;

			//add time to increase to third ring
			if (!GetInputParameter2(dat, dat2, "Time to increase vaccination rings", "%lf", (void *) &(P.TimeToIncVaccRing), 1, 1, 0)) P.TimeToIncVaccRing = 1e10;
			//add time to efficiency for the third ring
			if (!GetInputParameter2(dat, dat2, "Delay from vaccination to full protection third ring", "%lf", (void *) &(P.VaccTimeToEfficacyThirdVaccRing), 1, 1, 0)) P.VaccTimeToEfficacyThirdVaccRing = 0;

			if(!GetInputParameter2(dat,dat2,"Proportion of ring to vaccinate","%lf",(void *) &(P.PropRingVacc),1,1,0)) P.PropRingVacc=1;
			if(P.PropRingVacc==1&&P.VaccPropScale!=1)
			{
				P.PropRingVacc=P.PropRingVacc*P.VaccPropScale;
			}
			if(!GetInputParameter2(dat,dat2,"Minimum vaccination age","%i",(void *) &(P.MinVaccAge),1,1,0)) P.MinVaccAge=0;
		}
	}
	if(!GetInputParameter2(dat,dat2,"Vaccination start time","%lf",(void *) &(P.VaccTimeStartBase),1,1,0)) P.VaccTimeStartBase=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat,dat2,"Proportion of population vaccinated","%lf",(void *) &(P.VaccProp),1,1,0)) P.VaccProp=0;
	if(!GetInputParameter2(dat,dat2,"Time taken to reach max vaccination coverage (in years)","%lf",(void *) &(P.VaccCoverageIncreasePeriod),1,1,0)) P.VaccCoverageIncreasePeriod=0;
	P.VaccCoverageIncreasePeriod*=DAYS_PER_YEAR;
	if (!GetInputParameter2(dat, dat2, "Do geographic vaccination", "%i", (void*)&(P.DoGeoVaccination), 1, 1, 0)) P.DoGeoVaccination = 0;
	if (P.DoGeoVaccination)
	{
		if (!GetInputParameter2(dat, dat2, "Only do geographic vaccination when ring not established", "%i", (void*)&(P.OnlyDoGeoVaccWhenNoRing), 1, 1, 0)) P.OnlyDoGeoVaccWhenNoRing = 0;
	}
	if (!GetInputParameter2(dat, dat2, "Output vaccine distributions", "%i", (void*)&(P.DoVaccOutput), 1, 1, 0)) P.DoVaccOutput = 0;
	if (!GetInputParameter2(dat, dat2, "Do clustered vaccination acceptance by household", "%i", (void*)&(P.DoClusterVaccAccept), 1, 1, 0)) P.DoClusterVaccAccept = 0;
	//if (!GetInputParameter2(dat, dat2, "Geo vaccination trigger incidence per cell", "%lf", (void*)&(P.GeoVaccCellIncThresh), 1, 1, 0)) P.GeoVaccCellIncThresh = 1000000000;
	//if (!GetInputParameter2(dat, dat2, "Geo vaccination start time", "%lf", (void*)&(P.GeoVaccTimeStartBase), 1, 1, 0)) P.GeoVaccTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Vaccination radius", "%lf", (void*)&(P.VaccRadius), 1, 1, 0)) P.VaccRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Vaccination radius", "%i", (void*)&(P.VaccRadius), 1, 1, 0)) P.VaccRadius = 0;
	//if (!GetInputParameter2(dat, dat2, "Scale geographic vaccination radius by population density", "%i", (void*)&(P.ScaleGeoVaccRadius), 1, 1, 0)) P.ScaleGeoVaccRing = 0;
	if (!GetInputParameter2(dat, dat2, "Limit vaccine doses per case per cell", "%i", (void*)&(P.LimitGeoVaccDosesPerCase), 1, 1, 0)) P.LimitGeoVaccDosesPerCase = 0;
	if (P.LimitGeoVaccDosesPerCase)
	{
		if (!GetInputParameter2(dat, dat2, "Vaccine doses per case per cell", "%i", (void*)&(P.VaccDosesPerCasePerCell), 1, 1, 0)) P.VaccDosesPerCasePerCell = 100;
		if (!GetInputParameter2(dat, dat2, "Stop vaccination in cells once vaccine acceptance threshold reached", "%i", (void*)&(P.StopVaccinationPostThreshold), 1, 1, 0)) P.StopVaccinationPostThreshold = 0;
	}
	if (!P.LimitGeoVaccDosesPerCase)
	{
		if (!GetInputParameter2(dat, dat2, "Population of high density cell", "%i", (void*)&(P.PopHighDensityCell), 1, 1, 0)) P.PopHighDensityCell = 500;
		if (!GetInputParameter2(dat, dat2, "Vaccination radius of high density cell", "%lf", (void*)&(P.VaccRadiusHighDensity), 1, 1, 0)) P.VaccRadiusHighDensity = 500;
	}
	if(!GetInputParameter2(dat,dat2,"Minimum radius from case to vaccinate","%lf",(void *) &(P.VaccMinRadius),1,1,0)) P.VaccMinRadius=0;
	if(!GetInputParameter2(dat,dat2,"Maximum number of vaccine courses available","%lf",(void *) &(P.VaccMaxCoursesBase),1,1,0)) P.VaccMaxCoursesBase=1e20;
	if(!GetInputParameter2(dat,dat2,"Start time of additional vaccine production","%lf",(void *) &(P.VaccNewCoursesStartTimeBase),1,1,0)) P.VaccNewCoursesStartTimeBase=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat,dat2,"End time of additional vaccine production","%lf",(void *) &(P.VaccNewCoursesEndTime),1,1,0)) P.VaccNewCoursesEndTime=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat, dat2, "Do daily vaccine replenishment", "%i", (void*)&(P.DoVaccDailyReplenishment), 1, 1, 0)) P.DoVaccDailyReplenishment = 0;
	if (P.DoVaccDailyReplenishment)
	{
		if (!GetInputParameter2(dat, dat2, "Rate of additional vaccine production (courses per day)", "%lf", (void*)&(P.VaccNewCoursesRate), 1, 1, 0)) P.VaccNewCoursesRate = 0;

	}
	if (!GetInputParameter2(dat, dat2, "Do bulk vaccine replenishment", "%i", (void*)&(P.DoVaccBulkReplenishment), 1, 1, 0)) P.DoVaccBulkReplenishment = 0;
	if (P.DoVaccBulkReplenishment)
	{
		if (!GetInputParameter2(dat, dat2, "Number of vaccine doses per delivery", "%lf", (void*)&(P.VaccNewCoursesBulk), 1, 1, 0)) P.VaccNewCoursesBulk = 0;
		if (!GetInputParameter2(dat, dat2, "Time between vaccine doses deliveries (in days)", "%lf", (void*)&(P.VaccNewCoursesDelay), 1, 1, 0)) P.VaccNewCoursesDelay = 0;
	}
	if(!GetInputParameter2(dat,dat2,"Apply mass rather than reactive vaccination","%i",(void *) &(P.DoMassVacc),1,1,0)) P.DoMassVacc=0;
	if(!GetInputParameter2(dat,dat2,"Priority age range for mass vaccination","%i",(void *) P.VaccPriorityGroupAge,2,1,0)) {P.VaccPriorityGroupAge[0]=1;P.VaccPriorityGroupAge[1]=0;}
	if(P.DoAdUnits)
		{
		if(!GetInputParameter2(dat,dat2,"Vaccinate administrative units rather than rings","%i",(void *) &(P.VaccByAdminUnit),1,1,0)) P.VaccByAdminUnit=0;
		if(!GetInputParameter2(dat,dat2,"Administrative unit divisor for vaccination","%i",(void *) &(P.VaccAdminUnitDivisor),1,1,0)) P.VaccAdminUnitDivisor=1;
		if((P.VaccAdminUnitDivisor==0)||(P.VaccByAdminUnit==0)) P.VaccAdminUnitDivisor=1;
		}
	else
		{P.VaccAdminUnitDivisor=1;P.VaccByAdminUnit=0;}
	

	if(!GetInputParameter2(dat,dat2,"Movement restrictions trigger incidence per cell","%i",(void *) &(P.MoveRestrCellIncThresh),1,1,0)) P.MoveRestrCellIncThresh=1000000000;
	if(!GetInputParameter2(dat,dat2,"Delay to start movement restrictions","%lf",(void *) &(P.MoveDelayMean),1,1,0)) P.MoveDelayMean=0;
	if(!GetInputParameter2(dat,dat2,"Duration of movement restrictions","%lf",(void *) &(P.MoveRestrDuration),1,1,0)) P.MoveRestrDuration=7;
	if(!GetInputParameter2(dat,dat2,"Residual movements after restrictions","%lf",(void *) &(P.MoveRestrEffect),1,1,0)) P.MoveRestrEffect=0;
	if(!GetInputParameter2(dat,dat2,"Minimum radius of movement restrictions","%lf",(void *) &(P.MoveRestrRadius),1,1,0)) P.MoveRestrRadius=0;
	if(!GetInputParameter2(dat,dat2,"Movement restrictions start time","%lf",(void *) &(P.MoveRestrTimeStartBase),1,1,0)) P.MoveRestrTimeStartBase=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat,dat2,"Impose blanket movement restrictions","%i",(void *) &(P.DoBlanketMoveRestr),1,1,0)) P.DoBlanketMoveRestr=0;
	if(P.DoAdUnits)
		{
		if(!GetInputParameter2(dat,dat2,"Movement restrictions in administrative units rather than rings","%i",(void *) &(P.MoveRestrByAdminUnit),1,1,0)) P.MoveRestrByAdminUnit=0;
		if(!GetInputParameter2(dat,dat2,"Administrative unit divisor for movement restrictions","%i",(void *) &(P.MoveRestrAdminUnitDivisor),1,1,0)) P.MoveRestrAdminUnitDivisor=1;
		if((P.MoveRestrAdminUnitDivisor==0)||(P.MoveRestrByAdminUnit==0)) P.MoveRestrAdminUnitDivisor=1;
		}
	else
		{P.MoveRestrAdminUnitDivisor=1;P.MoveRestrByAdminUnit=0;}
	
	if (!GetInputParameter2(dat, dat2, "Include delay to case detection", "%i", (void*)&(P.DoDetectDelay), 1, 1, 0)) P.DoDetectDelay = 0;
	if (!GetInputParameter2(dat, dat2, "Do clustered case detection by household", "%i", (void*)&(P.DoClusterCaseDetection), 1, 1, 0)) P.DoClusterCaseDetection = 0;
	if (P.DoDetectDelay)
	{
		if(!GetInputParameter2(dat, dat2, "Mean detection delay", "%lf", (void*)&(P.DetectTime), 1, 1, 0)) P.DetectTime=0;
	}
	//Case detection - horrible and hard coded at the moment: ggilani 03/02/15
	//if(!GetInputParameter2(dat,dat2,"Include case detection","%i",(void *) &(P.DoCaseDetection),1,1,0)) P.DoCaseDetection=0;
	//if(P.DoCaseDetection)
	//{
		//if(!GetInputParameter2(dat,dat2,"Do case detection by admin unit","%i",(void *) &(P.DoCaseDetectionAdunit),1,1,0)) P.DoCaseDetectionAdunit=0;
		//if((P.DoAdUnits)&&(P.DoCaseDetectionAdunit))
		//{
		//	double AdunitCaseDetection[MAX_ADUNITS];

		//	if(!GetInputParameter(dat,dat2,"Case detection rate per admin unit","%lf",(void *) AdunitCaseDetection,P.NumAdunits,1,0))
		//	{
		//		for(i=0;i<P.NumAdunits;i++) AdunitCaseDetection[i]=0;
		//	}
		//	for(i=0;i<P.NumAdunits;i++)
		//	{
				//Keep this as legacy code which related to setting different case detection rates for different countries, i.e. for West Africa
				//if((AdUnits[i].id==P.AltAU1)||(AdUnits[i].id==P.AltAU2)) //detection rate scaling for specific adunits
				//{
				//	AdUnits[i].caseDetectRate=AdunitCaseDetection[i]*P.RRAlt;
				//}
				//else if((int)(AdUnits[i].id/P.CountryDivisor)==P.TargetCountry) //detection rate scaling for guinea
				//{
				//	AdUnits[i].caseDetectRate=AdunitCaseDetection[i]*P.RR1;
				//	AdUnits[i].caseDetectInit = AdunitCaseDetection[i]*P.RR1; //add to allow us to update case detection
				//}
				//else if((int)(AdUnits[i].id/P.CountryDivisor)==P.TargetCountry2) //detection rate scaling for liberia
				//{
				//	AdUnits[i].caseDetectRate=AdunitCaseDetection[i]*P.RR2;
				//}
				//else if((int)(AdUnits[i].id/P.CountryDivisor)==P.TargetCountry3) //detection rate scaling for sierra leone
				//{
				//	AdUnits[i].caseDetectRate=AdunitCaseDetection[i]*P.RR3;
				//}
				//else
				//{
				//	AdUnits[i].caseDetectRate=AdunitCaseDetection[i];
				//}
			//	AdUnits[i].caseDetectRate = AdunitCaseDetection[i] * P.RR1;
			//	AdUnits[i].caseDetectInit = AdunitCaseDetection[i] * P.RR1; //This is to do with updating case ascertainment, as for above.

			//}
		//}
		//else
		//{
		//	if(!GetInputParameter2(dat,dat2,"Case detection rate","%lf",(void *) &(P.CaseDetectionRate),1,1,0)) P.CaseDetectionRate=1;
		//}
		//if (!GetInputParameter2(dat, dat2, "Update case detection", "%i", (void *) &(P.DoUpdateCaseDetection), 1, 1, 0)) P.DoUpdateCaseDetection = 0;
		//if (P.DoUpdateCaseDetection)
		//{
		//	if (!GetInputParameter2(dat, dat2, "Update case detection by time", "%i", (void*)&(P.DoUpdateCaseDetectionByTime), 1, 1, 0)) P.DoUpdateCaseDetectionByTime = 0;
		//	if(P.DoUpdateCaseDetectionByTime)
		//	{
		//		GetInputParameter(dat, dat2, "Number of change points to update case detection", "%i", (void*)&(P.NUpdateCaseDetection), 1, 1, 0);
		//		GetInputParameter(dat, dat2, "Change points to update case detection", "%lf", (void*)&(P.TimeToUpdateCaseDetection), P.NUpdateCaseDetection, 1, 0);
		//		//GetInputParameter(dat, dat2, "Updated case detection", "%lf", (void *) &(P.UpdatedCaseDetectionRate), 1, 1, 0);
		//		!GetInputParameter2(dat, dat2, "List of update case detections", "%lf", (void*)P.ListUpdateCaseDetection, P.NUpdateCaseDetection, 1, 0);
		//		P.ListUpdateCaseDetection[0] = P.RR2;
		//		P.ListUpdateCaseDetection[1] = P.RR3;
		//	}
		//	if (!GetInputParameter2(dat, dat2, "Update case detection by case numbers", "%i", (void*)&(P.DoUpdateCaseDetectionByCases), 1, 1, 0)) P.DoUpdateCaseDetectionByCases = 0;
		//	if(P.DoUpdateCaseDetectionByCases)
		//	{
		//		P.UpdateCaseDetectionByCasesFlag = 0;
		//		GetInputParameter(dat, dat2, "Threshold case number to update case detection rate", "%i", (void*)&(P.CaseThresholdUntilUpdateCaseDetection), 1, 1, 0);
		//		GetInputParameter(dat, dat2, "Case detection rate after threshold reached", "%lf", (void*)&(P.CaseDetectionRateAfterThresholdReached), 1, 1, 0);
		//	}
		//}
		//if (!GetInputParameter2(dat, dat2, "Do clustered case detection", "%i", (void *) &(P.DoClusterCaseDetection), 1, 1, 0)) P.DoClusterCaseDetection = 0;
		//some parameters related to delay between onset and ascertainment
		//if (!GetInputParameter2(dat, dat2, "Include delay to case detection", "%i", (void*)&(P.DoDetectDelay), 1, 1, 0)) P.DoDetectDelay = 0;
		//if (P.DoDetectDelay)
		//{
		//	GetInputParameter(dat, dat2, "Mean detection delay", "%lf", (void*)&(P.DetectTime), 1, 1, 0);
		//	if (!GetInputParameter2(dat, dat2, "Delay to detection inverse CDF", "%lf", (void*)P.detect_icdf, CDF_RES + 1, 1, 0))
		//	{
		//		P.detect_icdf[CDF_RES] = 1e10;
		//		for (i = 0; i < CDF_RES; i++)
		//			P.detect_icdf[i] = -log(1 - ((double)i) / CDF_RES);
		//	}
		//	for (i = 0; i <= CDF_RES; i++)
		//		P.detect_icdf[i] = exp(-P.detect_icdf[i]);
		//}
	//}

	//Parameters for pseudo-contact tracing effect
	if(!GetInputParameter2(dat,dat2,"Include contact tracing","%i",(void *) &(P.DoContactTracing),1,1,0)) P.DoContactTracing=0;
	if(P.DoContactTracing)
	{
		if (!GetInputParameter2(dat, dat2, "Contact tracing trigger incidence per cell", "%lf", (void*)&(P.ContactTracingCellIncThresh), 1, 1, 0)) P.ContactTracingCellIncThresh = 1000000000;
		if (!GetInputParameter2(dat, dat2, "Contact tracing start time", "%lf", (void*)&(P.ContactTracingTimeStartBase), 1, 1, 0)) P.ContactTracingTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
		GetInputParameter(dat,dat2,"Relative infectiousness of contact traced case","%lf",(void *) &(P.RelativeInfectiousnessContactTraced),1,1,0);
		GetInputParameter(dat,dat2,"Duration of contact tracing","%lf",(void *) &(P.contactTraceDuration),1,1,0);
		//GetInputParameter(dat,dat2,"Cumulative cases before contact tracing begins","%i",(void *) &(P.contactTraceCaseThreshold),1,1,0);
		//if(!GetInputParameter2(dat,dat2,"Cumulative cases when contact tracing increases","%i",(void *) &(P.contactTraceCaseThresholdInc),1,1,0)) P.contactTraceCaseThresholdInc=1e9;
		//couple of flags to help with contact tracing
		P.GlobalContractTracingStarted=P.GlobalContractTracingIncreased=0;
		if(P.DoAdUnits)
		{
			//int AdunitCTCapacity[MAX_ADUNITS];
			//int AdunitCTCapacityInc[MAX_ADUNITS];
			//int AdunitCTThreshold[MAX_ADUNITS];
			if (!GetInputParameter(dat, dat2, "Contact tracing capacity per admin unit", "%i", (void*)&(P.AdunitCTCapacity), 1, 1, 0)) P.AdunitCTCapacity = 0;
			if (!GetInputParameter(dat, dat2, "Contact tracing threshold per admin unit", "%i", (void*)&(P.AdunitCTThreshold), 1, 1, 0)) P.AdunitCTThreshold = 0;
			if (!GetInputParameter2(dat, dat2, "Contact tracing increased capacity per admin unit", "%i", (void*)&(P.AdunitCTCapacityInc), 1, 1, 0)) P.AdunitCTCapacityInc = 0;
			if (!GetInputParameter2(dat, dat2, "Contact tracing capacity reached before increasing teams", "%lf", (void*)&(P.CapacityToMoreCT), 1, 1, 0)) P.CapacityToMoreCT = 1;
			if (!GetInputParameter2(dat, dat2, "Subsequent time to increase teams", "%lf", (void*)&(P.DelayToCT), 1, 1, 0)) P.DelayToCT = 0;

			for(i=0;i<P.NumAdunits;i++)
			{
				//Some terrible hard coding to assign different thresholds to Guinea and Liberia&Sierra Leone!! Replace as soon as possible!! ggilani: 19/11/14
				if((int)(AdUnits[i].id/P.CountryDivisor)==P.TargetCountry) //i.e. if the admin unit is in Guinea, assign contact tracing threshold 1
				{
					AdUnits[i].contactTraceCaseThreshold=P.CT_scale1*P.AdunitCTThreshold;
					AdUnits[i].contactTraceCapacity=P.CT_scale1*P.AdunitCTCapacity; //scaling up contact tracing capacity if needed
					AdUnits[i].contactTraceCapacityInc=P.CTinc_scale1*P.AdunitCTCapacityInc; //scaling up increased contact tracing capacity if needed
					AdUnits[i].nextTimeToCT = 0;
				}
				else //else assign contact tracing threshold 2
				{
					AdUnits[i].contactTraceCaseThreshold=P.CT_thresh2*P.AdunitCTThreshold;
					AdUnits[i].contactTraceCapacity=P.CT_scale2*P.AdunitCTCapacity; //scaling up contact tracing capacity if needed
					AdUnits[i].contactTraceCapacityInc=P.CTinc_scale2*P.AdunitCTCapacityInc;
				}
				//also set the flag for whether contact tracing has begun to zero, and set day on which contact tracing begins to zero for each admin unit
				//AdUnits[i].contactTraceThresholdCrossed=0;
				//AdUnits[i].contactTraceStartDay=1e9;
				//AdUnits[i].contactTraceCurrent=0; //added this to keep track of total number of people being contact traced at any given time - ggilani 07/06/17
			}
		}
		else
		{
			if(!GetInputParameter2(dat,dat2,"Countrywide capacity for contact tracing","%i",(void *) &(P.contactTraceCapacity),1,1,0)) P.contactTraceCapacity=0;
			if(!GetInputParameter(dat,dat2,"Cumulative cases per country before contact tracing begins","%i",(void *) &(P.contactTraceCaseThreshold),1,1,0)) P.contactTraceCaseThreshold=1;
		}
		if(P.DoContactTracing) //  check this!
		{
			if(!GetInputParameter2(dat,dat2,"Proportion of contacts to trace","%lf",(void *) &(P.propContactTraced),1,1,0)) P.propContactTraced=1; //so if we don't specify this, everyone will be contact traced
			if(!GetInputParameter2(dat, dat2, "Proportion of contacts lost to follow up", "%lf", (void*)&(P.propContactLost), 1, 1, 0)) P.propContactLost = 0;
			if(!GetInputParameter2(dat,dat2,"Time to hospitalisation for contact traced case","%lf",(void *) &(P.HospitalisationTime_contactTrace),1,1,0)) P.HospitalisationTime_contactTrace=1;
			if (!GetInputParameter2(dat, dat2, "Number of times to hospitalisation contact traced", "%i", (void *) &(P.NMeanTimeToHospCT), 1, 1, 0)) P.NMeanTimeToHosp = 0;
			if (P.NMeanTimeToHospCT > 0)
			{
				if (P.NMeanTimeToHospCT >= MAX_CHANGE_POINTS) ERR_CRITICAL("MAX_CHANGE_POINTS too small\n");
				GetInputParameter(dat, dat2, "Change points times to hospitalisation contact traced", "%lf", (void *)P.ChangePointMeanTimeToHospCT, P.NMeanTimeToHospCT, 1, 0);
				GetInputParameter(dat, dat2, "Times to hospitalisation contact traced", "%lf", (void *)P.MeanTimeToHospCT, P.NMeanTimeToHospCT, 1, 0);
			}
			P.CurrIndMeanTimeToHospCT = 0;
		}
	}
	//Moved the number of rings to here as it should be included for both ring vaccination and new contact tracing: ggilani 06/06/17
	if(P.DoRingVaccination||P.DoContactTracing)
	{
		if(!GetInputParameter2(dat,dat2,"Number of rings to vaccinate","%i",(void *) &(P.NVaccRings),1,1,0)) P.NVaccRings=1;
		if(P.VaccRingScale>1&&P.NVaccRings==1)
		{
			P.NVaccRings=P.NVaccRings*P.VaccRingScale;
		}
		if (!GetInputParameter2(dat, dat2, "Probability of establishing vaccination ring/contact tracing", "%lf", (void*)&(P.ProbEstablishRing), 1, 1, 0)) P.ProbEstablishRing = 1;
	}

	//Capital City effect
	if(!GetInputParameter2(dat,dat2,"Include capital city effect","%i",(void *) &(P.DoCapitalCityEffect),1,1,0)) P.DoCapitalCityEffect=0;
	if((P.DoCapitalCityEffect)&&(P.DoAdUnits)) //requires admin units
	{
		//read in target capital cities
		if(!GetInputParameter2(dat,dat2,"Capital city admin unit","%i",(void *) &(P.CapitalCityAdunit),1,1,0)) P.CapitalCityAdunit=0;
		if(!GetInputParameter2(dat,dat2,"Capital city admin unit 2","%i",(void *) &(P.CapitalCityAdunit2),1,1,0)) P.CapitalCityAdunit2=0;
		if(!GetInputParameter2(dat,dat2,"Capital city admin unit 3","%i",(void *) &(P.CapitalCityAdunit3),1,1,0)) P.CapitalCityAdunit3=0;
		if(!GetInputParameter2(dat,dat2,"Do capital city distance effect","%i",(void *) &(P.DoCapitalCityDistanceEffect),1,1,0)) P.DoCapitalCityDistanceEffect=0;
		if(!GetInputParameter2(dat,dat2,"Capital city distance effect","%lf",(void *) &(P.CapitalCityDistanceEffect),1,1,0)) P.CapitalCityDistanceEffect=1.0;
		if(!GetInputParameter2(dat,dat2,"Do capital city population effect","%i",(void *) &(P.DoCapitalCityPopEffect),1,1,0)) P.DoCapitalCityPopEffect=0;
		if(!GetInputParameter2(dat,dat2,"Capital city population effect","%lf",(void *) &(P.CapitalCityPopEffect),1,1,0)) P.CapitalCityPopEffect=1.0;
		if(!GetInputParameter2(dat,dat2,"Do capital city additive effect","%i",(void *) &(P.DoCapitalCityAddEffect),1,1,0)) P.DoCapitalCityAddEffect=0;
		//if(!GetInputParameter2(dat,dat2,"Capital city additive effect","%lf",(void *) &(P.CapitalCityAddEffect),1,1,0)) P.CapitalCityAddEffect=0.0;

	}

	if(!GetInputParameter2(dat,dat2,"Trigger incidence per cell for place closure","%i",(void *) &(P.PlaceCloseCellIncThresh),1,1,0)) P.PlaceCloseCellIncThresh=1000000000;
	if(!GetInputParameter2(dat,dat2,"Delay to start place closure","%lf",(void *) &(P.PlaceCloseDelayMean),1,1,0)) P.PlaceCloseDelayMean=0;
	if(!GetInputParameter2(dat,dat2,"Duration of place closure","%lf",(void *) &(P.PlaceCloseDurationBase),1,1,0)) P.PlaceCloseDurationBase=7;
	if(!GetInputParameter2(dat,dat2,"Duration of second place closure","%lf",(void *) &(P.PlaceCloseDuration2),1,1,0)) P.PlaceCloseDuration2=7;
	if(P.DoPlaces) 
		{
		if(!GetInputParameter2(dat,dat2,"Proportion of places remaining open after closure by place type","%lf",(void *) P.PlaceCloseEffect,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++) P.PlaceCloseEffect[i]=1;
		}
	if(P.DoHouseholds)
		if(!GetInputParameter2(dat,dat2,"Relative houshold contact rate after closure","%lf",(void *) &P.PlaceCloseHouseholdRelContact,1,1,0)) P.PlaceCloseHouseholdRelContact=1;
	if(!GetInputParameter2(dat,dat2,"Relative spatial contact rate after closure","%lf",(void *) &P.PlaceCloseSpatialRelContact,1,1,0))
		{
		P.PlaceCloseSpatialRelContact=1;
		}
	if(!GetInputParameter2(dat,dat2,"Minimum radius for place closure","%lf",(void *) &(P.PlaceCloseRadius),1,1,0)) P.PlaceCloseRadius=0;
	if(!GetInputParameter2(dat,dat2,"Place closure start time","%lf",(void *) &(P.PlaceCloseTimeStartBase),1,1,0)) P.PlaceCloseTimeStartBase=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat,dat2,"Place closure second start time","%lf",(void *) &(P.PlaceCloseTimeStartBase2),1,1,0)) P.PlaceCloseTimeStartBase2=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat,dat2,"Places close only once","%i",(void *) &(P.DoPlaceCloseOnceOnly),1,1,0)) P.DoPlaceCloseOnceOnly=0;
	if(!GetInputParameter2(dat,dat2,"Place closure incidence threshold","%i",(void *) &(P.PlaceCloseIncTrig),1,1,0)) P.PlaceCloseIncTrig=1;
	if(!GetInputParameter2(dat,dat2,"Place closure fractional incidence threshold","%lf",(void *) &(P.PlaceCloseFracIncTrig),1,1,0)) P.PlaceCloseFracIncTrig=0;
	if((P.DoAdUnits)&&(P.DoPlaces))
		{
		if(!GetInputParameter2(dat,dat2,"Place closure in administrative units rather than rings","%i",(void *) &(P.PlaceCloseByAdminUnit),1,1,0)) P.PlaceCloseByAdminUnit=0;
		if(!GetInputParameter2(dat,dat2,"Administrative unit divisor for place closure","%i",(void *) &(P.PlaceCloseAdminUnitDivisor),1,1,0)) P.PlaceCloseAdminUnitDivisor=1;
		if(!GetInputParameter2(dat,dat2,"Place types to close for admin unit closure (0/1 array)","%i",(void *) &(P.PlaceCloseAdunitPlaceTypes),P.PlaceTypeNum,1,0))
			for(i=0;i<P.PlaceTypeNum;i++) P.PlaceCloseAdunitPlaceTypes[i]=0;
		if(!GetInputParameter2(dat,dat2,"Cumulative proportion of place members needing to become sick for admin unit closure","%lf",(void *) &(P.PlaceCloseCasePropThresh),1,1,0)) P.PlaceCloseCasePropThresh=2;
		if(!GetInputParameter2(dat,dat2,"Proportion of places in admin unit needing to pass threshold for place closure","%lf",(void *) &(P.PlaceCloseAdunitPropThresh),1,1,0)) P.PlaceCloseAdunitPropThresh=2;
		if((P.PlaceCloseAdminUnitDivisor<1)||(P.PlaceCloseByAdminUnit==0)) P.PlaceCloseAdminUnitDivisor=1;
		}
	else
		{P.PlaceCloseAdminUnitDivisor=1;P.PlaceCloseByAdminUnit=0;}
	
	if(!GetInputParameter2(dat,dat2,"Trigger incidence per cell for social distancing","%i",(void *) &(P.SocDistCellIncThresh),1,1,0)) P.SocDistCellIncThresh=1000000000;
	if(!GetInputParameter2(dat,dat2,"Duration of social distancing","%lf",(void *) &(P.SocDistDuration),1,1,0)) P.SocDistDuration=7;
	if(P.DoPlaces)
		{
		if(!GetInputParameter2(dat,dat2,"Relative place contact rate given social distancing by place type","%lf",(void *) P.SocDistPlaceEffect,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++) P.SocDistPlaceEffect[i]=1;
		if(!GetInputParameter2(dat,dat2,"Relative place contact rate given enhanced social distancing by place type","%lf",(void *) P.ESocDistPlaceEffect,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++) P.ESocDistPlaceEffect[i]=1;
		}
	if(P.DoHouseholds)
		{
		if(!GetInputParameter2(dat,dat2,"Relative houshold contact rate given social distancing","%lf",(void *) &P.SocDistHouseholdEffect,1,1,0)) P.SocDistHouseholdEffect=1;
		if(!GetInputParameter2(dat,dat2,"Relative houshold contact rate given enhanced social distancing","%lf",(void *) &P.ESocDistHouseholdEffect,1,1,0)) P.ESocDistHouseholdEffect=1;
		}
	if(!GetInputParameter2(dat,dat2,"Relative spatial contact rate given social distancing","%lf",(void *) &P.SocDistSpatialEffect,1,1,0)) P.SocDistSpatialEffect=1;
	if(!GetInputParameter2(dat,dat2,"Minimum radius for social distancing","%lf",(void *) &(P.SocDistRadius),1,1,0)) P.SocDistRadius=0;
	if(!GetInputParameter2(dat,dat2,"Social distancing start time","%lf",(void *) &(P.SocDistTimeStartBase),1,1,0)) P.SocDistTimeStartBase=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat,dat2,"Proportion compliant with enhanced social distancing","%lf",(void *) &(P.ESocProportionCompliant),1,1,0)) P.ESocProportionCompliant=0;
	if(!GetInputParameter2(dat,dat2,"Relative spatial contact rate given enhanced social distancing","%lf",(void *) &P.ESocDistSpatialEffect,1,1,0)) P.ESocDistSpatialEffect=1;


	if(!GetInputParameter2(dat,dat2,"Airport closure effectiveness","%lf",(void *) &(P.AirportCloseEffectiveness),1,1,0)) P.AirportCloseEffectiveness=0;
	P.AirportCloseEffectiveness=1.0-P.AirportCloseEffectiveness;
	if(!GetInputParameter2(dat,dat2,"Airport closure start time","%lf",(void *) &(P.AirportCloseTimeStartBase),1,1,0)) P.AirportCloseTimeStartBase=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat,dat2,"Airport closure duration","%lf",(void *) &(P.AirportCloseDuration),1,1,0)) P.AirportCloseDuration=USHRT_MAX/P.TimeStepsPerDay;

	if(P.DoHouseholds)
		{
		if(!GetInputParameter2(dat,dat2,"Household quarantine start time","%lf",(void *) &(P.HQuarantineTimeStartBase),1,1,0)) P.HQuarantineTimeStartBase=USHRT_MAX/P.TimeStepsPerDay;
		if(!GetInputParameter2(dat,dat2,"Delay to start household quarantine","%lf",(void *) &(P.HQuarantineHouseDelay),1,1,0)) P.HQuarantineHouseDelay=0;
		if(!GetInputParameter2(dat,dat2,"Length of time households are quarantined","%lf",(void *) &(P.HQuarantineHouseDuration),1,1,0)) P.HQuarantineHouseDuration=0;
		if(!GetInputParameter2(dat,dat2,"Duration of household quarantine policy","%lf",(void *) &(P.HQuarantinePolicyDuration),1,1,0)) P.HQuarantinePolicyDuration=USHRT_MAX/P.TimeStepsPerDay;
		if(!GetInputParameter2(dat,dat2,"Relative household contact rate after quarantine","%lf",(void *) &(P.HQuarantineHouseEffect),1,1,0)) P.HQuarantineHouseEffect=1;
		if(P.DoPlaces)
			{
			if(!GetInputParameter2(dat,dat2,"Residual place contacts after household quarantine by place type","%lf",(void *) P.HQuarantinePlaceEffect,P.PlaceTypeNum,1,0))
				for(i=0;i<NUM_PLACE_TYPES;i++) P.HQuarantinePlaceEffect[i]=1;
			}
		if(!GetInputParameter2(dat,dat2,"Residual spatial contacts after household quarantine","%lf",(void *) &(P.HQuarantineSpatialEffect),1,1,0)) P.HQuarantineSpatialEffect=1;
		if(!GetInputParameter2(dat,dat2,"Household level compliance with quarantine","%lf",(void *) &(P.HQuarantinePropHouseCompliant),1,1,0)) P.HQuarantinePropHouseCompliant=1;
		if(!GetInputParameter2(dat,dat2,"Individual level compliance with quarantine","%lf",(void *) &(P.HQuarantinePropIndivCompliant),1,1,0)) P.HQuarantinePropIndivCompliant=1;
		}
	else
		P.HQuarantineTimeStartBase=1e10;
	if(!GetInputParameter2(dat,dat2,"Case isolation start time","%lf",(void *) &(P.CaseIsolationTimeStartBase),1,1,0)) P.CaseIsolationTimeStartBase=USHRT_MAX/P.TimeStepsPerDay;
	if(!GetInputParameter2(dat,dat2,"Proportion of detected cases isolated","%lf",(void *) &(P.CaseIsolationProp),1,1,0)) P.CaseIsolationProp=0;
	if(!GetInputParameter2(dat,dat2,"Delay to start case isolation","%lf",(void *) &(P.CaseIsolationDelay),1,1,0)) P.CaseIsolationDelay=0;
	if(!GetInputParameter2(dat,dat2,"Duration of case isolation","%lf",(void *) &(P.CaseIsolationDuration),1,1,0)) P.CaseIsolationDuration=0;
	if(!GetInputParameter2(dat,dat2,"Duration of case isolation policy","%lf",(void *) &(P.CaseIsolationPolicyDuration),1,1,0)) P.CaseIsolationPolicyDuration=1e10;
	if(!GetInputParameter2(dat,dat2,"Residual contacts after case isolation","%lf",(void *) &(P.CaseIsolationEffectiveness),1,1,0)) P.CaseIsolationEffectiveness=1;
	if(P.DoHouseholds)
		{
		if(!GetInputParameter2(dat,dat2,"Residual household contacts after case isolation","%lf",(void *) &(P.CaseIsolationHouseEffectiveness),1,1,0))
			P.CaseIsolationHouseEffectiveness=P.CaseIsolationEffectiveness;
		}
	if(P.DoPlaces)
		{
		if(!GetInputParameter2(dat,dat2,"Number of key workers randomly distributed in the population","%i",(void *) &(P.KeyWorkerPopNum),1,1,0)) P.KeyWorkerPopNum=0;
		if(!GetInputParameter2(dat,dat2,"Number of key workers in different places by place type","%i",(void *) P.KeyWorkerPlaceNum,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++) P.KeyWorkerPlaceNum[i]=0;
		if(!GetInputParameter2(dat,dat2,"Proportion of staff who are key workers per chosen place by place type","%lf",(void *) P.KeyWorkerPropInKeyPlaces,P.PlaceTypeNum,1,0))
			for(i=0;i<NUM_PLACE_TYPES;i++) P.KeyWorkerPropInKeyPlaces[i]=1.0;
		if(!GetInputParameter2(dat,dat2,"Trigger incidence per cell for key worker prophylaxis","%i",(void *) &(P.KeyWorkerProphCellIncThresh),1,1,0)) P.KeyWorkerProphCellIncThresh=1000000000;
		if(!GetInputParameter2(dat,dat2,"Key worker prophylaxis start time","%lf",(void *) &(P.KeyWorkerProphTimeStartBase),1,1,0)) P.KeyWorkerProphTimeStartBase=USHRT_MAX/P.TimeStepsPerDay;
		if(!GetInputParameter2(dat,dat2,"Duration of key worker prophylaxis","%lf",(void *) &(P.KeyWorkerProphDuration),1,1,0)) P.KeyWorkerProphDuration=0;
		if(!GetInputParameter2(dat,dat2,"Time interval from start of key worker prophylaxis before policy restarted","%lf",(void *) &(P.KeyWorkerProphRenewalDuration),1,1,0)) P.KeyWorkerProphRenewalDuration=P.KeyWorkerProphDuration;
		if(P.DoHouseholds)
			{
			if(!GetInputParameter2(dat,dat2,"Proportion of key workers whose households are also treated as key workers","%lf",(void *) &(P.KeyWorkerHouseProp),1,1,0)) P.KeyWorkerHouseProp=0;
			}
		if(!GetInputParameter2(dat,dat2,"Minimum radius for key worker prophylaxis","%lf",(void *) &(P.KeyWorkerProphRadius),1,1,0)) P.KeyWorkerProphRadius=0;
		}
	else
		{		
		P.KeyWorkerPopNum=0;
		P.KeyWorkerProphTimeStartBase=1e10;
		}
	if(!GetInputParameter2(dat,dat2,"Initial rate of importation of infections","%lf",(void *) &(P.InfectionImportRate1),1,1,0)) P.InfectionImportRate1=0;
	if(!GetInputParameter2(dat,dat2,"Changed rate of importation of infections","%lf",(void *) &(P.InfectionImportRate2),1,1,0)) P.InfectionImportRate2=0;
	if(!GetInputParameter2(dat,dat2,"Time when infection rate changes","%lf",(void *) &(P.InfectionImportChangeTime),1,1,0)) P.InfectionImportChangeTime=1e10;
	if(!GetInputParameter2(dat,dat2,"Imports via air travel","%i",(void *) &(P.DoImportsViaAirports),1,1,0)) P.DoImportsViaAirports=0;
	if(!GetInputParameter2(dat,dat2,"Length of importation time profile provided","%i",(void *) &(P.DurImportTimeProfile),1,1,0)) P.DurImportTimeProfile=0;
	if(P.DurImportTimeProfile>0)
		{
		if(P.DurImportTimeProfile>=MAX_DUR_IMPORT_PROFILE) ERR_CRITICAL("MAX_DUR_IMPORT_PROFILE too small\n");
		GetInputParameter(dat,dat2,"Daily importation time profile","%lf",(void *) P.ImportInfectionTimeProfile,P.DurImportTimeProfile,1,0);
		}
	//added some more code to allow us to specify where to import to
	if(!GetInputParameter2(dat,dat2,"Import to specific location","%i",(void *) &(P.DoImportToSpecLocation),1,1,0)) P.DoImportToSpecLocation=0;
	if(!GetInputParameter2(dat,dat2,"Location for imported cases","%lf",(void *) &(P.ImportLocation[0]),2,1,0))
	{
		P.ImportLocation[0]=-1000;
		P.ImportLocation[1]=-1000;
	}
	//Added this to parameter list so that recording infection events (and the number to record) can easily be turned off and on: ggilani - 10/10/2014
	if(!GetInputParameter2(dat,dat2,"Record infection events","%i",(void *) &(P.DoRecordInfEvents),1,1,0)) P.DoRecordInfEvents=0;
	if(P.DoRecordInfEvents)
	{
		if(!GetInputParameter2(dat,dat2,"Max number of infection events to record","%i",(void *) &(P.MaxInfEvents),1,1,0)) P.MaxInfEvents=1000;
		if (!GetInputParameter2(dat, dat2, "Record infection events per run", "%i", (void*) & (P.RecordInfEventsPerRun), 1, 1, 0)) P.RecordInfEventsPerRun = 0;
	}
	else
	{
		P.MaxInfEvents=0;
	}
	//Include a limit to the number of infections to simulate, if this happens before time runs out
	if(!GetInputParameter2(dat,dat2,"Limit number of infections","%i",(void *) &(P.LimitNumInfections),1,1,0)) P.LimitNumInfections=0;
	if(P.LimitNumInfections)
	{
		if(!GetInputParameter2(dat,dat2,"Max number of infections","%i",(void *) &(P.MaxNumInfections),1,1,0)) P.MaxNumInfections=60000;
	}
	//Add cross border contact parameter
	if(!GetInputParameter2(dat,dat2,"Proportion of infections allowed across country borders","%lf", (void *) &(P.PropCrossBorderInf),1,1,0)) P.PropCrossBorderInf=1;
	if(P.BC_scale!=1.0)
	{
		P.PropCrossBorderInf*=P.BC_scale; //scale cross border contact if necessary
	}
	//Add origin-destination matrix parameter
	if(!GetInputParameter2(dat,dat2,"Output origin destination matrix","%i", (void *) &(P.DoOriginDestinationMatrix),1,1,0)) P.DoOriginDestinationMatrix=0;
	//Some parameters relating to road networks: ggilani 12/02/15
	if(P.DoRoadNetwork)
	{
		if(!GetInputParameter2(dat,dat2,"Do distance road effect","%i", (void *) &(P.DoRoadDistanceEffect),1,1,0)) P.DoRoadDistanceEffect=0;
		if(!GetInputParameter2(dat,dat2,"Distance scaling for road accessibility","%lf", (void *) &(P.RoadAccessDistance),1,1,0)) P.RoadAccessDistance=1;
		if(!GetInputParameter2(dat,dat2,"Do population road effect","%i", (void *) &(P.DoRoadPopEffect),1,1,0)) P.DoRoadPopEffect=0;
		if(!GetInputParameter2(dat,dat2,"Population scaling for road accessibility","%lf", (void *) &(P.RoadAccessPop),1,1,0)) P.RoadAccessPop=1;
		if(!GetInputParameter2(dat,dat2,"Maximum road type to consider","%i", (void *) &(P.MaxRoadType),1,1,0)) P.MaxRoadType=2; //2 corresponds to primary roads in the road network file
		if(!GetInputParameter2(dat,dat2,"Maximum cells neighbouring road","%i", (void *) &(P.MaxRoadNeighbour),1,1,0)) P.MaxRoadNeighbour=2;
	}
	
	fclose(dat);

	if(P.DoOneGen!=0) P.DoOneGen=1;
	P.ColourPeriod=2000;
	P.MoveRestrRadius2=P.MoveRestrRadius*P.MoveRestrRadius;
	P.SocDistRadius2=P.SocDistRadius*P.SocDistRadius;
	P.VaccRadius2=P.VaccRadius*P.VaccRadius;
	if (!P.LimitGeoVaccDosesPerCase)
	{
		P.VaccRadiusHighDensity2 = P.VaccRadiusHighDensity * P.VaccRadiusHighDensity;
	}
	P.VaccMinRadius2=P.VaccMinRadius*P.VaccMinRadius;
	P.TreatRadius2=P.TreatRadius*P.TreatRadius;
	P.PlaceCloseRadius2=P.PlaceCloseRadius*P.PlaceCloseRadius;
	P.KeyWorkerProphRadius2=P.KeyWorkerProphRadius*P.KeyWorkerProphRadius;
	if(P.TreatRadius2==0) P.TreatRadius2=-1;
	if(P.VaccRadius2==0) P.VaccRadius2=-1;
	if(P.PlaceCloseRadius2==0) P.PlaceCloseRadius2=-1;
	if(P.MoveRestrRadius2==0) P.MoveRestrRadius2=-1;
	if(P.SocDistRadius2==0) P.SocDistRadius2=-1;
	if(P.KeyWorkerProphRadius2==0) P.KeyWorkerProphRadius2=-1;
	if(P.TreatCellIncThresh<1) P.TreatCellIncThresh=1;
	if(P.MoveRestrCellIncThresh<1) P.MoveRestrCellIncThresh=1;
	if(P.PlaceCloseCellIncThresh<1) P.PlaceCloseCellIncThresh=1;
	if(P.KeyWorkerProphCellIncThresh<1) P.KeyWorkerProphCellIncThresh=1;
	P.usHQuarantineHouseDuration=((unsigned short int) (P.HQuarantineHouseDuration*P.TimeStepsPerDay));
	P.usVaccTimeToEfficacy=((unsigned short int) (P.VaccTimeToEfficacy*P.TimeStepsPerDay));
	P.usVaccTimeToEfficacyThirdRing=((unsigned short int) (P.VaccTimeToEfficacyThirdVaccRing*P.TimeStepsPerDay));
	P.usVaccTimeEfficacySwitch=((unsigned short int) (P.VaccTimeEfficacySwitch*P.TimeStepsPerDay));
	P.usCaseIsolationDelay=((unsigned short int) (P.CaseIsolationDelay*P.TimeStepsPerDay));
	P.usCaseIsolationDuration=((unsigned short int) (P.CaseIsolationDuration*P.TimeStepsPerDay));
	P.usCaseAbsenteeismDuration=((unsigned short int) (P.CaseAbsenteeismDuration*P.TimeStepsPerDay));
	P.usCaseAbsenteeismDelay=((unsigned short int) (P.CaseAbsenteeismDelay*P.TimeStepsPerDay));
	if(P.DoUTM_coords)
		{
		for(i=0;i<=1000;i++)
			{
			asin2sqx[i]=asin(sqrt(((double) (i))/1000));
			asin2sqx[i]=asin2sqx[i]*asin2sqx[i];
			}
		for(t=0;t<=360;t++)
			{
			sinx[(int) t]=sin(PI*t/180);
			cosx[(int) t]=cos(PI*t/180);
			}
		}
	fprintf(stderr,"Parameters read\n");
}

void ReadInterventions(char *IntFile)
{
	FILE *dat,*dat2;
	double r,s,t,dt,startt,stopt;
	int n,i,j,k,au,ni,f,nsr;
	char buf[65536],txt[65536];
	intervention CurInterv;

	fprintf(stderr,"Reading intervention file.\n");
	if(!(dat=fopen(IntFile,"r"))) ERR_CRITICAL("Unable to open intervention file\n");
	fscanf(dat,"%*[^<]"); // needs to be separate line because start of file
	fscanf(dat,"<%[^>]",txt);
	if(strcmp(txt,"\?xml version=\"1.0\" encoding=\"ISO-8859-1\"\?")!=0) ERR_CRITICAL("Intervention file not XML.\n");
	fscanf(dat,"%*[^<]<%[^>]",txt);
	if(strcmp(txt,"InterventionSettings")!=0) ERR_CRITICAL("Intervention has no top level.\n");
	ni=0;
	while(!feof(dat))
		{
		fscanf(dat,"%*[^<]<%[^>]",txt);
		if(strcmp(txt,"intervention")==0)
			{
			ni++;
			fscanf(dat,"%*[^<]<%[^>]",txt);
			if(strcmp(txt,"parameters")!=0) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			if(!GetXMLNode(dat,"Type","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			if(strcmp(txt,"Treatment")==0)
				CurInterv.InterventionType=0;
			else if(strcmp(txt,"Vaccination")==0)
				CurInterv.InterventionType=1;
			else if(strcmp(txt,"ITN")==0)
				CurInterv.InterventionType=2;
			else if(strcmp(txt,"IRS")==0)
				CurInterv.InterventionType=3;
			else if(strcmp(txt,"GM")==0)
				CurInterv.InterventionType=4;
			else if(strcmp(txt,"MSAT")==0)
				CurInterv.InterventionType=5;
			else
				sscanf(txt,"%i",&CurInterv.InterventionType);
			if(!GetXMLNode(dat,"AUThresh","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt,"%i",&CurInterv.DoAUThresh);
			if(!GetXMLNode(dat,"StartTime","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt,"%lg",&CurInterv.StartTime);
			startt=CurInterv.StartTime;
			if(!GetXMLNode(dat,"StopTime","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt,"%lg",&CurInterv.StopTime);
			stopt=CurInterv.StopTime;
			if(!GetXMLNode(dat,"MinDuration","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt,"%lg",&CurInterv.MinDuration);
			CurInterv.MinDuration*=DAYS_PER_YEAR;
			if(!GetXMLNode(dat,"RepeatInterval","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt,"%lg",&CurInterv.RepeatInterval);
			CurInterv.RepeatInterval*=DAYS_PER_YEAR;
			if(!GetXMLNode(dat,"MaxPrevAtStart","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt,"%lg",&CurInterv.StartThresholdHigh);
			if(!GetXMLNode(dat,"MinPrevAtStart","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt,"%lg",&CurInterv.StartThresholdLow);
			if(!GetXMLNode(dat,"MaxPrevAtStop","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt,"%lg",&CurInterv.StopThreshold);
			if(GetXMLNode(dat,"NoStartAfterMinDur","parameters",txt,1))
				sscanf(txt,"%i",&CurInterv.NoStartAfterMin);
			else
				CurInterv.NoStartAfterMin=0;
			if(!GetXMLNode(dat,"Level","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt,"%lg",&CurInterv.Level);
			if(GetXMLNode(dat,"LevelCellVar","parameters",txt,1))
				sscanf(txt,"%lg",&CurInterv.LevelCellVar);
			else
				CurInterv.LevelCellVar=0;
			if(GetXMLNode(dat,"LevelAUVar","parameters",txt,1))
				sscanf(txt,"%lg",&CurInterv.LevelAUVar);
			else
				CurInterv.LevelCellVar=0;
			if(GetXMLNode(dat,"LevelCountryVar","parameters",txt,1))
				sscanf(txt,"%lg",&CurInterv.LevelCountryVar);
			else
				CurInterv.LevelCellVar=0;
			if(GetXMLNode(dat,"LevelClustering","parameters",txt,1))
				sscanf(txt,"%lg",&CurInterv.LevelClustering);
			else
				CurInterv.LevelClustering=0;
			if(GetXMLNode(dat,"ControlParam","parameters",txt,1))
				sscanf(txt,"%lg",&CurInterv.ControlParam);
			else
				CurInterv.ControlParam=0;
			if(GetXMLNode(dat,"TimeOffset","parameters",txt,1))
				sscanf(txt,"%lg",&CurInterv.TimeOffset);
			else
				CurInterv.TimeOffset=0;

			if(!GetXMLNode(dat,"MaxRounds","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt,"%u",&CurInterv.MaxRounds);
			if(!GetXMLNode(dat,"MaxResource","parameters",txt,1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt,"%u",&CurInterv.MaxResource);
			if(GetXMLNode(dat,"NumSequentialReplicas","parameters",txt,1))
				sscanf(txt,"%i",&nsr);
			else
				nsr=0;
			do {fscanf(dat,"%*[^<]<%[^>]",txt);} while((strcmp(txt,"/intervention")!=0)&&(strcmp(txt,"/parameters")!=0)&&(!feof(dat)));
			if(strcmp(txt,"/parameters")!=0) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			fscanf(dat,"%*[^<]<%[^>]",txt);
			if((strcmp(txt,"adunits")!=0)&&(strcmp(txt,"countries")!=0)) ERR_CRITICAL("Incomplete adunits/countries specification in intervention file\n");
			if(strcmp(txt,"adunits")==0)
				{
				while(GetXMLNode(dat,"A","adunits",buf,0))
					{
					sscanf(buf,"%s",txt);
					j=atoi(txt);
					if(j==0)
						{
						f=1;au=-1;
						do
							{au++;f=strcmp(txt,AdUnits[au].ad_name);}
						while ((f)&&(au<P.NumAdunits));
						if(!f) 
							{
							r=fabs(CurInterv.Level)+(2.0*ranf()-1)*CurInterv.LevelAUVar;
							if((CurInterv.Level<1)&&(r>1)) 
								r=1;
							else if(r<0)
								r=0;
							for(k=0;k<=nsr;k++)
								{
								AdUnits[au].InterventionList[AdUnits[au].NI]=CurInterv;
								AdUnits[au].InterventionList[AdUnits[au].NI].Level=r;
								AdUnits[au].InterventionList[AdUnits[au].NI].StartTime=startt+((double) k)*(stopt-startt);
								AdUnits[au].InterventionList[AdUnits[au].NI].StopTime=stopt+((double) k)*(stopt-startt);
								AdUnits[au].NI++;
								}
							}
						}
					else
						{
						k=(j%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor;
						au=P.AdunitLevel1Lookup[k];
						if((au>=0)&&(AdUnits[au].id/P.AdunitLevel1Divisor==j/P.AdunitLevel1Divisor))
							{
							r=CurInterv.Level+(2.0*ranf()-1)*CurInterv.LevelAUVar;
							if((CurInterv.Level<1)&&(r>1)) 
								r=1;
							else if(r<0)
								r=0;
							for(k=0;k<=nsr;k++)
								{
								AdUnits[au].InterventionList[AdUnits[au].NI]=CurInterv;
								AdUnits[au].InterventionList[AdUnits[au].NI].Level=r;
								AdUnits[au].InterventionList[AdUnits[au].NI].StartTime=startt+((double) k)*(stopt-startt);
								AdUnits[au].InterventionList[AdUnits[au].NI].StopTime=stopt+((double) k)*(stopt-startt);
								AdUnits[au].NI++;
								}
							}
						}
					}
				}
			else
				{
				while(GetXMLNode(dat,"C","countries",buf,0)) 
					{
					s=(2.0*ranf()-1)*CurInterv.LevelCountryVar;
					sscanf(buf,"%s",txt);
					j=atoi(txt);
					for(au=0;au<P.NumAdunits;au++)
						if(((j==0)&&(strcmp(txt,AdUnits[au].cnt_name)==0))||((j>0)&&(j==AdUnits[au].cnt_id)))
							{
							r=CurInterv.Level+(2.0*ranf()-1)*CurInterv.LevelAUVar+s;
							if((CurInterv.Level<1)&&(r>1)) 
								r=1;
							else if(r<0)
								r=0;
							for(k=0;k<=nsr;k++)
								{
								AdUnits[au].InterventionList[AdUnits[au].NI]=CurInterv;
								AdUnits[au].InterventionList[AdUnits[au].NI].Level=r;
								AdUnits[au].InterventionList[AdUnits[au].NI].StartTime=startt+((double) k)*(stopt-startt);
								AdUnits[au].InterventionList[AdUnits[au].NI].StopTime=stopt+((double) k)*(stopt-startt);
								AdUnits[au].NI++;
								}
							}
					}
				}
			fscanf(dat,"%*[^<]<%[^>]",txt);
			if(strcmp(txt,"/intervention")!=0) ERR_CRITICAL("Incorrect intervention specification in intervention file\n");
			}
		}
	if(strcmp(txt,"/InterventionSettings")!=0) ERR_CRITICAL("Intervention has no top level closure.\n");
	fprintf(stderr,"%i interventions read\n",ni);
	fclose(dat);
}

int GetXMLNode(FILE *dat,char *NodeName,char *ParentName,char *Value,int ResetFilePos)
{
// ResetFilePos=1 leaves dat cursor in same position as when function was called. 0 leaves it at end of NodeName closure
// GetXMLNode returns 1 if NodeName found, 0 otherwise. If NodeName not found, ParentName closure must be

	char buf[65536],CloseNode[2048],CloseParent[2048];
	int CurPos,ret;

	sprintf(CloseParent,"/%s",ParentName);
	CurPos=ftell(dat);
	do 
		{fscanf(dat,"%*[^<]<%[^>]",buf);}
	while((strcmp(buf,CloseParent)!=0)&&(strcmp(buf,NodeName)!=0)&&(!feof(dat)));
	if(strcmp(buf,CloseParent)==0)
		ret=0;
	else
		{
		if(strcmp(buf,NodeName)!=0) ERR_CRITICAL("Incomplete node specification in XML file\n");
		fscanf(dat,">%[^<]",buf);
		if(strlen(buf)<2048) strcpy(Value,buf);
//		fprintf(stderr,"# %s=%s\n",NodeName,Value);
		fscanf(dat,"<%[^>]",buf);
		sprintf(CloseNode,"/%s",NodeName);
		if(strcmp(buf,CloseNode)!=0) ERR_CRITICAL("Incomplete node specification in XML file\n");
		ret=1;
		}
	if(ResetFilePos) fseek(dat,CurPos,0);
	return ret;
}

void ReadAirTravel(char *AirTravelFile)
{
	int i,j,k,l;
	float sc,t,t2;
	float *buf;
	double traf;
	char outname[1024];
	FILE *dat;

	fprintf(stderr,"Reading airport data...\nAirports with no connections = ");
	if(!(dat=fopen(AirTravelFile,"r"))) ERR_CRITICAL("Unable to open airport file\n");
	fscanf(dat,"%i %i",&P.Nairports,&P.Air_popscale);
	sc=((double)P.N)/((double) P.Air_popscale);
	if(P.Nairports>MAX_AIRPORTS) ERR_CRITICAL("Too many airports\n");
	if(P.Nairports<2) ERR_CRITICAL("Too few airports\n");
	if(!(buf=(float *) calloc(P.Nairports+1,sizeof(float)))) ERR_CRITICAL("Unable to allocate airport storage\n");
	if(!(Airports=(airport *) calloc(P.Nairports,sizeof(airport)))) ERR_CRITICAL("Unable to allocate airport storage\n");
	for(i=0;i<P.Nairports;i++)
		{
		fscanf(dat,"%f %f %lg",&(Airports[i].loc_x),&(Airports[i].loc_y),&traf);
		traf*=(P.AirportTrafficScale*sc);
		if((Airports[i].loc_x<P.SpatialBoundingBox[0])||(Airports[i].loc_x>P.SpatialBoundingBox[2])
			||(Airports[i].loc_y<P.SpatialBoundingBox[1])||(Airports[i].loc_y>P.SpatialBoundingBox[3]))
			{
			Airports[i].loc_x=Airports[i].loc_y=-1;
			Airports[i].total_traffic=0;
			}
		else
			{
			//fprintf(stderr,"(%f,%f) ",Airports[i].loc_x,Airports[i].loc_y);
			Airports[i].loc_x-=P.SpatialBoundingBox[0];
			Airports[i].loc_y-=P.SpatialBoundingBox[1];
			Airports[i].total_traffic=traf;
			}
		t=0;
		for(j=k=0;j<P.Nairports;j++)
			{
			fscanf(dat,"%f",buf+j);
			if(buf[j]>0) {k++;t+=buf[j];}
			}
		Airports[i].num_connected=k;
		if(Airports[i].num_connected>0)
			{
			if(!(Airports[i].prop_traffic=(float *) calloc(Airports[i].num_connected,sizeof(float)))) ERR_CRITICAL("Unable to allocate airport storage\n");
			if(!(Airports[i].conn_airports =(unsigned short int *) calloc(Airports[i].num_connected,sizeof(unsigned short int)))) ERR_CRITICAL("Unable to allocate airport storage\n");
			for(j=k=0;j<P.Nairports;j++)
				if(buf[j]>0)
					{
					Airports[i].conn_airports[k]=j;
					Airports[i].prop_traffic[k]=buf[j]/t;
					k++;
					}
			}
		else
			{
			if(Airports[i].total_traffic>0)
				fprintf(stderr,"#%i# ",i);
			else
				fprintf(stderr,"%i ",i);
			}
		}
	fclose(dat);
	free(buf);
	fprintf(stderr,"\nAirport data read OK.\n");
	for(i=0;i<P.Nairports;i++)
		{
/*		fprintf(stderr,"(%f %i|",Airports[i].total_traffic,Airports[i].num_connected);
*/		t=0;k=0;
		for(j=Airports[i].num_connected-1;j>=0;j--)
			{
			if((Airports[i].prop_traffic[j]>0)&&(Airports[Airports[i].conn_airports[j]].total_traffic==0))
				{
				t+=Airports[i].prop_traffic[j];
				Airports[i].num_connected--;
				if(j<Airports[i].num_connected)
					{
					Airports[i].prop_traffic[j]=Airports[i].prop_traffic[Airports[i].num_connected];
					Airports[i].conn_airports[j]=Airports[i].conn_airports[Airports[i].num_connected];
					}
				Airports[i].prop_traffic[Airports[i].num_connected]=0;
				Airports[i].conn_airports[Airports[i].num_connected]=0;
				}
			else if(Airports[i].prop_traffic[j]>0)
				k=1;
			}
/*		fprintf(stderr,"%f %i ",t,k);
*/		t=1.0-t;
		if(k)
			{
			Airports[i].total_traffic*=t;
			t2=0;
			for(j=0;j<Airports[i].num_connected;j++)
				{
				Airports[i].prop_traffic[j]=t2+Airports[i].prop_traffic[j];
				t2=Airports[i].prop_traffic[j];
				}
			for(j=0;j<Airports[i].num_connected;j++)
				Airports[i].prop_traffic[j]/=t2;
/*			if((Airports[i].num_connected>0)&&(Airports[i].prop_traffic[Airports[i].num_connected-1]!=1))
				fprintf(stderr,"<%f> ",Airports[i].prop_traffic[Airports[i].num_connected-1]);
*/			}
		else
			{Airports[i].total_traffic=0;Airports[i].num_connected=0;}
		if(Airports[i].num_connected>0)
			{
			for(j=k=0;k<128;k++)
				{
				t=((double) k)/128;
				while(Airports[i].prop_traffic[j]<t) j++;
				Airports[i].Inv_prop_traffic[k]=j;
				}
			Airports[i].Inv_prop_traffic[128]=Airports[i].num_connected-1;
			}
/*		fprintf(stderr,"%f) ",Airports[i].total_traffic);
*/		}
	fprintf(stderr,"Airport data clipped OK.\n");
	for(i=0;i<MAX_DIST;i++) AirTravelDist[i]=0;
	for(i=0;i<P.Nairports;i++)
		if(Airports[i].total_traffic>0)
			{
			for(j=0;j<Airports[i].num_connected;j++)
				{
				k=(int) Airports[i].conn_airports[j];
				traf=floor(sqrt(dist2_raw(Airports[i].loc_x,Airports[i].loc_y,Airports[k].loc_x,Airports[k].loc_y))/OUTPUT_DIST_SCALE);
				l=(int) traf;
				//fprintf(stderr,"%(%i) ",l);
				if(l<MAX_DIST)
					AirTravelDist[l]+=Airports[i].total_traffic*Airports[i].prop_traffic[j];	
				}
			}
	sprintf(outname,"%s.airdist.xls",OutFile);
	if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open air travel output file\n");
	fprintf(dat,"dist\tfreq\n");
	for(i=0;i<MAX_DIST;i++)
		fprintf(dat,"%i\t%lg\n",i,AirTravelDist[i]);
	fclose(dat);
}

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
					sscanf(buf,"%lg,%lg,%lg,%i,%i",&x,&y,&t,&i2,&l);
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
			i=(int) (((double) P.N)*ranf(refseed));
			if(Hosts[i].keyworker)
				l++;
			else
				{
				Hosts[i].keyworker=1;
				m++;
				P.KeyWorkerNum++;
				P.KeyWorkerIncHouseNum++;
				l=0;
				if(ranf(ranf_seed)<P.KeyWorkerHouseProp)
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
				k=(int) (((double) P.Nplace[j])*ranf(ranf_seed));
				for(i2=0;(m<P.KeyWorkerPlaceNum[j])&&(i2<Places[j][k].n);i2++)
					{
					i=Places[j][k].members[i2];
					if((i<0)||(i>=P.N)) fprintf(stderr,"## %i # ",i);
					if((Hosts[i].keyworker)||(ranf(ranf_seed)>=P.KeyWorkerPropInKeyPlaces[j]))
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
							if((!Hosts[j2].keyworker)&&(ranf(ranf_seed)<P.KeyWorkerHouseProp))
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
			Hosts[i].infectiousness=P.AgeInfectiousness[HOST_AGE_GROUP(i)]*gen_gamma_mt(ranf_seed, P.InfectiousnessGamA,P.InfectiousnessGamR,tn); //made this multi-threaded: 28/11/14
			//Hosts[i].infectiousness=P.AgeInfectiousness[HOST_AGE_GROUP(i)]*gen_beta_mt(P.InfectiousnessBetaA,P.InfectiousnessBetaB,tn);
		q=P.ProportionSymptomatic[HOST_AGE_GROUP(i)];
		if(ranf_mt(ranf_seed, tn)<q) //made this multi-threaded: 28/11/14
			Hosts[i].infectiousness=-P.SymptInfectiousness*Hosts[i].infectiousness;
		j=(int) floor((q=ranf_mt(ranf_seed, tn)*CDF_RES)); //made this multi-threaded: 28/11/14
		q-=((double) j);
		Hosts[i].recovery_time=(unsigned short int) floor(0.5-(P.InfectiousPeriod*log(q*P.infectious_icdf[j+1]+(1.0-q)*P.infectious_icdf[j])/P.TimeStep));
		//adding a step here to see, if we are doing funeral transmission, the host dies: ggilani 14/11/14
		if(P.DoFuneralTransmission)
		{
			if(P.DoMortality)
			{
				if (P.DoEventMortality)
				{
					if (ranf_mt(ranf_seed, tn) <= P.RecoveryProb[(int)ceil(Hosts[i].recovery_time * P.TimeStep)]) Hosts[i].to_die = 0; //made this multi-threaded: 28/11/14
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
					if (ranf_mt(ranf_seed, tn) < probMort)
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
			Hosts[i].hospitalised = (ranf_mt(ranf_seed,tn) < P.PropHospSeekPreOutbreak);
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
	if(P.OutputBitmap)
		{
		InitBMHead();
#ifdef FRESSCA
		if(P.OutputBitmap>=1)
			{
			sprintf(buf,"%s.kml",OutFile);
			KMLFile = fopen(buf,"w");
			}
#endif
		}
#ifdef FRESSCA
	SetupDistrNet();
#endif
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


#ifdef FRESSCA
void SetupDistrNet(void)
{
	int i,j,k,l;
	double x,y;
	char filename[1024];

  //Setting Up Distribution Network
  	if(P.DoDistributionVaccination)
		{
    	fprintf(stderr,"Setting up Distribution Network\n");
    	VaccineDistributionNetwork = new GEO_DISTR_NETWORK(DistribNetworkFile,((int) (P.SampleTime/P.TimeStep))+1);
    	fprintf(stderr,"\nNumber of health centers = %d\n",VaccineDistributionNetwork->get_num_centers(HEALTH));
#ifdef FRESSCA_DEBUG
    	VaccineDistributionNetwork->Print(stderr);
#endif
		VaccineDistributionNetwork->InitGE(KMLFile);
    // Assign each microcell a health center at the bottom of the network
     //   fprintf(stderr,"\nP.Spatial = %f %f",P.SpatialBoundingBox[0],P.SpatialBoundingBox[1]);
		for(k=0;k<P.NMCP;k++)
			{
			l=McellLookup[k]-Mcells;
			if((P.DistribNetCountry==-1)||(Mcells[l].country==P.DistribNetCountry))
				{
				x=P.mcwidth*(((double) (l/P.nmch))+0.5)+P.SpatialBoundingBox[0];
				y=P.mcheight*(((double) (l%P.nmch))+0.5)+P.SpatialBoundingBox[1];
	     		Mcells[l].VaccineCenter = VaccineDistributionNetwork->find_closest_health_center((float) y,(float) x);
             					//fprintf(stderr,"\nCenter Name = %s",Mcells[l].VaccineCenter->get_name());
	     		Mcells[l].VaccineCenter->addDemanders(1);
#ifdef FRESSCA_DEBUG
           			fprintf(stderr,"MCell[%d] x = %f y = %f Center = %s\n",l,y,x,Mcells[l].VaccineCenter->get_name());
#endif
				}
			else
				Mcells[l].VaccineCenter=0;
    		}	
    	int count = 0;
		k=VaccineDistributionNetwork->get_num_centers(HEALTH);
		fprintf(stderr,"Assigned %i health centers to microcells\n",k);
 		MCellDistrIndex = (int **)calloc(k,sizeof(int*));
		int *mca;
		if(!(mca=(int *) calloc(P.NMCP,sizeof(int)))) ERR_CRITICAL("Unable to allocate host storage\n");
  		for(i=0;i<k;i++)
			{
			j=VaccineDistributionNetwork->get_center_from_level_index(HEALTH,i)->get_ndemanders();
			count+=j;
    		MCellDistrIndex[i] = mca;
			mca+=j;
  			}	
		MCellDistrIndex[k]=mca;
     	fprintf(stderr,"\nNumber of demanders = %d Number of MCells = %d\n",count,P.NMCP);
 		int *iIndex = (int*)calloc(k+1,sizeof(int));
		int *nCtr = (int*)calloc(k+1,sizeof(int));
   		for(i=0;i<P.NMCP;i++)
			{ 
			l=McellLookup[i]-Mcells;
			if(!Mcells[l].VaccineCenter)
				{
       			MCellDistrIndex[k][iIndex[k]] = l; 
       			iIndex[k]++; 
				}
			else
				{
       			j = Mcells[l].VaccineCenter->get_id_level();
				if(j>k-1) fprintf(stderr,"*** Uh Oh - %i %i\n",j,k);
       			MCellDistrIndex[j][iIndex[j]] = l; 
       			iIndex[j]++; 
				nCtr[j]+=Mcells[l].n;
				}
   			}
		noCenterMcellNum=iIndex[k];
		for(i=0;i<k;i++)
			{
			fprintf(stderr,"Centre %i (%i): %i demanders %i mcells %i people\n",VaccineDistributionNetwork->get_center_from_level_index(HEALTH,i)->get_id(),
			i,VaccineDistributionNetwork->get_center_from_level_index(HEALTH,i)->get_ndemanders(),iIndex[i],nCtr[i]);
			}
  		count = 0;
		int *queuearray;
		unsigned short int *queuearrayt;
		if(!(State.dvacc_queue=(int *) calloc(P.N,sizeof(int)))) ERR_CRITICAL("Unable to allocate host storage\n");
		if(!(State.dvacc_expiry_time=(unsigned short int *) calloc(2*P.N,sizeof(unsigned short int)))) ERR_CRITICAL("Unable to allocate host storage\n");
		queuearray=State.dvacc_queue;
		queuearrayt=State.dvacc_expiry_time;
  		for(i=0;i<P.NMCP;i++)
			{
    		McellLookup[i]->dvacc_queue = queuearray;
			queuearray+=McellLookup[i]->n;
    		McellLookup[i]->dvacc_min_vacc_time =queuearrayt;
    		McellLookup[i]->dvacc_expiry_time =queuearrayt+P.N;
			queuearrayt+=McellLookup[i]->n;
    		McellLookup[i]->dvacc_count = McellLookup[i]->ndvacc_queue = 0;
			}
	}
}
#endif

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
					sscanf(buf,"%lg,%lg,%lg,%i,%i",&x,&y,&t,&i2,&j2);
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
					sscanf(buf,"%lg %lg %lg %i",&x,&y,&t,&i2);
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
		m+=(Mcells[i].n=(int) ignbin(ranf_seed, (long) (P.N-m), s));
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
				s=ranf(ranf_seed);
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
#ifdef NEW_AGE_MODEL
	P.ts_age=-((int) DEMOG_EQUILIB_YEARS*P.TimeStepsPerYear);
#endif
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
			xh=P.mcwidth*(ranf_mt(ranf_seed, tn)+x);
			yh=P.mcheight*(ranf_mt(ranf_seed, tn)+y);
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

#ifdef NEW_AGE_MODEL
	if(!P.DoLoadSnapshot) 
		{
		fprintf(stderr,"Running demographic model for %i years\n",DEMOG_EQUILIB_YEARS);
		P.ts_age=0;
#pragma omp parallel for private(tn,i) schedule(static,1)
		for(tn=0;tn<P.NumThreads;tn++)
		  for(i=tn;i<P.N;i+=P.NumThreads)
			EquilibPersonDemog(i,tn);
		fprintf(stderr,"Demography updated\n");
		}
	for(i=0;i<P.N;i++)
		{
		Hosts[i].init_birth_time=Hosts[i].birth_time;
		AgeDist[HOST_AGE_GROUP(i)]++;
		AgeDist2[((int) Hosts[i].life_expectancy)/(((int) P.TimeStepsPerYear)*AGE_GROUP_WIDTH)]++;
		}
#endif

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
					m+=(Mcells[i].np[j2]=(int) ignbin_mt(ranf_seed, (long) (P.Nplace[j2]-m), s, tn));
				t-=mcell_dens[i]/maxd;
				if(Mcells[i].np[j2]>0)
					{
					if(!(Mcells[i].places[j2]=(int *) malloc(Mcells[i].np[j2]*sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
					x=(double) (i/P.nmch);
					y=(double) (i%P.nmch);
					for(j=0;j<Mcells[i].np[j2];j++)
						{
						s=ranf_mt(ranf_seed, tn);
						xh=P.mcwidth*(ranf_mt(ranf_seed, tn)+x);
						yh=P.mcheight*(ranf_mt(ranf_seed, tn)+y);
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


#ifdef COUNTRY_THAILAND
#define MEAN_CHILD_AGE_GAP 2 /* Mean number of years separating adjacent children */
#define MIN_ADULT_AGE 17
#define MAX_MF_PARTNER_AGE_GAP 7  /*  Max number of years older male partner is than female */
#define MAX_FM_PARTNER_AGE_GAP 2  /* Max number of years younger male partner is than female */
#define MIN_PARENT_AGE_GAP 17 /* Min number of years older a mother is than their child */
#define MAX_PARENT_AGE_GAP 33  /* Max number of years older a mother is than their child */
#define MAX_CHILD_AGE 30
#define ONE_CHILD_TWO_PERS_PROB 0.03
#define TWO_CHILD_THREE_PERS_PROB 0.0
#define TWO_PERS_HOUSE_PROB_OLD 0.4
#define ONE_PERS_HOUSE_PROB_OLD 0.25
#define ONE_PERS_HOUSE_PROB_YOUNG 0
#define TWO_PERS_HOUSE_PROB_YOUNG 0
#define ONE_CHILD_PROB_YOUNGEST_CHILD_UNDER_FIVE 0.3
#define PROB_YOUNGEST_CHILD_UNDER_FIVE 0
#define ZERO_CHILD_THREE_PERS_PROB 0.0
#define ONE_CHILD_FOUR_PERS_PROB 0.1
#define YOUNG_AND_SINGLE 25
#define YOUNG_AND_SINGLE_SLOPE 0
#define NOCHILD_PERS_AGE 52
#define THREE_CHILD_FIVE_PERS_PROB 0.5
#endif
#ifdef COUNTRY_INDONESIA
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
#ifdef COUNTRY_UK
#define MEAN_CHILD_AGE_GAP 2 
#define MIN_ADULT_AGE 19
#define MAX_MF_PARTNER_AGE_GAP 5  
#define MAX_FM_PARTNER_AGE_GAP 5  
#define MIN_PARENT_AGE_GAP 19 
#define MAX_PARENT_AGE_GAP 44  
#define MAX_CHILD_AGE 20
#define ONE_CHILD_TWO_PERS_PROB 0.08
#define TWO_CHILD_THREE_PERS_PROB 0.11
#define ONE_PERS_HOUSE_PROB_OLD 0.5
#define TWO_PERS_HOUSE_PROB_OLD 0.5
#define ONE_PERS_HOUSE_PROB_YOUNG 0.23
#define TWO_PERS_HOUSE_PROB_YOUNG 0.23
#define ONE_CHILD_PROB_YOUNGEST_CHILD_UNDER_FIVE 0.5
#define TWO_CHILDREN_PROB_YOUNGEST_UNDER_FIVE 0.0
#define PROB_YOUNGEST_CHILD_UNDER_FIVE 0
#define ZERO_CHILD_THREE_PERS_PROB 0.25
#define ONE_CHILD_FOUR_PERS_PROB 0.2
#define YOUNG_AND_SINGLE_SLOPE 0.7
#define YOUNG_AND_SINGLE 36
#define NOCHILD_PERS_AGE 44
#define OLD_PERS_AGE 60
#define THREE_CHILD_FIVE_PERS_PROB 0.5
#define OLDER_GEN_GAP 19
#endif
#ifdef COUNTRY_CHINA
#define MEAN_CHILD_AGE_GAP 2 
#define MIN_ADULT_AGE 19
#define MAX_MF_PARTNER_AGE_GAP 5  
#define MAX_FM_PARTNER_AGE_GAP 5  
#define MIN_PARENT_AGE_GAP 19 
#define MAX_PARENT_AGE_GAP 44  
#define MAX_CHILD_AGE 20
#define ONE_CHILD_TWO_PERS_PROB 0.08
#define TWO_CHILD_THREE_PERS_PROB 0.11
#define ONE_PERS_HOUSE_PROB_OLD 0.5
#define TWO_PERS_HOUSE_PROB_OLD 0.5
#define ONE_PERS_HOUSE_PROB_YOUNG 0.23
#define TWO_PERS_HOUSE_PROB_YOUNG 0.23
#define ONE_CHILD_PROB_YOUNGEST_CHILD_UNDER_FIVE 0.5
#define TWO_CHILDREN_PROB_YOUNGEST_UNDER_FIVE 0.0
#define PROB_YOUNGEST_CHILD_UNDER_FIVE 0
#define ZERO_CHILD_THREE_PERS_PROB 0.25
#define ONE_CHILD_FOUR_PERS_PROB 0.2
#define YOUNG_AND_SINGLE_SLOPE 0.7
#define YOUNG_AND_SINGLE 36
#define NOCHILD_PERS_AGE 44
#define OLD_PERS_AGE 60
#define THREE_CHILD_FIVE_PERS_PROB 0.5
#define OLDER_GEN_GAP 19
#endif
#ifdef COUNTRY_US
#define MEAN_CHILD_AGE_GAP 2 
#define MIN_ADULT_AGE 19
#define MAX_MF_PARTNER_AGE_GAP 5  
#define MAX_FM_PARTNER_AGE_GAP 5  
#define MIN_PARENT_AGE_GAP 19 
#define MAX_PARENT_AGE_GAP 45
#define MAX_CHILD_AGE 20
#define ONE_CHILD_TWO_PERS_PROB 0.08
#define TWO_CHILD_THREE_PERS_PROB 0.11
#define ONE_PERS_HOUSE_PROB_OLD 0.42
#define TWO_PERS_HOUSE_PROB_OLD 0.42
#define ONE_PERS_HOUSE_PROB_YOUNG 0.26
#define TWO_PERS_HOUSE_PROB_YOUNG 0.26
#define ONE_CHILD_PROB_YOUNGEST_CHILD_UNDER_FIVE 0.69
#define TWO_CHILDREN_PROB_YOUNGEST_UNDER_FIVE 0.0
#define PROB_YOUNGEST_CHILD_UNDER_FIVE 0
#define ZERO_CHILD_THREE_PERS_PROB 0.15
#define ONE_CHILD_FOUR_PERS_PROB 0.12
#define YOUNG_AND_SINGLE_SLOPE 0.75
#define YOUNG_AND_SINGLE 36
#define NOCHILD_PERS_AGE 45
#define OLD_PERS_AGE 59
#define THREE_CHILD_FIVE_PERS_PROB 0.5
#define OLDER_GEN_GAP 19
#endif
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
			a[i]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];
		}
	else
		{
		if(n==1)
			{
			if(ranf_mt(ranf_seed, tn)<ONE_PERS_HOUSE_PROB_OLD)
				{
				do
					{a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];}
#ifdef COUNTRY_THAILAND
				while(a[0]<NOCHILD_PERS_AGE);
#else
				while((a[0]<NOCHILD_PERS_AGE)
					||(ranf_mt(ranf_seed, tn)>(((double) a[0])-NOCHILD_PERS_AGE+1)/(OLD_PERS_AGE-NOCHILD_PERS_AGE+1)));
#endif
				}
			else if((ONE_PERS_HOUSE_PROB_YOUNG>0)&&(ranf_mt(ranf_seed, tn)<ONE_PERS_HOUSE_PROB_YOUNG/(1-ONE_PERS_HOUSE_PROB_OLD)))
				{
				do
					{a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];}
				while((a[0]>YOUNG_AND_SINGLE)||(a[0]<MIN_ADULT_AGE)
					||(ranf_mt(ranf_seed, tn)>1-YOUNG_AND_SINGLE_SLOPE*(((double) a[0])-MIN_ADULT_AGE)/(YOUNG_AND_SINGLE-MIN_ADULT_AGE)));
				}
			else
				while((a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))])<MIN_ADULT_AGE);
			}
		else if(n==2)
			{
			if(ranf_mt(ranf_seed, tn)<TWO_PERS_HOUSE_PROB_OLD)
				{
				do
					{a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];}
#ifdef COUNTRY_THAILAND
				while(a[0]<NOCHILD_PERS_AGE);
#else
				while((a[0]<NOCHILD_PERS_AGE)
					||(ranf_mt(ranf_seed, tn)>(((double) a[0])-NOCHILD_PERS_AGE+1)/(OLD_PERS_AGE-NOCHILD_PERS_AGE+1)));
#endif
				do
					{a[1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];}
#ifdef COUNTRY_THAILAND
				while((a[1]>=a[0]+MAX_MF_PARTNER_AGE_GAP)||(a[1]<a[0]-MAX_FM_PARTNER_AGE_GAP)||(a[1]<MIN_ADULT_AGE));
#else
				while((a[1]>a[0]+MAX_MF_PARTNER_AGE_GAP)||(a[1]<a[0]-MAX_FM_PARTNER_AGE_GAP)||(a[1]<NOCHILD_PERS_AGE)
					||(ranf_mt(ranf_seed, tn)>(((double) a[1])-NOCHILD_PERS_AGE+1)/(OLD_PERS_AGE-NOCHILD_PERS_AGE+1)));
#endif
				}
			else if(ranf_mt(ranf_seed, tn)<ONE_CHILD_TWO_PERS_PROB/(1-TWO_PERS_HOUSE_PROB_OLD))
				{
				while((a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))])>MAX_CHILD_AGE);
				do
					{a[1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];}
#ifdef COUNTRY_THAILAND
				while((a[1]>=a[0]+MAX_PARENT_AGE_GAP)||(a[1]<a[0]+MIN_PARENT_AGE_GAP));
#else
				while((a[1]>a[0]+MAX_PARENT_AGE_GAP)||(a[1]<a[0]+MIN_PARENT_AGE_GAP)||(a[1]<MIN_ADULT_AGE));
#endif
				}
			else if((TWO_PERS_HOUSE_PROB_YOUNG>0) &&(ranf_mt(ranf_seed,tn)<TWO_PERS_HOUSE_PROB_YOUNG/(1-TWO_PERS_HOUSE_PROB_OLD-ONE_CHILD_TWO_PERS_PROB)))
				{
				do
					{a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];}
				while((a[0]<MIN_ADULT_AGE)||(a[0]>YOUNG_AND_SINGLE)
					||(ranf_mt(ranf_seed, tn)>1-YOUNG_AND_SINGLE_SLOPE*(((double) a[0])-MIN_ADULT_AGE)/(YOUNG_AND_SINGLE-MIN_ADULT_AGE)));
				do
					{a[1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed,tn))];}
#ifdef COUNTRY_THAILAND
				while((a[1]>=a[0]+MAX_MF_PARTNER_AGE_GAP)||(a[1]<a[0]-MAX_FM_PARTNER_AGE_GAP)||(a[1]<MIN_ADULT_AGE));
#else
				while((a[1]>a[0]+MAX_MF_PARTNER_AGE_GAP)||(a[1]<a[0]-MAX_FM_PARTNER_AGE_GAP)||(a[1]<MIN_ADULT_AGE));
#endif
				}
			else
				{
				do
					{a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];}
				while(a[0]<MIN_ADULT_AGE);
				do
					{a[1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];}
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
					nc=(ranf_mt(ranf_seed, tn)<ZERO_CHILD_THREE_PERS_PROB)?0:((ranf_mt(ranf_seed, tn)<TWO_CHILD_THREE_PERS_PROB)?2:1);
				else
					nc=1;
				}
			else if(n==4)
				nc=(ranf_mt(ranf_seed, tn)<ONE_CHILD_FOUR_PERS_PROB)?1:2;
			else if(n==5)
				nc=(ranf_mt(ranf_seed, tn)<THREE_CHILD_FIVE_PERS_PROB)?3:2;
			else
#ifdef COUNTRY_INDONESIA
				do
					{
					nc=ignpoi_mt(ranf_seed, FRAC_CHILDREN_BIG_HOUSEHOLDS*((double) n),tn);
					}
				while((nc>n-2)||(nc<2));
#else
				nc=n-2-(int) (3*ranf_mt(ranf_seed, tn));
#endif
			if(nc==0)
				{
				do
					{
					a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];
					a[1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];
					}
#ifdef COUNTRY_THAILAND
				while((a[0]<MIN_ADULT_AGE)||(a[1]>=a[0]+MAX_PARENT_AGE_GAP)
					||(a[1]<a[0]+MIN_PARENT_AGE_GAP)||(a[1]>NOCHILD_PERS_AGE));
#else
				while((a[1]<MIN_ADULT_AGE)||(a[0]<MIN_ADULT_AGE));
#endif
				do
					{
					a[2]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];
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
					a[i]=a[i-1]+1+((int) ignpoi_mt(ranf_seed, MEAN_CHILD_AGE_GAP-1,tn));
				j=a[nc-1]-(MAX_PARENT_AGE_GAP-MIN_PARENT_AGE_GAP);
				if(j>0)
					j+=MAX_PARENT_AGE_GAP;
				else
					j=MAX_PARENT_AGE_GAP;
				do
					{
					a[nc]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];
					a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];
					k=((nc>1)?a[nc-1]:0)+a[0];
					l=k-MAX_CHILD_AGE;
					if(ranf_mt(ranf_seed, tn)<ONE_CHILD_PROB_YOUNGEST_CHILD_UNDER_FIVE) l=k-5;
					}
				while((l>0)||(a[nc]>a[0]+j)||(a[nc]<k+MIN_PARENT_AGE_GAP));
				for(i=1;i<nc;i++) a[i]+=a[0];
				if((n>nc+1)&&(ranf_mt(ranf_seed, tn)>PROP_OTHER_PARENT_AWAY))
					{
					do
						{a[nc+1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];}
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
						a[i]=a[i-1]+1+((int) ignpoi_mt(ranf_seed, MEAN_CHILD_AGE_GAP-1,tn));
					a[0]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))]-a[(int) (ranf_mt(ranf_seed, tn)*((double) nc))];
					for(i=1;i<nc;i++) a[i]+=a[0];
					k=(((nc==1)&&(ranf_mt(ranf_seed, tn)<ONE_CHILD_PROB_YOUNGEST_CHILD_UNDER_FIVE))||((nc==2)&&(ranf_mt(ranf_seed, tn)<TWO_CHILDREN_PROB_YOUNGEST_UNDER_FIVE))
						||((nc>2)&&(ranf_mt(ranf_seed, tn)<PROB_YOUNGEST_CHILD_UNDER_FIVE)))?5:MAX_CHILD_AGE;
					}
				while((a[0]<0)||(a[0]>k)||(a[nc-1]>MAX_CHILD_AGE));
				j=a[nc-1]-a[0]-(MAX_PARENT_AGE_GAP-MIN_PARENT_AGE_GAP);
				if(j>0)
					j+=MAX_PARENT_AGE_GAP;
				else
					j=MAX_PARENT_AGE_GAP;
				do
					{
					a[nc]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];
					k=a[nc-1];
					l=k-MAX_CHILD_AGE;
					}
				while((a[nc]>a[0]+j)||(a[nc]<k+MIN_PARENT_AGE_GAP)||(a[nc]<MIN_ADULT_AGE));
				if((n>nc+1)&&(ranf_mt(ranf_seed, tn)>PROP_OTHER_PARENT_AWAY))
					{
					do
						{a[nc+1]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))];}
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
						while((a[i]=State.InvAgeDist[ad][(int) (1000.0*ranf_mt(ranf_seed, tn))])<j);
					}
				}
			}
		}
	for(i=0;i<n;i++)
		{
#ifdef NEW_AGE_MODEL
		Hosts[pers+i].birth_time=P.ts_age-(int) ((((double) a[i])+ranf_mt(ranf_seed, tn))*P.TimeStepsPerYear);
		if(P.DoDeath)
			{
			f=P.CumulPropDead[0][HOST_AGE_YEAR(pers+i)];
			f=P.InvLifeExpecDist[0][(int) ceil(1000.0*(f+ranf_mt(ranf_seed, tn)*(0)))];
			Hosts[pers+i].life_expectancy=((unsigned short int) (f*P.TimeStepsPerYear));
//			if(((int) f)<a[i]) fprintf(stderr,"##%lg %i %i %lg %lg## ",f,a[i],HOST_AGE_YEAR(pers+i),P.CumulPropDead[0][HOST_AGE_YEAR(pers+i)],P.InvLifeExpecDist[0][(int) ceil(1000*P.CumulPropDead[0][HOST_AGE_YEAR(pers+i)])]);
			}
#else
		Hosts[pers+i].age=(unsigned short int) a[i];
#pragma omp atomic
		AgeDist[HOST_AGE_GROUP(pers+i)]++;
#endif
		}
}

#define NO_EMP_WEIGHT 100

//void AssignPeopleToPlaces(void)
//{
//	int i,i2,j,j2,k,k2,l,m,m2,tp,f,f2,f3,f4,ic,mx,my,a,cnt,tn,ca,nt,nn;
//	int *PeopleArray;
//	int *NearestPlaces[MAX_NUM_THREADS];
//	double s,t,s2,*NearestPlacesProb[MAX_NUM_THREADS];
//	cell *ct;
//	int dbg;
//	int npt;
//	int nh_assigned,i3; //nh_assigned keeps track of the number of people in each household assigned to the same place, i3 is just another counting variable - ggilani 10/02/17
//	int g,g1,g2,maxph; //added this variables to keep track of place variables when assigning members of the same household to the same place - ggilani 13/02/17
//
//	npt=NUM_PLACE_TYPES;
//
//	if(P.DoPlaces)
//		{
//		fprintf(stderr,"Assigning people to places....\n");
//		for(i=0;i<P.NC;i++)
//			{
//			Cells[i].infected=Cells[i].susceptible;
//			if(!(Cells[i].susceptible=(int *) calloc(Cells[i].n,sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");	
//			Cells[i].cumTC=Cells[i].n;
//			}
///*
// *	PropPlaces initialisation is only valid for non-overlapping places.
// */
//		for(i=0;i<P.N;i++)
//		{
//			for(tp=0;tp<npt;tp++)
//			{
//				Hosts[i].PlaceLinks[tp]=-1;
//				Hosts[i].PlacePresence[tp]=0;
//			}
//		}
//
//		for(tp=0;tp<P.PlaceTypeNum;tp++)
//			if(tp!=HOTEL_PLACE_TYPE)
//				{
//				cnt=0;
//				for(a=0;a<P.NCP;a++)
//					{
//					i=(int) (CellLookup[a]-Cells);
//					Cells[i].n=0;
//					for(j=0;j<Cells[i].cumTC;j++)
//						{
//						k=HOST_AGE_YEAR(Cells[i].members[j]);
//						f=((PropPlaces[k][tp]>0)&&(ranf()<PropPlaces[k][tp]));
//						if(f)
//							for(k=0;(k<tp)&&(f);k++)
//								if(Hosts[Cells[i].members[j]].PlaceLinks[k]>=0) f=0; /*(ranf()<P.PlaceExclusivityMatrix[tp][k]); */
//						/* Am assuming people can only belong to 1 place (and a hotel) at present */
//						if(f)
//							{
//							Cells[i].susceptible[Cells[i].n]=Cells[i].members[j];
//							Cells[i].n++;
//							cnt++;
//							Hosts[Cells[i].members[j]].PlacePresence[tp]=1; //marking that this person does actually belong to this place type - will use this later to ensure we can assign some household members to similar places
//							}
//						}
//					Cells[i].S=Cells[i].n;
//					Cells[i].I=0;
//					}
//				if(!(PeopleArray=(int *) calloc(cnt,sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
//				j2=0;
//				for(a=0;a<P.NCP;a++)
//					{
//					i=(int) (CellLookup[a]-Cells);
//					for(j=0;j<Cells[i].n;j++)
//						{
//						PeopleArray[j2]=Cells[i].susceptible[j];
//						j2++;
//						}
//					}
//				for(i2=0;i2<2;i2++)
//					for(k=0;k<cnt;k++)
//						{
//						while((l=(int) (((double) cnt)*ranf()))==k);
//						j2=PeopleArray[k];
//						PeopleArray[k]=PeopleArray[l];
//						PeopleArray[l]=j2;
//						}
//				m=0;
//				if(tp<P.nsp)
//					{
//					for(i=0;i<P.Nplace[tp];i++)
//						{
//						m+=(int) (Places[tp][i].treat_end_time=(float) Places[tp][i].n);
//						Places[tp][i].n=0;
//						}
//					}
//				else if(P.PlaceTypeSizePower[tp]==0 && P.PlaceTypeSizeSD[tp]==0)
//					{
//					for(i=0;i<P.Nplace[tp];i++)
//						{
//						Places[tp][i].n=0;
//						j=1+((int) ignpoi(P.PlaceTypeMeanSize[tp]-1));
//						if(j>USHRT_MAX-1) j=USHRT_MAX-1;
//						m+=(int) (Places[tp][i].treat_end_time=(float) j);
//						}
//					}
//				//added this code to allow a place size to be specified according to a lognormal distribution - ggilani 09/02/17
//				else if(P.PlaceTypeSizePower[tp]==0 && P.PlaceTypeSizeSD[tp]>0)
//					{
//					for(i=0;i<P.Nplace[tp];i++)
//						{
//						Places[tp][i].n=0;
//						j=(int) gen_lognormal(P.PlaceTypeMeanSize[tp],P.PlaceTypeSizeSD[tp]);
//						if(j>USHRT_MAX-1) j=USHRT_MAX-1;
//						m+=(int) (Places[tp][i].treat_end_time=(float) j);
//						}
//					}
//				else
//					{
//					s=pow(P.PlaceTypeSizeOffset[tp]/(P.PlaceTypeSizeOffset[tp]+P.PlaceTypeSizeMax[tp]-1),P.PlaceTypeSizePower[tp]);
//					for(i=0;i<P.Nplace[tp];i++)
//						{
//						j=(int) floor(P.PlaceTypeSizeOffset[tp]*pow((1-s)*ranf()+s,-1/P.PlaceTypeSizePower[tp])+1-P.PlaceTypeSizeOffset[tp]);
//						if(j>USHRT_MAX-1) j=USHRT_MAX-1;
//						m+=(int) (Places[tp][i].treat_end_time=(float) j);
//						Places[tp][i].n=0;
//						}
//					}
//				if(tp<P.nsp)
//					{
//					t=((double) m)/((double) P.Nplace[tp]);
//					fprintf(stderr,"Adjusting place weights by cell (Capacity=%i Demand=%i  Av place size=%lg)\n",m,cnt,t);
//					for(i=0;i<P.Nplace[tp];i++)
//						if(Places[tp][i].treat_end_time>0)
//							{
//							j=(int) (Places[tp][i].loc_x/P.cwidth);
//							k=j*P.nch+((int) (Places[tp][i].loc_y/P.cheight));
//							Cells[k].I+=(int) Places[tp][i].treat_end_time;
//							}
//					for(k=0;k<P.NC;k++)
//						{
//						i=k%P.nch;
//						j=k/P.nch;
//						f2=Cells[k].I;f3=Cells[k].S;
//						if((i>0)&&(j>0))
//							{f2+=Cells[(j-1)*P.nch+(i-1)].I;f3+=Cells[(j-1)*P.nch+(i-1)].S;}
//						if(i>0)
//							{f2+=Cells[j*P.nch+(i-1)].I;f3+=Cells[j*P.nch+(i-1)].S;}
//						if((i>0)&&(j<P.ncw-1))
//							{f2+=Cells[(j+1)*P.nch+(i-1)].I;f3+=Cells[(j+1)*P.nch+(i-1)].S;}
//						if(j>0)
//							{f2+=Cells[(j-1)*P.nch+i].I;f3+=Cells[(j-1)*P.nch+i].S;}
//						if(j<P.ncw-1)
//							{f2+=Cells[(j+1)*P.nch+i].I;f3+=Cells[(j+1)*P.nch+i].S;}
//						if((i<P.nch-1)&&(j>0))
//							{f2+=Cells[(j-1)*P.nch+(i+1)].I;f3+=Cells[(j-1)*P.nch+(i+1)].S;}
//						if(i<P.nch-1)
//							{f2+=Cells[j*P.nch+(i+1)].I;f3+=Cells[j*P.nch+(i+1)].S;}
//						if((i<P.nch-1)&&(j<P.ncw-1))
//							{f2+=Cells[(j+1)*P.nch+(i+1)].I;f3+=Cells[(j+1)*P.nch+(i+1)].S;}
//						Cells[k].L=f3;Cells[k].R=f2;
//						}
//					m=f2=f3=f4=0;
//					for(k=0;k<P.NC;k++)
//						if((Cells[k].S>0)&&(Cells[k].I==0))
//							{
//							f2+=Cells[k].S;f3++;
//							if(Cells[k].R==0) f4+=Cells[k].S;
//							}
//					fprintf(stderr,"Demand in cells with no places=%i in %i cells\nDemand in cells with no places <=1 cell away=%i\n",f2,f3,f4);
//					for(i=0;i<P.Nplace[tp];i++)
//						if(Places[tp][i].treat_end_time>0)
//							{
//							j=(int) (Places[tp][i].loc_x/P.cwidth);
//							k=j*P.nch+((int) (Places[tp][i].loc_y/P.cheight));
//							if((Cells[k].L>0)&&(Cells[k].R>0))
//								{
//								s=((double) Cells[k].L)/((double) Cells[k].R);
//								Places[tp][i].treat_end_time=ceil(Places[tp][i].treat_end_time*s);
//								}
//							m+=((int) Places[tp][i].treat_end_time);
//							}
//					for(i=0;i<P.NC;i++) Cells[i].L=Cells[i].I=Cells[i].R=0;
//					}
//				t=((double) m)/((double) P.Nplace[tp]);
//				fprintf(stderr,"Adjusting place weights (Capacity=%i Demand=%i  Av place size=%lg)\n",m,cnt,t);
//				for(i=m=0;i<P.Nplace[tp];i++)
//					{
//					s=((double) Places[tp][i].treat_end_time)*43/40-1;
//					m+=(int) (Places[tp][i].treat_end_time=1+((float) ignpoi(s)));
//					}
//				if(tp<P.nsp)
//					s=((double) cnt)*1.075;
//				else
//					s=((double) cnt)*1.125;
//				j2=((int) s)-m;
//				for(i=0;i<j2;i++) Places[tp][(int) (((double) P.Nplace[tp])*ranf())].treat_end_time++;
//				j2=-j2;
//				for(i=0;i<j2;i++)
//					{
//					while(Places[tp][j=(int) (((double) P.Nplace[tp])*ranf())].treat_end_time<2);
//					Places[tp][j].treat_end_time--;
//					}
//				if(P.PlaceTypeNearestNeighb[tp]==0)
//					{
//					for(i=0;i<P.NC;i++) Cells[i].S=0;
//					for(j=0;j<P.Nplace[tp];j++)
//						{
//						i=P.nch*((int) (Places[tp][j].loc_x/P.cwidth))+((int) (Places[tp][j].loc_y/P.cheight));
//						Cells[i].S+=(int) Places[tp][j].treat_end_time;
//						}
//					for(i=0;i<P.NC;i++)
//						{
//						if(Cells[i].S>Cells[i].cumTC)
//							{
//							free(Cells[i].susceptible);
//							if(!(Cells[i].susceptible=(int *) calloc(Cells[i].S,sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
//							}
//						Cells[i].S=0;
//						}
//					for(j=0;j<P.Nplace[tp];j++)
//						{
//						i=P.nch*((int) (Places[tp][j].loc_x/P.cwidth))+((int) (Places[tp][j].loc_y/P.cheight));
//						k=(int) Places[tp][j].treat_end_time;
//						for(j2=0;j2<k;j2++)
//							{
//							Cells[i].susceptible[Cells[i].S]=j;
//							Cells[i].S++;
//							}
//						}
//					}
//				for(i=0;i<P.NumThreads;i++)
//					{
//					if(!(NearestPlaces[i]=(int *) calloc(P.PlaceTypeNearestNeighb[tp]+CACHE_LINE_SIZE,sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
//					if(!(NearestPlacesProb[i]=(double *) calloc(P.PlaceTypeNearestNeighb[tp]+CACHE_LINE_SIZE,sizeof(double)))) ERR_CRITICAL("Unable to allocate cell storage\n");
//					}
//				P.KernelType=P.PlaceTypeKernelType[tp];
//				P.KernelScale=P.PlaceTypeKernelScale[tp];
//				P.KernelShape=P.PlaceTypeKernelShape[tp];
//				P.KernelP3=P.PlaceTypeKernelP3[tp];
//				P.KernelP4=P.PlaceTypeKernelP4[tp];
//				InitKernel(1,1.0);
//				UpdateProbs(1);
//				ca=0;
//				fprintf(stderr,"Allocating people to place type %i\n",tp);
//				a=cnt;
//				nt=P.NumThreads;
//				nn=P.PlaceTypeNearestNeighb[tp];
//				if(P.PlaceTypeNearestNeighb[tp]>0)
//					{
////#pragma omp parallel for private(i,i2,j,j2,k,k2,l,m,m2,f,f2,ic,cnt,tn,s,t,mx,my,i3,nh_assigned) firstprivate(a,nt,nn) reduction(+:ca) schedule(static,1) //added i3, nh_assigned to private
//					for(tn=0;tn<P.NumThreads;tn++)
//						{
//						for(j=tn;j<a;j+=nt)
//							{
//							if(j%1000==0) fprintf(stderr,"(%i) %i      \r",tp,j);
//							for(i2=0;i2<nn;i2++)
//								NearestPlacesProb[tn][i2]=0;
//							l=1;k=m=m2=f2=0;
//							i=PeopleArray[j];
//							ic=Hosts[i].mcell;
//							mx=ic/P.nmch;
//							my=ic%P.nmch;
//							if(Hosts[i].PlaceLinks[tp]<0) //added this so that if any hosts have already be assigned due to their household membership, they will not be reassigned
//								{	
//								while(((k<nn)||(l<4))&&(l<P.nmcw))
//									{
//									if((mx>=0)&&(my>=0)&&(mx<P.nmcw)&&(my<P.nmch))
//										{
//										ic=mx*P.nmch+my;
//										if(Mcells[ic].country==Mcells[Hosts[i].mcell].country)
//											{
//											for(cnt=0;cnt<Mcells[ic].np[tp];cnt++)
//												{
//												if(Mcells[ic].places[tp][cnt]>=P.Nplace[tp]) 
//													fprintf(stderr,"#%i %i %i  ",tp,ic,cnt);
//												t=dist2_raw(Households[Hosts[i].hh].loc_x,Households[Hosts[i].hh].loc_y,
//													Places[tp][Mcells[ic].places[tp][cnt]].loc_x,Places[tp][Mcells[ic].places[tp][cnt]].loc_y);
//												s=numKernel(t);
//												if(tp<P.nsp)
//													{
//													t=((double) Places[tp][Mcells[ic].places[tp][cnt]].treat_end_time);
//													if(HOST_AGE_YEAR(i)<P.PlaceTypeMaxAgeRead[tp])
//														{
//														if((t>0)&&(Places[tp][Mcells[ic].places[tp][cnt]].AvailByAge[HOST_AGE_YEAR(i)]>0))
//															s*=t;
//														else
//															s=0;
//														}
//													else if(t>0)
//														s*=t;
//													}
//												k2=0;j2=0;t=1e10;
//												if(s>0)
//													{
//													if(k<nn)
//														{
//														NearestPlaces[tn][k]=Mcells[ic].places[tp][cnt];
//														NearestPlacesProb[tn][k]=s;
//														k++;
//														}
//													else
//														{
//														for(i2=0;i2<nn;i2++)
//															{
//															if(NearestPlacesProb[tn][i2]<t)
//																{t=NearestPlacesProb[tn][i2];j2=i2;}
//															}
//														if(s>t)
//															{
//															NearestPlacesProb[tn][j2]=s;
//															NearestPlaces[tn][j2]=Mcells[ic].places[tp][cnt];
//															}
//														}
//													}
//												}
//											}
//										}
//									if(m2==0)
//										mx=mx+1;
//									else if(m2==1)
//										my=my-1;
//									else if(m2==2)
//										mx=mx-1;
//									else if(m2==3)
//										my=my+1;
//									f2=(f2+1)%l;
//									if(f2==0)
//										{
//										m2=(m2+1)%4;
//										m=(m+1)%2;
//										if(m==0) l++;
//										}
//									}
//							}
//							s=0;
//							if(k>nn) fprintf(stderr,"*** k>P.PlaceTypeNearestNeighb[tp] ***\n");
//							if(k==0)
//								{
//								fprintf(stderr,"# %i %i     \r",i,j); //Something going wrong here!
//								Hosts[i].PlaceLinks[tp]=-1;
//								}
//							else
//								{
//								for(i2=1;i2<k;i2++)
//									NearestPlacesProb[tn][i2]+=NearestPlacesProb[tn][i2-1];
//								s=NearestPlacesProb[tn][k-1];
//								t=ranf_mt(tn);
//								f=0;
//								for(i2=0;(i2<k)&&(!f);i2++)
//									if((f=(t<NearestPlacesProb[tn][i2]/s)))
//										{
//										Hosts[i].PlaceLinks[tp]=NearestPlaces[tn][i2];
//										ca++;
//
//										//add something about assigning family members here - ggilani 10/02/17
//										if((P.PlaceTypeNum==1)&&(P.PlaceHouseholdOverlap==1))
//										{
//											nh_assigned=1;
//											i3=0;
//											maxph=(int) ceil(Households[Hosts[i].hh].nh/P.PlaceHouseholdDivisor);
//											while((nh_assigned<maxph)&&(i3<Households[Hosts[i].hh].nh)) //so that households are roughly divided between places
//											{
//												if((Hosts[Households[Hosts[i].hh].FirstPerson+i3].PlaceLinks[tp]<0)&&(Hosts[Households[Hosts[i].hh].FirstPerson+i3].PlacePresence[tp]==1))
//												{
//													//then assign this person to the same place
//													Hosts[Households[Hosts[i].hh].FirstPerson+i3].PlaceLinks[tp]=NearestPlaces[tn][i2];
//													nh_assigned++;
//													ca++; //think this should be here - need to check! ggilani 10/02/2017
//												}
//												i3++;
//											}
//										}
//
//										if(tp<P.nsp)
//											{
////#pragma omp critical (places_treat_time)
//											Places[tp][Hosts[i].PlaceLinks[tp]].treat_end_time--;
//											}
//										}
//									if(!f) Hosts[i].PlaceLinks[tp]=-1;
//									if(NearestPlaces[tn][i2]>=P.Nplace[tp]) fprintf(stderr,"@%i %i %i  ",tp,i,j);
//								}
//							}
//						}
//					}
//				else
//					{
//					k2=cnt-ca;
//					m2=cnt;
//					s2=0;
//					a=k2/1000;
//					f=k2;
//					for(ic=0;ic<=30;ic++)
//						{
//						UpdateProbs(1);
//						m2=f-1;
//						if(ic<9)
//							f=100*(9-ic)*a;
//						else if(ic<18)
//							f=10*(18-ic)*a;
//						else if(ic<27)
//							f=(27-ic)*a;
//						else
//							{
//							m2=k2-1;f=0;
//							}
////#pragma omp parallel for private(i,i2,j,k,l,m,f2,f3,t,ct,s,g,g1,g2,i3,nh_assigned) reduction(+:ca) /* schedule(dynamic,500)*/ //add s to private variables, added g,g1,g2,i3 and nh_assigned to private variables
//						for(i2=m2;i2>=f;i2--)
//							{
//							if(i2%1000==0) 
//								fprintf(stderr,"(%i) %i            \r",tp,i2);
//							k=PeopleArray[i2];
//							i=Hosts[k].pcell;
//							f2=1;
//							f3=(HOST_AGE_YEAR(k)>=P.PlaceTypeMaxAgeRead[tp]);
//							if(Hosts[k].PlaceLinks[tp]<0)
//								while((f2>0) && (f2<1000))
//									{
//									do
//										{
//										s=ranf();
//										l=Cells[i].InvCDF[(int) floor(s*1024)];
//										while(Cells[i].cum_trans[l]<s) l++;
//										ct=CellLookup[l];
//										m=(int) (ranf()*((double) ct->S));
//										j=-1;
////#pragma omp critical
//											{
//											if(ct->susceptible[m]>=0)
//												if((f3)||(Places[tp][ct->susceptible[m]].AvailByAge[HOST_AGE_YEAR(k)]>0))
//													{
//													j=ct->susceptible[m];
//													ct->susceptible[m]=-1;
//													}
//											}
//										}
//									while(j<0);
//									t=dist2_raw(Households[Hosts[k].hh].loc_x,Households[Hosts[k].hh].loc_y,Places[tp][j].loc_x,Places[tp][j].loc_y);
//									s=((double) ct->S)/((double) ct->S0)*numKernel(t)/Cells[i].max_trans[l];
//									if((P.DoAdUnits)&&(P.InhibitInterAdunitPlaceAssignment[tp]>0))
//										{
//										if(Mcells[Hosts[k].mcell].adunit!=Mcells[Places[tp][j].mcell].adunit) s*=(1-P.InhibitInterAdunitPlaceAssignment[tp]);
//										}
//									if(ranf()<s)
//										{
////#pragma omp critical
//										l=(--ct->S);
//										if(m<l) ct->susceptible[m]=ct->susceptible[l];
////#pragma omp critical (places_treat_time)
//										Places[tp][j].treat_end_time--;
//										ca++;
//										Hosts[k].PlaceLinks[tp]=j;
//
//										//add something about assigning family members here - ggilani 10/02/17
//										if((P.PlaceTypeNum==1)&&(P.PlaceHouseholdOverlap==1))
//										{
//											nh_assigned=1;
//											i3=0;
//											while((nh_assigned<=(int) ceil(Households[Hosts[k].hh].nh/P.PlaceHouseholdDivisor))&&(i3<Households[Hosts[k].hh].nh)) //so that households are roughly divided between places
//											{
//												//also need to check for space in this scenario???
//												if((Hosts[Households[Hosts[k].hh].FirstPerson+i3].PlaceLinks[tp]<0)&&(Hosts[Households[Hosts[k].hh].FirstPerson+i3].PlacePresence[tp]==1)&&(Places[tp][j].treat_end_time>0))
//												{
//													//need to find a value of m for which place is equal to j
//													g=-1; //this is flag variable
//													g1=m; //g1, g2 are counting variables that allow us to move away from m and find the next closest 
//													g2=m;
//													while(g<0)
//													{
//														g1++;
//														g2--;
//														if(g1<(ct->S))
//														{
//															if(ct->susceptible[g1]==j)
//															{
//																m=g1;
//																g=0;
//															}
//														}
//														else if(g2>=0)
//														{
//															if(ct->susceptible[g2]==j)
//															{
//																m=g2;
//																g=0;
//															}
//														}
//													}
//													//then assign this person to the same place
////#pragma omp critical
//													l=(--ct->S);
//													if(m<l) ct->susceptible[m]=ct->susceptible[l];
//													Hosts[Households[Hosts[k].hh].FirstPerson+i3].PlaceLinks[tp]=j;
//													nh_assigned++;
////#pragma omp critical (places_treat_time)
//													Places[tp][j].treat_end_time--; 
//													ca++; 
//												}
//												i3++;
//											}
//										}
//
//
//
//										f2=0;
//										}
//									else
//										{
//										ct->susceptible[m]=j;
//										f2++;
//										}
//									}
//							}
//						}
//					}
//				fprintf(stderr,"%i hosts assigned to placetype %i\n",ca,tp);
//				free(PeopleArray);
//				for(i=0;i<P.Nplace[tp];i++)
//					{
//					Places[tp][i].treat_end_time=0;
//					Places[tp][i].n=0;
//					}
//				for(i=0;i<P.NumThreads;i++)
//					{
//					free(NearestPlacesProb[i]);
//					free(NearestPlaces[i]);
//					}
//				}
//		for(i=0;i<P.NC;i++)
//			{
//			Cells[i].n=Cells[i].cumTC;
//			Cells[i].cumTC=0;
//			Cells[i].S=Cells[i].I=Cells[i].L=Cells[i].R=0;
//			free(Cells[i].susceptible);
//			Cells[i].susceptible=Cells[i].infected;
//			}
//		P.KernelScale=P.MoveKernelScale;
//		P.KernelShape=P.MoveKernelShape;
//		P.KernelType=P.MoveKernelType;
//		P.KernelP3=P.MoveKernelP3;
//		P.KernelP4=P.MoveKernelP4;
//		}
//}

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
									t = ranf_mt(ranf_seed, tn);
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
							Places[j][i].ng=1+(int) ignpoi_mt(ranf_seed,t,tn);
						if(!(Places[j][i].group_start=(int *) calloc(Places[j][i].ng,sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
						if(!(Places[j][i].group_size=(int *) calloc(Places[j][i].ng,sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
						m=Places[j][i].n-Places[j][i].ng;
						for(k=l=0;k<Places[j][i].ng;k++)
							{
							t=1/((double) (Places[j][i].ng-k));
							Places[j][i].group_start[k]=l;
							Places[j][i].group_size[k]=1+ignbin_mt(ranf_seed, (long) m,t,tn);
							m-=(Places[j][i].group_size[k]-1);
							l+=Places[j][i].group_size[k];
							}
						for(k=0;k<Places[j][i].n;k++)
							{
							l=(int) (((double) Places[j][i].n)*ranf_mt(ranf_seed, tn));
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

				nhcws = (int)(((double)Places[P.HospPlaceTypeNum][i].n / (double)1000) * P.HCWPerThousand);
				nflws = (int)(((double)Places[P.HospPlaceTypeNum][i].n / (double)1000) * P.FLWPerThousand);

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


void InitModel(int run) //passing run number so we can save run number in the infection event log: ggilani - 15/10/2014
{
	int i,j,k,l,m,i2,j2,tn,nim,stt,stp,b,nhcws,nflws;
	double t,s;
	char buf[200];

	if(P.OutputBitmap)
		{
#ifdef WIN32_BM
		if(P.OutputBitmap==1)
			{
			sprintf(buf,"%s.avi",OutFile);
			avi = CreateAvi(buf,P.BitmapMovieFrame,NULL);
			}
#endif
		for(i=0;i<bmh->imagesize;i++) 
			{
			bmi2[i]=bmi3[i]=bmi4[i]=0;
			}
		}

#ifdef NEW_AGE_MODEL
	P.ts_age=0;
#pragma omp parallel for private(i) schedule(static,500)
	for(i=0;i<P.N;i++)
		Hosts[i].birth_time=Hosts[i].init_birth_time;
#endif

	ns=0;
	State.S=P.N;
	State.L=State.I=State.R=0;
	State.cumI=State.cumR=State.cumC=State.cumFC=State.cumFI=State.cumETU=State.cumH=State.cumCT=State.cumCC=State.cumTC=State.cumD=State.cumDC=State.trigDC=State.cumSDB=State.cumDD
		=State.cumInf_h=State.cumInf_n=State.cumInf_s=State.cumHQ
		=State.cumAC=State.cumAH=State.cumAA=State.cumACS
		=State.cumAPC=State.cumAPA=State.cumAPCS=0;
	State.cumT=State.cumUT=State.cumTP=State.cumV=State.sumRad2=State.maxRad2=State.cumV_daily=State.cumVG=0; //added State.cumVG
	State.mvacc_cum=0;
	State.NumBeds = 0;
	State.nringvacc_queue = State.ngeovacc_queue = State.nvacc_queue = State.ringvacc_cum = State.geovacc_cum = State.ringvacc_ind= State.geovacc_ind= State.vacc_ind=State.vacc_cum=0; //reset ring and geo vaccination variables, 21/08/19
#ifdef FRESSCA
	State.dvacc_wastage=0;
#endif
	for(i=0;i<NUM_AGE_GROUPS;i++) State.cumCa[i]=State.cumIa[i]=State.cumDa[i]= State.cumDCa[i] = State.cumETUa[i] = State.cumHa[i] = State.cumVa[i] = 0; //adding det case, hosp, vacc by age: ggilani 22/02/22
	for(i=0;i<P.EvolResistNumTypes;i++) State.cumC_resist[i]=State.cumI_resist[i]=State.cumT_resist[i]=0;
	for(i=0;i<2;i++) State.cumC_keyworker[i]=State.cumI_keyworker[i]=State.cumT_keyworker[i]=0;
	for(i=0;i<NUM_PLACE_TYPES;i++) State.NumPlacesClosed[i]=0;
	for(i=0;i<INFECT_TYPE_MASK;i++) State.cumItype[i]=0;
	//initialise cumulative case counts per country to zero: ggilani 12/11/14
	for(i=0;i<MAX_COUNTRIES;i++) State.cumC_country[i]=0;
	if(P.DoAdUnits)
		for(i=0;i<=P.NumAdunits;i++) 
			{
			State.cumI_adunit[i]=State.cumC_adunit[i]=State.cumT_adunit[i]=State.cumETU_adunit[i]= State.cumH_adunit[i]=State.cumDC_adunit[i]= State.cumDD_adunit[i] = State.cumSDB_adunit[i]=State.cumDR_adunit[i] = State.cumCT_adunit[i]=State.cumV_adunit[i]=State.cumVG_adunit[i]=State.cumC_adunit[i]=State.cumCC_adunit[i]=0; //added hospitalisation, added detected cases, contact tracing per adunit, cases who are contacts: ggilani 03/02/15, 15/06/17
			State.ETU_adunit[i] = State.H_adunit[i] = State.NumBeds_adunits[i]=0;
			AdUnits[i].place_close_trig=0;
			AdUnits[i].revacc = 0;
			AdUnits[i].currentETUBeds=0; //reset occupied beds to zero;
			AdUnits[i].ETUbedsActive = 0;
			AdUnits[i].totalETUBeds = 0;
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
			AdUnits[i].contactTraceThresholdCrossed=0;
			AdUnits[i].nct=0; //no-one being contact traced at beginning of run
			AdUnits[i].nct_queue=0; //no-one in contact tracing queue at beginning of run
			AdUnits[i].nh_queue=0; //no-one in hospital queue at beginning of run
			AdUnits[i].contactTraceStartDay=1e6;
			AdUnits[i].contactTraceCapacity = P.AdunitCTCapacity;
			AdUnits[i].nextTimeToSDB = 0;
			AdUnits[i].maxSDB = P.AdunitSDBCapacity;
			}
	for(j=0;j<MAX_NUM_THREADS;j++)
		{
		StateT[j].L=StateT[j].I=StateT[j].R=0;
		StateT[j].cumI = StateT[j].cumR = StateT[j].cumC = StateT[j].cumFC = StateT[j].cumFI = StateT[j].cumETU = StateT[j].cumH=StateT[j].cumCT = StateT[j].cumCC = StateT[j].cumTC = StateT[j].cumD = StateT[j].cumDC = StateT[j].cumDD = StateT[j].cumSDB=0; //added setting funeral infections cumFI to zero: ggilani 24/10/14
		StateT[j].cumInf_h=StateT[j].cumInf_n=StateT[j].cumInf_s=StateT[j].cumHQ=StateT[j].cumAC=StateT[j].cumACS=StateT[j].cumAH=StateT[j].cumAA=StateT[j].cumAPC=StateT[j].cumAPA=StateT[j].cumAPCS=0;
		StateT[j].cumT=StateT[j].cumUT=StateT[j].cumTP=StateT[j].cumV=StateT[j].sumRad2=StateT[j].maxRad2=StateT[j].cumV_daily=0;
		StateT[j].nringvacc_queue = StateT[j].nvacc_queue = StateT[j].vacc_cum = StateT[j].ngeovacc_queue = StateT[j].ringvacc_cum = StateT[j].geovacc_cum = 0; //reset ring and geo vaccination variables, 21/08/19
		for(i=0;i<NUM_AGE_GROUPS;i++) StateT[j].cumCa[i]=StateT[j].cumIa[i]=StateT[j].cumDa[i]= StateT[j].cumDCa[i] = StateT[j].cumETUa[i] = StateT[j].cumHa[i]= StateT[j].cumVa[i] = 0; //adding det cases, hosp, vacc by age: ggilani 22/02/22
		for(i=0;i<P.EvolResistNumTypes;i++) StateT[j].cumC_resist[i]=StateT[j].cumI_resist[i]=StateT[j].cumT_resist[i]=0;
		for(i=0;i<2;i++) StateT[j].cumC_keyworker[i]=StateT[j].cumI_keyworker[i]=StateT[j].cumT_keyworker[i]=0;
		for(i=0;i<NUM_PLACE_TYPES;i++) StateT[j].NumPlacesClosed[i]=0;
		for(i=0;i<INFECT_TYPE_MASK;i++) StateT[j].cumItype[i]=0;
		//initialise cumulative case counts per country per thread to zero: ggilani 12/11/14
		for(i=0;i<MAX_COUNTRIES;i++) StateT[j].cumC_country[i]=0;
		if(P.DoAdUnits)
			for(i=0;i<=P.NumAdunits;i++) 
				StateT[j].cumI_adunit[i]=StateT[j].cumC_adunit[i]=StateT[j].cumT_adunit[i]=StateT[j].cumETU_adunit[i]=StateT[j].cumH_adunit[i]=StateT[j].ETU_adunit[i]= StateT[j].H_adunit[i] = StateT[j].cumDC_adunit[i]= StateT[j].cumD_adunit[i]=StateT[j].cumDD_adunit[i]= StateT[j].cumSDB_adunit[i]=StateT[j].cumDR_adunit[i]=StateT[j].cumCT_adunit[i]=StateT[j].cumV_adunit[i]=StateT[j].cumV_adunit[i]=StateT[j].cumC_adunit[i]=StateT[j].cumCC_adunit[i]= StateT[j].nct_queue[i] = 0; //added hospitalisation, detected cases, contact tracing per adunit, cases who are contacts: ggilani 03/02/15, 15/06/17
		}
	nim=0;

	//added this to reorder the susceptible array to be in the same order as the member array - ggilani 12/03/17
	//put this back in to try and fix a bug - it seems to be working but I don't think I should need this - ggilani  22/08/19
	//for (i = 0; i < P.NC; i++)
	//{
	//	for (j = 0; j < Cells[i].n; j++) Cells[i].susceptible[j] = Cells[i].members[j];
	//}


#pragma omp parallel for private(i,b,j,k,l,m,tn,stt,stp) reduction(+:nim) schedule(static,1)
	for(tn=0;tn<P.NumThreads;tn++)
		{
		stp=P.NC/P.NumThreads+1;
		stt=(tn+1)*stp;
		if(stt>P.NC) stt=P.NC;
		for(i=tn*stp;i<stt;i++)
			{
			if((Cells[i].tot_treat!=0)||(Cells[i].tot_vacc!=0)||(Cells[i].S!=Cells[i].n)||(Cells[i].D>0)||(Cells[i].R>0))
				{
				Cells[i].S=Cells[i].n;
				Cells[i].L=Cells[i].I=Cells[i].R=Cells[i].cumTC=Cells[i].D=0;
				Cells[i].infected=Cells[i].latent=Cells[i].susceptible+Cells[i].S;
				Cells[i].tot_treat=Cells[i].tot_vacc=0;

				for(l=0;l<MAX_INTERVENTION_TYPES;l++) Cells[i].CurInterv[l]=-1;
				for(j=0;j<Cells[i].n;j++)
					{
					k=Cells[i].members[j];
					Cells[i].susceptible[j] = k; //added this in here instead
					if(Households[Hosts[k].hh].stockpile!=0) Households[Hosts[k].hh].stockpile=1;
					if(P.DoAirports) Hosts[k].PlaceLinks[HOTEL_PLACE_TYPE]=-1;
					Hosts[k].quar_start_time=
						Hosts[k].isolation_start_time=Hosts[k].absent_start_time=USHRT_MAX-1;
					Hosts[k].absent_stop_time=0;
					Hosts[k].quar_comply=2;
					if((AdUnits[Mcells[Hosts[k].mcell].adunit].id/P.CountryDivisor)==P.TargetCountry)
					{
						Hosts[k].susc=P.RelativeSusceptibilityGuinea;
					}
					else if ((AdUnits[Mcells[Hosts[k].mcell].adunit].id/P.CountryDivisor)==P.TargetCountry2)
					{	
						Hosts[k].susc=P.RelativeSusceptibilityLiberia;
					}
					else
					{
						Hosts[k].susc=1.0;
					}
					Hosts[k].to_die=0;
					Hosts[k].Travelling=0;
					Hosts[k].detected=0; //set detected to zero initially: ggilani - 19/02/15
					Hosts[k].dayDetected=USHRT_MAX-1; //set day on which each case is detected to zero initially: ggilani - 23/06/15
					Hosts[k].inf=0;
					Hosts[k].listpos=j;
					Hosts[k].treat_stop_time=Hosts[k].num_treats=Hosts[k].contactTraced_end_time=0;
					Hosts[k].vacc_start_time=Hosts[k].treat_start_time=Hosts[k].contactTraced_start_time=USHRT_MAX-1;
					Hosts[k].revacc = 0;
					Hosts[k].latent_time=Hosts[k].recovery_time=Hosts[k].hospital_time=Hosts[k].detect_time=0; //also set hospitalisation time to zero: ggilani 28/10/2014
					Hosts[k].hospitalised=Hosts[k].etu=0; //set hospitalised flag to zero: ggilani 28/10/14
					Hosts[k].contactTraced=0;
					Hosts[k].vaccRing = 0; //set flag for vacc ring to zero: ggilani 29/05/19
					Hosts[k].ringCase = -1; //set flag for ring case to -1;
					Hosts[k].resist=0;
					Hosts[k].infector=-1;
					Hosts[k].infect_type=0;
					Hosts[k].infectiousMult = 1; //reset to 1 - this is changed when funeral transmission temporarily increases infectiousness
#ifdef FRESSCA
					Hosts[k].vacc_queue_pos=-1;
#endif
					}
	// Next loop needs to count down for DoImmune host list reordering to work
				for(j=Cells[i].n-1;j>=0;j--)
					{
					k=Cells[i].members[j];
#ifdef NEW_AGE_MODEL
#ifdef FRESSCA
					if(P.DoDistributionVaccination) UpdateVaccStatus(k,tn);
#endif
#endif
					if (P.DoWholeHouseholdImmunity)
					{
						if (P.InitialImmunity[0] != 0)
						{
							if (Households[Hosts[k].hh].FirstPerson == k)
							{
								for (m = 0; m < Households[Hosts[k].hh].nh; m++)
									Hosts[k + m].inf = 0;
								if ((P.InitialImmunity[0] == 1) || (ranf_mt(ranf_seed, tn) < P.InitialImmunity[0]))
									for (m = Households[Hosts[k].hh].nh - 1; m >= 0; m--)
										DoImmune(k + m);
							}
						}
					}
					Hosts[k].inf=0;
					if(P.DoInitEquilib)
						{
						if(P.DoSIS)
							{
							if(P.SuscReductionFactorPerInfection>0)
								Hosts[k].susc=exp(-P.MeanAnnualDeathRate*((double) HOST_AGE_YEAR(k)));
							else
								Hosts[k].susc=(ranf_mt(ranf_seed, tn)>exp(-P.MeanAnnualDeathRate*((double) HOST_AGE_YEAR(k))))?0:1;
							}
						else if (ranf_mt(ranf_seed, tn)>exp(-P.MeanAnnualDeathRate*((double) HOST_AGE_YEAR(k))))
							DoImmune(k);
						}
					else
						{
						m=HOST_AGE_GROUP(k);
						if(P.DoSIS)
							{
							if(P.SuscReductionFactorPerInfection>0)
								Hosts[k].susc=1-P.InitialImmunity[m];
							else
								nim+=(Hosts[k].susc=(ranf_mt(ranf_seed,tn)<P.InitialImmunity[m])?0:1);
							}
						else
							{
							if((P.InitialImmunity[m]==1)||((P.InitialImmunity[m]>0)&&(ranf_mt(ranf_seed, tn)<P.InitialImmunity[m]))) {DoImmune(k);nim+=1;}
							}
						}
					}
				}
			}
	   }

	fprintf(stderr,"Finished cell init - %i people assigned as immune.\n",nim);


#pragma omp parallel for private(i,j,k,i2,j2,l) schedule(static,500)
	for(l=0;l<P.NMCP;l++)
		{
		i=(int) (McellLookup[l]-Mcells);
		Mcells[i].vacc_start_time=Mcells[i].treat_start_time=USHRT_MAX-1;
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
		Mcells[i].treat_end_time=0;
		Mcells[i].treat_trig=Mcells[i].vacc_trig=Mcells[i].vacc=Mcells[i].treat=0;
		Mcells[i].place_trig=Mcells[i].move_trig=Mcells[i].socdist_trig=Mcells[i].keyworkerproph_trig=
			Mcells[i].placeclose=Mcells[i].moverest=Mcells[i].socdist=Mcells[i].keyworkerproph=0;
		Mcells[i].move_start_time=USHRT_MAX-1;
		Mcells[i].place_end_time=Mcells[i].move_end_time=
			Mcells[i].socdist_end_time=Mcells[i].keyworkerproph_end_time=0;
		if(P.DoPlaces)
			for(j=0;j<P.PlaceTypeNum;j++)
				for(k=0;k<Mcells[i].np[j];k++)
					{
					j2=Mcells[i].places[j][k];
					Places[j][j2].treat=Places[j][j2].control_trig=0;
					Places[j][j2].treat_end_time=Places[j][j2].close_end_time=0;
					Places[j][j2].close_start_time=USHRT_MAX-1;
#ifdef ABSENTEEISM_PLACE_CLOSURE
					Places[j][j2].AbsentLastUpdateTime=0;
					for(i2=0;i2<MAX_ABSENT_TIME;i2++) Places[j][j2].Absent[i2]=0;
#endif
					}
		}


	for(i=0;i<MAX_NUM_THREADS;i++)
		{
		for(j=0;j<MAX_NUM_THREADS;j++)
			StateT[i].n_queue[j]=0;
		for(j=0;j<P.PlaceTypeNum;j++)
			StateT[i].np_queue[j]=0;
		}
	if(DoInitUpdateProbs) 
		{
		UpdateProbs(0);
		DoInitUpdateProbs=0;
		}
	//initialise event log to zero at the beginning of every run: ggilani - 10/10/2014. UPDATE: 15/10/14 - we are now going to store all events from all realisations in one file
	if((P.DoRecordInfEvents)&&(P.RecordInfEventsPerRun))
	{
		*nEvents=0;
		for(i=0;i<P.MaxInfEvents;i++)
		{
			InfEventLog[i].t=InfEventLog[i].infectee_x=InfEventLog[i].infectee_y=InfEventLog[i].t_infector=0.0;
			InfEventLog[i].infectee_ind=InfEventLog[i].infector_ind=0;
			InfEventLog[i].infectee_adunit = InfEventLog[i].listpos = InfEventLog[i].infectee_cell = InfEventLog[i].infector_cell = InfEventLog[i].thread = 0;
		}
	}

	SeedInfection(0,P.NumInitialInfections,0,run);
	P.ControlPropCasesId = P.PostAlertControlPropCasesId;
	//P.ControlPropCasesId=P.PreAlertControlPropCasesId;
	P.TreatTimeStart=1e10;

//  Mass Vacc now starts after outbreak alert trigger
//	if(P.DoMassVacc)
//		P.VaccTimeStart=P.VaccTimeStartBase;
//	else
	P.VaccTimeStart = 1e10;
	P.MoveRestrTimeStart=1e10;
	P.PlaceCloseTimeStart=1e10;
	P.PlaceCloseTimeStart2=1e10;
	P.SocDistTimeStart=1e10;
	P.AirportCloseTimeStart=1e10;
	P.CaseIsolationTimeStart=1e10;
	P.HQuarantineTimeStart=1e10;
	P.KeyWorkerProphTimeStart=1e10;
	//P.GeoVaccTimeStart = 1e10; //reset additional start times - ggilani 13/09/23
	//P.RingVaccTimeStart = 1e10;
	P.FuneralControlTimeStart = 1e10;
	P.ContactTracingTimeStart = 1e10;
	P.ETUTimeStart = 1e10;
	P.TreatMaxCourses=P.TreatMaxCoursesBase;
	P.VaccMaxCourses=P.VaccMaxCoursesBase;
	P.PlaceCloseDuration=P.PlaceCloseDurationBase;
	P.VaccNewCoursesStartTime = 1e10;
	//if do hospitalisation, reset a couple of hospitalisation parameters
	P.CurrIndETUBeds=0;
	P.CurrIndMeanTimeToHosp=0;
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

	fprintf(stderr,"Finished InitModel.\n");
}

void SeedInfection(double t,int *nsi,int rf, int run) //adding run number to pass it to event log
{
   int i,j,k,l,m,f,f2,n;

   f2=((t>=0)&&(P.DoImportsViaAirports));
   n=((rf==0)?P.NumSeedLocations:1);
   for(i=0;i<n;i++)
	{
	if((!P.DoRandomInitialInfectionLoc)||((P.DoAllInitialInfectioninSameLoc)&&(rf)))
		{
			//if we are importing cases to a specific location, and this is not initial seeding
			if((rf==1)&&(P.DoImportToSpecLocation))
			{
				k=(int) (P.ImportLocation[0]/P.mcwidth);
				l=(int) (P.ImportLocation[1]/P.mcheight);
			}
			else
			{
				k=(int) (P.LocationInitialInfection[i][0]/P.mcwidth);
				l=(int) (P.LocationInitialInfection[i][1]/P.mcheight);
			}
		j=k*P.nmch+l;
		m=0;
		for(k=0;(k<nsi[i])&&(m<10000);k++)
			{
			l=Mcells[j].members[(int) (ranf()*((double) Mcells[j].n))];
			if(Hosts[l].inf==0)
				{
				if(CalcPersonSusc(l,0,0,0)>0)
					{
					//only reset the initial location if rf==0, i.e. when initial seeds are being set, not when imported cases are being set
						if(rf==0)
						{
							P.LocationInitialInfection[i][0]=Households[Hosts[l].hh].loc_x;
							P.LocationInitialInfection[i][1]=Households[Hosts[l].hh].loc_y;
						}
					Hosts[l].infector=-2;Hosts[l].infect_type=INFECT_TYPE_MASK-1;
					Hosts[l].base_inf_level=1.0;
					Hosts[l].resist=CalcSeedResist();
					DoInfect(l,t,0,run);
					m=0;
					}
				}
			else
				{k--;m++;}
			}
		}
	else if(P.DoAllInitialInfectioninSameLoc)
		{
		f=0;
		do
			{
			m=0;
			do
				{
				l=(int) (ranf()*((double) P.N));
				j=Hosts[l].mcell;
				//fprintf(stderr,"%i ",AdUnits[Mcells[j].adunit].id);
				}
			while((Mcells[j].n<nsi[i])||(Mcells[j].n>P.MaxPopDensForInitialInfection)
				||(Mcells[j].n<P.MinPopDensForInitialInfection)||(((int)(AdUnits[Mcells[j].adunit].id/P.CountryDivisor)!=P.TargetCountry)&&(P.TargetCountry>=0)) // from Mcells[j].country to Mcells[j].adunit/P.CountryDivisor
				||((P.InitialInfectionsAdminUnit[i]>0)&&((AdUnits[Mcells[j].adunit].id%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor!=(P.InitialInfectionsAdminUnit[i]%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor)));
			for(k=0;(k<nsi[i])&&(m<10000);k++)
				{
				l=Mcells[j].members[(int) (ranf()*((double) Mcells[j].n))];
				if(Hosts[l].inf==0)
					{
					if(CalcPersonSusc(l,0,0,0)>0)
						{
						P.LocationInitialInfection[i][0]=Households[Hosts[l].hh].loc_x;
						P.LocationInitialInfection[i][1]=Households[Hosts[l].hh].loc_y;
						Hosts[l].infector=-2;Hosts[l].infect_type=INFECT_TYPE_MASK-1;
						Hosts[l].base_inf_level=1.0;
						Hosts[l].resist=CalcSeedResist();
						DoInfect(l,t,0,run);
						m=0;
						}
					}
				else
					{k--;m++;}
				}
			if(m)
				f++;
			else
				f=0;
			}
		while((f>0)&&(f<1000));
		}
	else
		{
		m=0;
		for(k=0;(k<nsi[i])&&(m<10000);k++)
			{
			do
				{
				l=(int) (ranf()*((double) P.N));
				j=Hosts[l].mcell;
				//fprintf(stderr,"%i ",AdUnits[Mcells[j].adunit].id);
				}
			while((Mcells[j].n==0)||(Mcells[j].n>P.MaxPopDensForInitialInfection)
				||(Mcells[j].n<P.MinPopDensForInitialInfection)||(((int)(AdUnits[Mcells[j].adunit].id/P.CountryDivisor)!=P.TargetCountry)&&(P.TargetCountry>=0)) // from Mcells[j].country to Mcells[j].adunit/P.CountryDivisor
				||((P.InitialInfectionsAdminUnit[i]>0)&&((AdUnits[Mcells[j].adunit].id%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor!=(P.InitialInfectionsAdminUnit[i]%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor)));
			l=Mcells[j].members[(int) (ranf()*((double) Mcells[j].n))];
			if(Hosts[l].inf==0)
				{
				if(CalcPersonSusc(l,0,0,0)>0)
					{
					P.LocationInitialInfection[i][0]=Households[Hosts[l].hh].loc_x;
					P.LocationInitialInfection[i][1]=Households[Hosts[l].hh].loc_y;
					Hosts[l].infector=-2;Hosts[l].infect_type=INFECT_TYPE_MASK-1;
					Hosts[l].base_inf_level=1.0;
					Hosts[l].resist=CalcSeedResist();
					DoInfect(l,t,0,run);
					m=0;
					}
				else
					{k--;m++;}
				}
			else
				{k--;m++;}
			}
		}
   }
   if(m>0) fprintf(stderr,"### Seeding error ###\n");
}

int CalcSeedResist(void)
{
	int i,j;
	double t;

	if(MAX_NUM_RESIST_TYPES==1)
		i=0;
	else
		{
		t=ranf();
		for(i=0;(i<MAX_NUM_RESIST_TYPES)&&(t>=0);i++)
			{
			t-=P.EvolResistSeedProp[i];
			}
		i--;
		}
	return i;
}

void RunModel(int run) //added run number as parameter
{
	int i,j,k,l,fs,fs2,ni,i2,nu,in;
	double ir,t,cI,lcI,t2;
	unsigned short int ts;
	int continueEvents=1;

	lcI=1;
	if(P.DoLoadSnapshot)
		{
		P.ts_age=(int) (P.SnapshotLoadTime*P.TimeStepsPerDay);
		t=((double) P.ts_age)*P.TimeStep;
		}
	else
		{
		t=0;
		P.ts_age=0;
		}
	fs=1;
	fs2=0;
	nu=0;
#ifdef FRESSCA
    //if(P.DoDistributionVaccination) VaccineDistributionNetwork->Update(30);
	//VaccineDistributionNetwork->ResetTimeSteps();  //Needed to sync it all up.
#endif

	for(ns=1;((ns<P.NumSamples)&&(!InterruptRun));ns++) //&&(continueEvents) <-removed this
		{
		if (continueEvents)
		{
			RecordSample(t, ns - 1);
			//update hospitalisation parameters at the beginning of every time step? ggilani - 11/03/2017
			if ((P.DoHospitalisation)&&(t>=P.ETUTimeStart))
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
			if ((P.DoRingVaccination)&&(t>P.RingVaccTimeStart))
			{
				UpdateVaccination(t,ns-1);
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
		if(P.LimitNumInfections) continueEvents=(State.cumI<P.MaxNumInfections);
		fprintf(stderr,"\r    t=%lg   %i    %i|%i    %i     %i   %i (%lg %lg %lg)   %lg    ",t,State.S,State.L,State.I,State.R,State.D,State.cumD,State.cumT,State.cumV,State.cumVG,sqrt(State.maxRad2)/1000); //added State.cumVG
		for(j=0;((j<P.UpdatesPerSample)&&(!InterruptRun)&&(continueEvents));j++)
			{
			ts=(unsigned short int) (P.TimeStepsPerDay*t);

			//check to see if interruptions to interventions are in place
			if (P.DoInterruptIntervention)
			{
				P.InterruptIntervention = 0;
				for (in = 0; in < P.NDaysInterrupt; in++)
				{
					if (P.DaysInterruptIntervention[in]== (int) t)
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
				if ((P.ResetSeedsFlag == 0) && (ts >= (P.TimeToResetSeeds*P.TimeStepsPerDay)))
				{
					P.newseed1 = (int)(ranf() * 1e8);
					P.newseed2 = (int)(ranf() * 1e8);
					setall(P.newseed1, P.newseed2);
					P.ResetSeedsFlag = 1;
				}
			}

			if(fs)
				{
				if(P.DoAirports) TravelDepartSweep(t);
				k=(int) t;
				if(P.DurImportTimeProfile>0)
					{
					if(k<P.DurImportTimeProfile)
						ir=P.ImportInfectionTimeProfile[k]*((t>P.InfectionImportChangeTime)?(P.InfectionImportRate2/P.InfectionImportRate1):1.0);
					else
						ir=0;
					}
				else
					ir=(t>P.InfectionImportChangeTime)?P.InfectionImportRate2:P.InfectionImportRate1;
#ifdef IMPORT_POP_SIZE
				ir*=((double) P.N)/IMPORT_POP_SIZE;
#endif
				if(ir>0)
					{
					ni=(int) ignpoi(P.TimeStep*ir);
					if(ni>0)
						{
						SeedInfection(t,&ni,1,run);
						}
					}
				if(P.FalsePositivePerCapitaIncidence>0)
					{
					ni=(int) ignpoi(P.TimeStep*P.FalsePositivePerCapitaIncidence*((double) P.N));
					if(ni>0)
						{
							for(k=0;k<ni;k++)
								{
								do
									{l=(int) (((double) P.N)*ranf());}
								while((abs(Hosts[l].inf)==5)||(ranf()>P.FalsePositiveAgeRate[HOST_AGE_GROUP(l)]));
								DoFalseCase(l,t,ts,0);
								}
						}
					}
				InfectSweep(t, run); //adding run number as a parameter to infect sweep so we can track run number: ggilani - 15/10/14
				if(!P.DoSI) IncubRecoverySweep(t,run);
				// If we are doing hospitalisation by admin unit, process the hospitalisation queues filled during the incubation recovery sweep: ggilani - 24/11/14
				if((P.DoHospitalisation)&&(P.DoETUByAdUnit))
				{
					HospitalSweepAdunits(t);
				}
				// If doing new contact tracing, update numbers of people under contact tracing after each time step
				if((P.DoContactTracing)&&(t>=P.ContactTracingTimeStart))
				{
					ContactTracingSweep(t);
				}
				if (((P.DoRingVaccination)||(P.DoGeoVaccination))&&(t >= P.VaccTimeStart))
				{
					VaccSweep(t); //process ring vaccination queue
				}
				

#ifdef NEW_AGE_MODEL
				if((P.DoDeath)&&(!P.DoSIS)&&(nu%P.UpdatesPerDemogUpdate==0)) DemogSweep(t);
#endif

#ifdef FRESSCA
    			if(P.DoDistributionVaccination)
					{
					if((ts==P.usSIAStartTime)||((ts>P.usSIAStartTime)&&((ts-P.usSIAStartTime)%P.usSIARepeatInterval==0)))
						SIASweep(ts);
					if(nu%P.UpdatesPerDemogUpdate==0)  //will only work with daily timestep
        				{
						DistrVaccSweep(t);
						VaccineDistributionNetwork->PrintGEFrame(KMLFile,P.UpdatesPerDemogUpdate);
						VaccineDistributionNetwork->Update(P.UpdatesPerDemogUpdate);
        				}
					}
#endif
				nu++;
				fs2=((P.DoDeath)||(P.DoSIS)||(State.L+State.I>0)||(ir>0)||(P.FalsePositivePerCapitaIncidence>0));
				if(!TreatSweep(t))
					{
					if((!fs2)&&(State.L+State.I==0)&&(P.FalsePositivePerCapitaIncidence==0)) {if((ir==0)&&(((int) t)>P.DurImportTimeProfile)) fs=0;}
					}
				if(P.DoAirports) TravelReturnSweep(t);
				}
			t+=P.TimeStep;
			if(P.DoDeath) P.ts_age++;
			if((P.DoSaveSnapshot)&&(t<=P.SnapshotSaveTime)&&(t+P.TimeStep>P.SnapshotSaveTime)) SaveSnapshot();
			if(t>P.TreatNewCoursesStartTime) P.TreatMaxCourses+=P.TimeStep*P.TreatNewCoursesRate;
			if ((t > P.VaccNewCoursesStartTime) && (t < P.VaccNewCoursesEndTime))
			{
				if (P.DoVaccDailyReplenishment)
				{
					P.VaccMaxCourses += P.TimeStep * P.VaccNewCoursesRate;
				}
				else if (P.DoVaccBulkReplenishment)
				{
					P.VaccMaxCourses += P.VaccNewCoursesBulk;
					P.VaccNewCoursesStartTime += P.VaccNewCoursesDelay;
				}
			}
			cI=((double) (State.S))/((double) P.N);
#ifndef NEW_AGE_MODEL			
			if(((lcI-cI)>0.2)&&(!P.DoSIS)) 
				{
				lcI=cI;
				UpdateProbs(0);
				DoInitUpdateProbs=1;
				}
#endif
			}

		}
	RecordSample(t,P.NumSamples-1);
	fprintf(stderr,"\nEnd of run\n");
#ifdef NEW_AGE_MODEL
	if((P.DoDeath)&&(!P.DoSIS))
		DemogSweep(t);
	else
		for(i2=0;i2<P.N;i2++) UpdatePersonDemog(i2,0);
#endif
	t2=t+P.SampleTime;
	while(fs)
		{
		fs=TreatSweep(t2);
		t2+=P.SampleStep;
		}
//	fprintf(stderr,"End RunModel\n");
	if(P.DoAirports)
		{
		t2=t;
		for(t2=t;t2<=t+MAX_TRAVEL_TIME;t2+=P.TimeStep)
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

void SaveAgeDistrib(void)
	{
	int i;
	FILE *dat;
	char outname[1024];

	sprintf(outname,"%s.agedist.xls",OutFile);
	if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
	if(P.DoDeath)
		{
		fprintf(dat,"age\tfreq\tlifeexpect\n");
		for(i=0;i<NUM_AGE_GROUPS;i++)
			fprintf(dat,"%i\t%lg\t%lg\n",i,AgeDist[i],AgeDist2[i]);
		fprintf(dat,"\np\tlife_expec\tage\n");
		for(i=0;i<=1000;i++)
			fprintf(dat,"%lg\t%lg\t%i\n",((double) i)/1000,P.InvLifeExpecDist[0][i],State.InvAgeDist[0][i]);
		}
	else
		{
		fprintf(dat,"age\tfreq\n");
		for(i=0;i<NUM_AGE_GROUPS;i++)
			fprintf(dat,"%i\t%lg\n",i,AgeDist[i]);
		}

	fclose(dat);
}

void SaveAgeDistrib2(void)
	{
#ifdef NEW_AGE_MODEL
	int i,j;
	FILE *dat;
	char outname[1024];

	for(i=0;i<NUM_AGE_GROUPS;i++) AgeDist[i]=AgeDist2[i]=0;
	for(i=0;i<P.N;i++)
		{
		j=HOST_AGE_GROUP(i);
		if(j>=NUM_AGE_GROUPS) j=NUM_AGE_GROUPS-1;
		AgeDist[j]++;
		if(P.DoDeath) AgeDist2[Hosts[i].life_expectancy/((int) P.TimeStepsPerYear)/AGE_GROUP_WIDTH]++;
		}

	sprintf(outname,"%s.agedist2.xls",OutFile);
	if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
	if(P.DoDeath)
		{
		fprintf(dat,"age\tfreq\tlifeexpect\n");
		for(i=0;i<NUM_AGE_GROUPS;i++)
			fprintf(dat,"%i\t%lg\t%lg\n",i,AgeDist[i],AgeDist2[i]);
		}
	else
		{
		fprintf(dat,"age\tfreq\n");
		for(i=0;i<NUM_AGE_GROUPS;i++)
			fprintf(dat,"%i\t%lg\n",i,AgeDist[i]);
		}

	fclose(dat);
#endif
}

void SaveDistribs(void)
{
	int i,j,k;
	FILE *dat;
	char outname[1024];
	double s;

	if(P.DoPlaces)
		{
		for(j=0;j<P.PlaceTypeNum;j++)
			if(j!=HOTEL_PLACE_TYPE)
				{
				for(i=0;i<P.Nplace[j];i++)
					Places[j][i].n=0;
				for(i=0;i<P.N;i++)
					{
					if(Hosts[i].PlaceLinks[j]>=P.Nplace[j])
						fprintf(stderr,"*%i %i: %i %i",i,j,Hosts[i].PlaceLinks[j],P.Nplace[j]);
					else if(Hosts[i].PlaceLinks[j]>=0)
						Places[j][Hosts[i].PlaceLinks[j]].n++;
					}
				}
		for(j=0;j<P.PlaceTypeNum;j++)
			for(i=0;i<MAX_DIST;i++)
				PlaceDistDistrib[j][i]=0;
		for(i=0;i<P.N;i++)
			for(j=0;j<P.PlaceTypeNum;j++)
				if((j!=HOTEL_PLACE_TYPE)&&(Hosts[i].PlaceLinks[j]>=0))
					if(Hosts[i].PlaceLinks[j]>=P.Nplace[j])
						fprintf(stderr,"*%i %i: %i ",i,j,Hosts[i].PlaceLinks[j]);
					else if((!P.DoOutputPlaceDistForOneAdunit)||
						((AdUnits[Mcells[Hosts[i].mcell].adunit].id%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor==(P.OutputPlaceDistAdunit%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor))
						{
						k=Hosts[i].PlaceLinks[j];
						s=sqrt(dist2_raw(Households[Hosts[i].hh].loc_x,Households[Hosts[i].hh].loc_y,Places[j][k].loc_x,Places[j][k].loc_y))/OUTPUT_DIST_SCALE;
						k=(int) s;
						if(k<MAX_DIST) PlaceDistDistrib[j][k]++;
						}
		for(j=0;j<P.PlaceTypeNum;j++)
			for(i=0;i<MAX_PLACE_SIZE;i++)
				PlaceSizeDistrib[j][i]=0;
		for(j=0;j<P.PlaceTypeNum;j++)
			if(j!=HOTEL_PLACE_TYPE)
				for(i=0;i<P.Nplace[j];i++)
					if(Places[j][i].n<MAX_PLACE_SIZE)
						PlaceSizeDistrib[j][Places[j][i].n]++;
		sprintf(outname,"%s.placedist.xls",OutFile);
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat,"dist");
		for(j=0;j<P.PlaceTypeNum;j++)
			if(j!=HOTEL_PLACE_TYPE)
				fprintf(dat,"\tfreq_p%i",j);
		fprintf(dat,"\n");
		for(i=0;i<MAX_DIST;i++)
			{
			fprintf(dat,"%i",i);
			for(j=0;j<P.PlaceTypeNum;j++)
				if(j!=HOTEL_PLACE_TYPE)
					fprintf(dat,"\t%i",PlaceDistDistrib[j][i]);
			fprintf(dat,"\n");
			}
		fclose(dat);
		sprintf(outname,"%s.placesize.xls",OutFile);
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat,"size");
		for(j=0;j<P.PlaceTypeNum;j++)
			if(j!=HOTEL_PLACE_TYPE)
				fprintf(dat,"\tfreq_p%i",j);
		fprintf(dat,"\n");
		for(i=0;i<MAX_PLACE_SIZE;i++)
			{
			fprintf(dat,"%i",i);
			for(j=0;j<P.PlaceTypeNum;j++)
				if(j!=HOTEL_PLACE_TYPE)
					fprintf(dat,"\t%i",PlaceSizeDistrib[j][i]);
			fprintf(dat,"\n");
			}
		fclose(dat);
		}
	if((P.DoAdunitBoundaries)&&(P.PrivateTreatPropCases>0)&&(P.DoHouseholds))
		{
		sprintf(outname,"%s.privatestock.xls",OutFile);
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat,"adunit_code\thouseholds\thouseholds_with_stock\n");
		for(j=0;j<P.NumAdunits;j++)
			fprintf(dat,"%i\t%i\t%i\n",AdUnits[j].id,P.HouseholdsByAdunit[j],P.PrivateStockByAdunit[j]);
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
	int i,j;
	FILE *dat;
	char outname[1024];

	sprintf(outname,"%s.origdestmat.csv",OutFile);
	if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat,"0,");
	for(i=0;i<P.NumAdunits;i++) fprintf(dat,"%i,",(AdUnits[i].id%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor);
	fprintf(dat,"\n");
	for(i=0;i<P.NumAdunits;i++)
	{
		fprintf(dat,"%i,",(AdUnits[i].id%P.AdunitLevel1Mask)/P.AdunitLevel1Divisor);
		for(j=0;j<P.NumAdunits;j++)
		{
			fprintf(dat,"%lg,",AdUnits[i].origin_dest[j]);
		}
		fprintf(dat,"\n");
	}
	fclose(dat);
}

void SaveResults(void)
{
	int i,j,k;
	double t;
	FILE *dat;
	char outname[1024];

	sprintf(outname,"%s.csv",OutFile);
	if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat,"t,S,L,I,R,D,incI,incR,incFC,incFI,incC,incDC,incD,incDD,incSDB,incTC,incETU,incH,incCT,incCC,cumT,cumTP,cumV,capV,cumVG,capVG,nBeds,Extinct,Detected,rmsRad,maxRad\n");//\t\t%lg\t%lg\t%lg\n",P.R0household,P.R0places,P.R0spatial);
	for(i=0;i<P.NumSamples;i++)
		{
		fprintf(dat,"%lg,%lf,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
			TimeSeries[i].t,TimeSeries[i].S,TimeSeries[i].L,TimeSeries[i].I,
			TimeSeries[i].R,TimeSeries[i].D,TimeSeries[i].incI,
			TimeSeries[i].incR,TimeSeries[i].incFC,TimeSeries[i].incFI,TimeSeries[i].incC,TimeSeries[i].incDC, TimeSeries[i].incD, TimeSeries[i].incDD, TimeSeries[i].incSDB,TimeSeries[i].incTC,TimeSeries[i].incETU,TimeSeries[i].incH,TimeSeries[i].incCT,TimeSeries[i].incCC, //added incidence funeral transmissions and hospitalisation
			TimeSeries[i].cumT,TimeSeries[i].cumTP,TimeSeries[i].cumV, TimeSeries[i].capV, TimeSeries[i].cumVG, TimeSeries[i].capVG, TimeSeries[i].nBeds,TimeSeries[i].extinct,TimeSeries[i].detected,TimeSeries[i].rmsRad,TimeSeries[i].maxRad);
		}
	fclose(dat);
	if (P.DoControlOutput)
	{
		sprintf(outname, "%s.controls.xls", OutFile);
		if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t\tS\tincC\tincTC\tincFC\tincFI\tincETU\tincH\cumT\tcumUT\tcumTP\tcumV\tincHQ\tincAC\tincAH\tincAA\tincACS\tincAPC\tincAPA\tincAPCS");
		for (j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\tprClosed_%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lg\t%lf\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg",
				TimeSeries[i].t, TimeSeries[i].S, TimeSeries[i].incC, TimeSeries[i].incTC, TimeSeries[i].incFC, TimeSeries[i].incFI, TimeSeries[i].incETU,TimeSeries[i].incH, //added incidence of funeral transmissions and hospitalisation
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
	if((P.DoAdUnits)&&(P.DoAdunitOutput))
		{
		sprintf(outname,"%s.adunit.csv",OutFile); //modifying to csv file
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat,"t,");
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,"I_%s,",AdUnits[i].ad_name); //"\tI%i"
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,"C_%s,",AdUnits[i].ad_name); //"\tC%i"
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,"DC_%s,",AdUnits[i].ad_name); //added detected cases: ggilani 03/02/15
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "D_%s,", AdUnits[i].ad_name); //added deaths: ggilani 05/10/23
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "DD_%s,", AdUnits[i].ad_name); //added detected deaths: ggilani 05/10/23
		//for(i=0;i<P.NumAdunits;i++) fprintf(dat,"T_%s,",AdUnits[i].ad_name); //"\tT%i"

		if (P.DoFuneralTransmission)
		{
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "SDB_%s,", AdUnits[i].ad_name); //added safe burials: ggilani 05/10/23
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "capSDB_%s,", AdUnits[i].ad_name); //added safe burials: ggilani 05/10/23
		}

		if((P.DoHospitalisation)&(P.DoETUByAdUnit))
		{
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "Beds_%s,", AdUnits[i].ad_name); //"\tT%i" //added number of beds
			for(i=0;i<P.NumAdunits;i++) fprintf(dat,"incETU_%s,",AdUnits[i].ad_name); //"\tT%i" //added incidence of hospitalisation
			for(i=0;i<P.NumAdunits;i++) fprintf(dat,"ETU_%s,",AdUnits[i].ad_name); //"\tT%i" //added hospitalisation
			if (P.DoOutputETUCapacity)
			{
				for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "capETU_%s,", AdUnits[i].ad_name); //"\tT%i" //added hospital capacity indicator: ggilani 25/04/22
			}
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "incH_%s,", AdUnits[i].ad_name); //"\tT%i" //added incidence of hospitalisation
			//for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "H_%s,", AdUnits[i].ad_name); //"\tT%i" //added hospitalisation
		}
		if(P.DoContactTracing)
		{
			for(i=0;i<P.NumAdunits;i++) fprintf(dat,"incCT_%s,",AdUnits[i].ad_name); //"\tT%i" //added contact tracing
			for(i=0;i<P.NumAdunits;i++) fprintf(dat,"CT_%s,",AdUnits[i].ad_name); //"\tT%i" //added contact tracing
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
		fprintf(dat,"\n");
		for(i=0;i<P.NumSamples;i++)
			{
			fprintf(dat,"%lg,",TimeSeries[i].t); //"%lg"
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,"%lg,",TimeSeries[i].incI_adunit[j]); //"\t%lg"
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,"%lg,",TimeSeries[i].incC_adunit[j]); //"\t%lg"
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,"%lg,",TimeSeries[i].incDC_adunit[j]); //added detected cases: ggilani 03/02/15 
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,"%lg,",TimeSeries[i].incD_adunit[j]); //"\t%lg"
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
			if(P.DoContactTracing)
			{
				for(j=0;j<P.NumAdunits;j++)
					fprintf(dat,"%lg,",TimeSeries[i].incCT_adunit[j]); //"\t%lg" //added contact tracing
				for(j=0;j<P.NumAdunits;j++)
					fprintf(dat,"%lg,",TimeSeries[i].CT_adunit[j]); //"\t%lg" //added contact tracing
				for(j = 0; j < P.NumAdunits; j++)
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
			fprintf(dat,"\n");
			}
		fclose(dat);
		}
	if(P.EvolResistNumTypes>1)
		{
		sprintf(outname,"%s.resist.xls",OutFile);
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat,"t");
		for(i=0;i<P.EvolResistNumTypes;i++) fprintf(dat,"\tI%i",i);
		for(i=0;i<P.EvolResistNumTypes;i++) fprintf(dat,"\tC%i",i);
		for(i=0;i<P.EvolResistNumTypes;i++) fprintf(dat,"\tT%i",i);
		fprintf(dat,"\n");
		for(i=0;i<P.NumSamples;i++)
			{
			fprintf(dat,"%lg",TimeSeries[i].t);
			for(j=0;j<P.EvolResistNumTypes;j++)
				fprintf(dat,"\t%lg",TimeSeries[i].incI_resist[j]);
			for(j=0;j<P.EvolResistNumTypes;j++)
				fprintf(dat,"\t%lg",TimeSeries[i].incC_resist[j]);
			for(j=0;j<P.EvolResistNumTypes;j++)
				fprintf(dat,"\t%lg",TimeSeries[i].cumT_resist[j]);
			fprintf(dat,"\n");
			}
		fclose(dat);
		}
	//if(P.KeyWorkerProphTimeStartBase<P.SampleTime)
	if(P.DoKeyworkerOutput)
		{
		sprintf(outname,"%s.keyworker.csv",OutFile);
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat,"t");
		for(i=0;i<2;i++) fprintf(dat,",I%i",i);
		for(i=0;i<2;i++) fprintf(dat,",C%i",i);
		for(i=0;i<2;i++) fprintf(dat,",T%i",i);
		for (i = 0; i < 2; i++) fprintf(dat, ",D%i", i);
		fprintf(dat,",%i,%i\n",P.KeyWorkerNum,P.KeyWorkerIncHouseNum);
		for(i=0;i<P.NumSamples;i++)
			{
			fprintf(dat,"%lg",TimeSeries[i].t);
			for(j=0;j<2;j++)
				fprintf(dat,",%lg",TimeSeries[i].incI_keyworker[j]);
			for(j=0;j<2;j++)
				fprintf(dat,",%lg",TimeSeries[i].incC_keyworker[j]);
			for(j=0;j<2;j++)
				fprintf(dat,",%lg",TimeSeries[i].cumT_keyworker[j]);
			for (j = 0; j < 2; j++)
				fprintf(dat, ",%lg", TimeSeries[i].incD_keyworker[j]);
			fprintf(dat,"\n");
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
	if(P.DoROutput)
	{ 
	sprintf(outname,"%s.R0.csv",OutFile);
	if (!(dat = fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
	for (i = 0; i < MAX_SEC_REC; i++)
		{
		fprintf(dat,"%i",i);
		for (j = 0; j < MAX_GEN_REC; j++)
			fprintf(dat,",%lg",indivR0[i][j]);
		fprintf(dat,"\n");
		}
	fclose(dat);
	}
	for(i=1;i<=MAX_HOUSEHOLD_SIZE;i++)
		{
		t=0;
		for(j=1;j<=MAX_HOUSEHOLD_SIZE;j++)
			t+=inf_household[i][j];
		inf_household[i][0]=denom_household[i]-t;
		}
	for(i=1;i<=MAX_HOUSEHOLD_SIZE;i++)
		{
		t=0;
		for(j=1;j<=MAX_HOUSEHOLD_SIZE;j++)
			t+=case_household[i][j];
		case_household[i][0]=denom_household[i]-t;
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

	if ((P.DoVaccOutput)&(P.OutbreakDetected))
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
				fprintf(dat, "%i,%i,%i,%i,%i,%lg,%i,%lg,%lg,%i,%lg", j, Mcells[j].n, Mcells[j].ntriggervacc,Mcells[j].totalvacc, Mcells[j].popvacc, Mcells[j].minvaccdist, Mcells[j].minvaccdist_dose, Mcells[j].minvaccdist_t, Mcells[j].maxvaccdist, Mcells[j].maxvaccdist_dose, Mcells[j].maxvaccdist_t);
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
	if(P.DoRecordInfEvents)
	{
		sprintf(outname,"%s.infevents.csv",OutFile);
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		//fprintf(dat,"t\tind_infectee\tx_infectee\ty_infectee\tt_infector\tind_infector\tx_infector\tiy_infector\n");
		for(i=0;i<*nEvents;i++)
			{
			fprintf(dat,"%lg,%i,%lg,%lg,%lg,%i,%lg,%lg\n",
				InfEventLog[i].t,InfEventLog[i].infectee_ind,InfEventLog[i].infectee_x,InfEventLog[i].infectee_y,InfEventLog[i].t_infector,InfEventLog[i].infector_ind,InfEventLog[i].infector_x,InfEventLog[i].infector_y);
			}
		fclose(dat);
	}
	if(P.DoInfectionTree)
		{
		sprintf(outname,"%s.tree.csv",OutFile);
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "id,onset,infector,infector_onset,district,detected,outcome\n");
		for (i = 0; i < P.N; i++)
			if (Hosts[i].infect_type % INFECT_TYPE_MASK > 0)
				if(Hosts[i].infector>0)
					fprintf(dat, "%i,%i,%i,%i,%i,%i,%i\n", i, Hosts[i].latent_time, Hosts[i].infector, Hosts[Hosts[i].infector].latent_time, AdUnits[Mcells[Hosts[i].mcell].adunit].id, Hosts[i].detected,Hosts[i].to_die);
				else
					fprintf(dat, "%i,%i,%i,%i,%i,%i,%i\n", i, Hosts[i].latent_time, Hosts[i].infector, -1, AdUnits[Mcells[Hosts[i].mcell].adunit].id, Hosts[i].detected, Hosts[i].to_die);
				//fprintf(dat,"%i,%i,%i,%i\n",i,Hosts[i].infector,Hosts[i].infect_type%INFECT_TYPE_MASK,(int) HOST_AGE_YEAR(i));
		fclose(dat);
		}
#if defined(WIN32_BM) || defined(IMAGE_MAGICK)
	static int dm[13]={0,31,28,31,30,31,30,31,31,30,31,30,31};
	int d,m,y,dml,f;
#ifdef WIN32_BM
	if(P.OutputBitmap==1) CloseAvi(avi);
	if((TimeSeries[P.NumSamples-1].extinct)&&(P.OutputNonExtinct))
		{
		sprintf(outname,"%s.avi",OutFile);
		DeleteFile(outname);
		}
#endif
	if(P.OutputBitmap>=1)
		{
	// Generate Google Earth .kml file
#ifdef WIN32_BM
		sprintf(outname,"%s.ge.kml",OutFile); //sprintf(outname,"%s.ge\\%s.kml",OutFileBase,OutFile);
#else	
		sprintf(outname,"%s.ge.kml",OutFile);
#endif
		if(!(dat=fopen(outname,"w"))) 
			{
				ERR_CRITICAL("Unable to open output kml file\n");
			}
		fprintf(dat,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://earth.google.com/kml/2.2\">\n<Document>\n");
		fprintf(dat,"<name>%s</name>\n",OutFile);
		y=2009;
		m=1;
		d=1;
		for(i=0;i<P.NumSamples;i++)
			{
			fprintf(dat,"<GroundOverlay>\n<name>Snapshot %i</name>\n",i+1);
			fprintf(dat,"<TimeSpan>\n<begin>%i-%02i-%02iT00:00:00Z</begin>\n",y,m,d);
			d+=(int) P.SampleStep; // SampleStep has to be an integer here.
			do
				{
				f=1;
				dml=dm[m];
				if((m==2)&&(y%4==0)) dml=29;
				if(d>dml)
					{
					m++;
					if(m>12)
						{
						m-=12;
						y++;
						}
					d-=dml;
					f=0;
					}
				}
			while(!f);
			fprintf(dat,"<end>%i-%02i-%02iT00:00:00Z</end>\n</TimeSpan>\n",y,m,d);
			sprintf(outname,"%s.%i.png",OutFile,i+1);
			fprintf(dat,"<Icon>\n<href>%s</href>\n</Icon>\n",outname);
			fprintf(dat,"<LatLonBox>\n<north>%lg</north>\n<south>%lg</south>\n<east>%lg</east>\n<west>%lg</west>\n</LatLonBox>\n",
				P.SpatialBoundingBox[3],P.SpatialBoundingBox[1],P.SpatialBoundingBox[2],P.SpatialBoundingBox[0]);
			fprintf(dat,"</GroundOverlay>\n");
			}
		fprintf(dat,"</Document>\n</kml>\n");
		fclose(dat);
		}
#endif
}


void SaveSummaryResults(void)
{
	int i,j;
	double c,t;
	FILE *dat;
	char outname[1024];

	c=1/((double) (P.NRactE+P.NRactNE));
	sprintf(outname,"%s.csv",OutFile);
	if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat,"t,S,L,I,R,D,incI,incR,incFC,incFI,incC,incDC,incTC,incETU,cumT,cumTmax,cumTP,cumV,cumVmax,Extinct,Detected,rmsRad,maxRad,vS,vI,vR,vD,vincI,vincR,vincFC,vincFI,tvincC,vincDC,vincTC,vincH,vrmsRad,vmaxRad,%i,%i,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
		P.NRactNE,P.NRactE,P.R0household,P.R0places,P.R0spatial,c*PeakHeightSum,c*PeakHeightSS-c*c*PeakHeightSum*PeakHeightSum,c*PeakTimeSum,c*PeakTimeSS-c*c*PeakTimeSum*PeakTimeSum);
	c=1/((double) P.NRactual);
	//added this for sake of test border control
	//c=1;
	for(i=0;i<P.NumSamples;i++)
		{
		fprintf(dat,"%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,",
			c*TSMean[i].t,c*TSMean[i].S,c*TSMean[i].L,c*TSMean[i].I,c*TSMean[i].R,
			c*TSMean[i].D,c*TSMean[i].incI,c*TSMean[i].incR,c*TSMean[i].incFC,c*TSMean[i].incFI,c*TSMean[i].incC,c*TSMean[i].incDC,c*TSMean[i].incTC,c*TSMean[i].incETU, //added incidence of funeral transmission and hospitalisation
			c*TSMean[i].cumT,TSMean[i].cumTmax,c*TSMean[i].cumTP,c*TSMean[i].cumV,TSMean[i].cumVmax,c*TSMean[i].extinct,c*TSMean[i].detected,c*TSMean[i].rmsRad,c*TSMean[i].maxRad);
		fprintf(dat,"%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
			c*TSVar[i].S-c*c*TSMean[i].S*TSMean[i].S,
			c*TSVar[i].I-c*c*TSMean[i].I*TSMean[i].I,
			c*TSVar[i].R-c*c*TSMean[i].R*TSMean[i].R,
			c*TSVar[i].D-c*c*TSMean[i].D*TSMean[i].D,
			c*TSVar[i].incI-c*c*TSMean[i].incI*TSMean[i].incI,
			c*TSVar[i].incR-c*c*TSMean[i].incR*TSMean[i].incR,
			c*TSVar[i].incD-c*c*TSMean[i].incFC*TSMean[i].incFC,
			c*TSVar[i].incFI-c*c*TSMean[i].incFI*TSMean[i].incFI, //added funeral transmissions
			c*TSVar[i].incC-c*c*TSMean[i].incC*TSMean[i].incC,
			c*TSVar[i].incDC-c*c*TSMean[i].incDC*TSMean[i].incDC, //added detected cases
			c*TSVar[i].incTC-c*c*TSMean[i].incTC*TSMean[i].incTC,
			c*TSVar[i].incETU-c*c*TSMean[i].incETU*TSMean[i].incETU, //added hospitalisation
			c*TSVar[i].rmsRad-c*c*TSMean[i].rmsRad*TSMean[i].rmsRad,
			c*TSVar[i].maxRad-c*c*TSMean[i].maxRad*TSMean[i].maxRad);
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
	if((P.DoAdUnits)&&(P.DoAdunitOutput))
		{
		sprintf(outname,"%s.adunit.csv",OutFile);
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat,"t");
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,",I_%s",AdUnits[i].ad_name);
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,",C_%s",AdUnits[i].ad_name);
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,",DC_%s",AdUnits[i].ad_name); //added detected cases: ggilani 03/02/15
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,",T_%s",AdUnits[i].ad_name);
		if((P.DoHospitalisation)&(P.DoETUByAdUnit))
		{
			for(i=0;i<P.NumAdunits;i++) fprintf(dat,",ETU_%s",AdUnits[i].ad_name); //added hospitalisation
		}
		if(P.DoContactTracing)
		{
			for(i=0;i<P.NumAdunits;i++) fprintf(dat,",CT_%s",AdUnits[i].ad_name); //added contact tracing
		}
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,",%lg",P.PopByAdunit[i][0]);
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,",%lg",P.PopByAdunit[i][1]);
		fprintf(dat,"\n");
		for(i=0;i<P.NumSamples;i++)
			{
			fprintf(dat,"%lg",c*TSMean[i].t);
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,",%lg",c*TSMean[i].incI_adunit[j]);
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,",%lg",c*TSMean[i].incC_adunit[j]);
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,",%lg",c*TSMean[i].incDC_adunit[j]); //added detected cases: ggilani 03/02/15
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,",%lg",c*TSMean[i].cumT_adunit[j]);
			if((P.DoHospitalisation)&(P.DoETUByAdUnit))
			{	
				for(j=0;j<P.NumAdunits;j++)
					fprintf(dat,",%lg",c*TSMean[i].incETU_adunit[j]); //added hospitalisation
			}
			if(P.DoContactTracing)
			{	
				for(j=0;j<P.NumAdunits;j++)
					fprintf(dat,",%lg",c*TSMean[i].incCT_adunit[j]); //added contact tracing
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, ",%lg", c*TSMean[i].incCC_adunit[j]); //added cases who are contacts
			}
			fprintf(dat,"\n");
			}
		fclose(dat);
		sprintf(outname,"%s.adunitVar.csv",OutFile);
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat,"t");
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,",I_%s",AdUnits[i].ad_name);
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,",C_%s",AdUnits[i].ad_name);
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,",DC_%s",AdUnits[i].ad_name); //added detected cases: ggilani 03/02/15
		for(i=0;i<P.NumAdunits;i++) fprintf(dat,",T_%s",AdUnits[i].ad_name);
		if((P.DoHospitalisation)&(P.DoETUByAdUnit))
		{
			for(i=0;i<P.NumAdunits;i++) fprintf(dat,",ETU_%s",AdUnits[i].ad_name); //added hospitalisation
		}
		if(P.DoContactTracing)
		{
			for(i=0;i<P.NumAdunits;i++) fprintf(dat,",CT_%s",AdUnits[i].ad_name); //added contact tracing
		}
		fprintf(dat,"\n");
		for(i=0;i<P.NumSamples;i++)
			{
			fprintf(dat,"%lg",c*TSMean[i].t);
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,",%lg",c*TSVar[i].incI_adunit[j]-c*c*TSMean[i].incI_adunit[j]*TSMean[i].incI_adunit[j]);
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,",%lg",c*TSVar[i].incC_adunit[j]-c*c*TSMean[i].incC_adunit[j]*TSMean[i].incC_adunit[j]);
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,",%lg",c*TSVar[i].incDC_adunit[j]-c*c*TSMean[i].incDC_adunit[j]*TSMean[i].incDC_adunit[j]); //added detected cases: ggilani 03/02/15
			for(j=0;j<P.NumAdunits;j++)
				fprintf(dat,",%lg",c*TSVar[i].cumT_adunit[j]-c*c*TSMean[i].cumT_adunit[j]*TSMean[i].cumT_adunit[j]);
			if((P.DoHospitalisation)&(P.DoETUByAdUnit))
			{
				for(j=0;j<P.NumAdunits;j++)
					fprintf(dat,",%lg",c*TSVar[i].incETU_adunit[j]-c*c*TSMean[i].incETU_adunit[j]*TSMean[i].incETU_adunit[j]); //added hospitalisation
			}
			if(P.DoContactTracing)
			{
				for(j=0;j<P.NumAdunits;j++)
					fprintf(dat,",%lg",c*TSVar[i].incCT_adunit[j]-c*c*TSMean[i].incCT_adunit[j]*TSMean[i].incCT_adunit[j]); //added contact tracing
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, ",%lg", c*TSVar[i].incCC_adunit[j] - c * c*TSMean[i].incCC_adunit[j] * TSMean[i].incCC_adunit[j]); //added cases who are contacts
			}
			fprintf(dat,"\n");
			}
		fclose(dat);
	}

	if(P.EvolResistNumTypes>1)
		{
		sprintf(outname,"%s.resist.xls",OutFile);
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat,"t");
		for(i=0;i<P.EvolResistNumTypes;i++) fprintf(dat,"\tI%i",i);
		for(i=0;i<P.EvolResistNumTypes;i++) fprintf(dat,"\tC%i",i);
		for(i=0;i<P.EvolResistNumTypes;i++) fprintf(dat,"\tT%i",i);
		for(i=0;i<P.EvolResistNumTypes;i++) fprintf(dat,"\tvI%i",i);
		for(i=0;i<P.EvolResistNumTypes;i++) fprintf(dat,"\tvC%i",i);
		for(i=0;i<P.EvolResistNumTypes;i++) fprintf(dat,"\tvT%i",i);
		fprintf(dat,"\n");
		for(i=0;i<P.NumSamples;i++)
			{
			fprintf(dat,"%lg",c*TSMean[i].t);
			for(j=0;j<P.EvolResistNumTypes;j++)
				fprintf(dat,"\t%lg",c*TSMean[i].incI_resist[j]);
			for(j=0;j<P.EvolResistNumTypes;j++)
				fprintf(dat,"\t%lg",c*TSMean[i].incC_resist[j]);
			for(j=0;j<P.EvolResistNumTypes;j++)
				fprintf(dat,"\t%lg",c*TSMean[i].cumT_resist[j]);
			for(j=0;j<P.EvolResistNumTypes;j++)
				fprintf(dat,"\t%lg",c*TSVar[i].incI_resist[j]-c*c*TSMean[i].incI_resist[j]*TSMean[i].incI_resist[j]);
			for(j=0;j<P.EvolResistNumTypes;j++)
				fprintf(dat,"\t%lg",c*TSVar[i].incC_resist[j]-c*c*TSMean[i].incC_resist[j]*TSMean[i].incC_resist[j]);
			for(j=0;j<P.EvolResistNumTypes;j++)
				fprintf(dat,"\t%lg",c*TSVar[i].cumT_resist[j]-c*c*TSMean[i].cumT_resist[j]*TSMean[i].cumT_resist[j]);
			fprintf(dat,"\n");
			}
		fclose(dat);
		}
	if(P.KeyWorkerProphTimeStartBase<P.SampleTime)
		{
		sprintf(outname,"%s.keyworker.xls",OutFile);
		if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat,"t");
		for(i=0;i<2;i++) fprintf(dat,"\tI%i",i);
		for(i=0;i<2;i++) fprintf(dat,"\tC%i",i);
		for(i=0;i<2;i++) fprintf(dat,"\tT%i",i);
		for(i=0;i<2;i++) fprintf(dat,"\tvI%i",i);
		for(i=0;i<2;i++) fprintf(dat,"\tvC%i",i);
		for(i=0;i<2;i++) fprintf(dat,"\tvT%i",i);
		fprintf(dat,"\t%i\t%i\n",P.KeyWorkerNum,P.KeyWorkerIncHouseNum);
		for(i=0;i<P.NumSamples;i++)
			{
			fprintf(dat,"%lg",c*TSMean[i].t);
			for(j=0;j<2;j++)
				fprintf(dat,"\t%lg",c*TSMean[i].incI_keyworker[j]);
			for(j=0;j<2;j++)
				fprintf(dat,"\t%lg",c*TSMean[i].incC_keyworker[j]);
			for(j=0;j<2;j++)
				fprintf(dat,"\t%lg",c*TSMean[i].cumT_keyworker[j]);
			for(j=0;j<2;j++)
				fprintf(dat,"\t%lg",c*TSVar[i].incI_keyworker[j]-c*c*TSMean[i].incI_keyworker[j]*TSMean[i].incI_keyworker[j]);
			for(j=0;j<2;j++)
				fprintf(dat,"\t%lg",c*TSVar[i].incC_keyworker[j]-c*c*TSMean[i].incC_keyworker[j]*TSMean[i].incC_keyworker[j]);
			for(j=0;j<2;j++)
				fprintf(dat,"\t%lg",c*TSVar[i].cumT_keyworker[j]-c*c*TSMean[i].cumT_keyworker[j]*TSMean[i].cumT_keyworker[j]);
			fprintf(dat,"\n");
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
	FILE *dat;
	char outname[1024];

	sprintf(outname,"%s.seeds.csv",OutFile);
	if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat,"%i,%i\n",P.newseed1,P.newseed2);
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
	FILE *dat;
	char outname[1024];

	sprintf(outname,"%s.infevents.csv",OutFile);
	if(!(dat=fopen(outname,"w"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat,"type,t,thread,ind_infectee,cell_infectee,listpos_infectee,adunit_infectee,x_infectee,y_infectee,t_infector,ind_infector,cell_infector\n");
	for(i=0;i<*nEvents;i++)
		{
		fprintf(dat,"%i,%lg,%i,%i,%i,%i,%i,%lg,%lg,%lg,%i,%i\n",
			InfEventLog[i].type,InfEventLog[i].t,InfEventLog[i].thread,InfEventLog[i].infectee_ind,InfEventLog[i].infectee_cell,InfEventLog[i].listpos,InfEventLog[i].infectee_adunit,InfEventLog[i].infectee_x,InfEventLog[i].infectee_y,InfEventLog[i].t_infector,InfEventLog[i].infector_ind,InfEventLog[i].infector_cell);
		}
	fclose(dat);
}

void LoadSnapshot(void)
{
	FILE *dat;
	int i,j,tsi,*CellMemberArray,*CellSuscMemberArray;
	long l;
	long long CM_offset,CSM_offset;
	double t,st;
	int **Array_InvCDF;
	float *Array_tot_prob,**Array_cum_trans,**Array_max_trans;

	if(!(dat=fopen(SnapshotLoadFile,"rb"))) ERR_CRITICAL("Unable to open snapshot file\n");
	fprintf(stderr,"Loading snapshot.");
	if(!(Array_InvCDF=(int **) malloc(P.NCP*sizeof(int *)))) ERR_CRITICAL("Unable to allocate temp cell storage\n");
	if(!(Array_max_trans=(float **) malloc(P.NCP*sizeof(float *)))) ERR_CRITICAL("Unable to temp allocate cell storage\n");
	if(!(Array_cum_trans=(float **) malloc(P.NCP*sizeof(float *)))) ERR_CRITICAL("Unable to temp allocate cell storage\n");
	if(!(Array_tot_prob=(float *) malloc(P.NCP*sizeof(float)))) ERR_CRITICAL("Unable to temp allocate cell storage\n");
	for(i=0;i<P.NCP;i++)
		{
		Array_InvCDF[i]=Cells[i].InvCDF;
		Array_max_trans[i]=Cells[i].max_trans;
		Array_cum_trans[i]=Cells[i].cum_trans;
		Array_tot_prob[i]=Cells[i].tot_prob;
		}

	fread_big((void *) &i,sizeof(int),1,dat); if(i!=P.N) {fprintf(stderr,"Incorrect N (%i %i) in snapshot file.\n",P.N,i);exit(1);}
	fread_big((void *) &i,sizeof(int),1,dat); if(i!=P.NH) ERR_CRITICAL("Incorrect NH in snapshot file.\n");
	fread_big((void *) &i,sizeof(int),1,dat); if(i!=P.NC) {fprintf(stderr,"## %i neq %i\n",i,P.NC); ERR_CRITICAL("Incorrect NC in snapshot file.\n");}
	fread_big((void *) &i,sizeof(int),1,dat); if(i!=P.NCP) ERR_CRITICAL("Incorrect NCP in snapshot file.\n");
	fread_big((void *) &i,sizeof(int),1,dat); if(i!=P.ncw) ERR_CRITICAL("Incorrect ncw in snapshot file.\n");
	fread_big((void *) &i,sizeof(int),1,dat); if(i!=P.nch) ERR_CRITICAL("Incorrect nch in snapshot file.\n");
	fread_big((void *) &l,sizeof(long),1,dat); if(l!=P.seed1) ERR_CRITICAL("Incorrect seed1 in snapshot file.\n");
	fread_big((void *) &l,sizeof(long),1,dat); if(l!=P.seed2) ERR_CRITICAL("Incorrect seed2 in snapshot file.\n");
	fread_big((void *) &t,sizeof(double),1,dat); if(t!=P.TimeStep) ERR_CRITICAL("Incorrect TimeStep in snapshot file.\n");
	fread_big((void *) &(P.SnapshotLoadTime),sizeof(double),1,dat);
	P.NumSamples=1+(int) ceil((P.SampleTime-P.SnapshotLoadTime)/P.SampleStep);
	fprintf(stderr,".");
	fread_big((void *) &CellMemberArray,sizeof(int *),1,dat);
	fprintf(stderr,".");
	fread_big((void *) &CellSuscMemberArray,sizeof(int *),1,dat);
	fprintf(stderr,".");
	CM_offset=State.CellMemberArray-CellMemberArray;
	CSM_offset=State.CellSuscMemberArray-CellSuscMemberArray;

	zfread_big((void *) Hosts,sizeof(person),(size_t) P.N,dat);
	fprintf(stderr,".");
	zfread_big((void *) Households,sizeof(household),(size_t) P.NH,dat);
	fprintf(stderr,".");
	zfread_big((void *) Cells,sizeof(cell),(size_t) P.NC,dat);
	fprintf(stderr,".");
	zfread_big((void *) Mcells,sizeof(microcell),(size_t) P.NMC,dat);
	fprintf(stderr,".");
	zfread_big((void *) State.CellMemberArray,sizeof(int),(size_t) P.N,dat);
	fprintf(stderr,".");
	zfread_big((void *) State.CellSuscMemberArray,sizeof(int),(size_t) P.N,dat);
	fprintf(stderr,".");
	for(i=0;i<P.NC;i++)
		{
		if(Cells[i].n>0)
			{
			Cells[i].members+=CM_offset;
			Cells[i].susceptible+=CSM_offset;
			Cells[i].latent+=CSM_offset;
			Cells[i].infected+=CSM_offset;
			}
		for(j=0;j<MAX_INTERVENTION_TYPES;j++) Cells[i].CurInterv[j]=-1; // turn interventions off in loaded image
		}
	for(i=0;i<P.NMC;i++)
		if(Mcells[i].n>0)
			Mcells[i].members+=CM_offset;

	for(i=0;i<P.NCP;i++)
		{
		Cells[i].InvCDF=Array_InvCDF[i];
		Cells[i].max_trans=Array_max_trans[i];
		Cells[i].cum_trans=Array_cum_trans[i];
		Cells[i].tot_prob=Array_tot_prob[i];
		}
	free(Array_tot_prob);
	free(Array_cum_trans);
	free(Array_max_trans);
	free(Array_InvCDF);
	fprintf(stderr,"\n");
	fclose(dat);
}

void SaveSnapshot(void)
{
	FILE *dat;
	int i=1,j,n,st;
	double sz;

	if(!(dat=fopen(SnapshotSaveFile,"wb"))) ERR_CRITICAL("Unable to open snapshot file\n");

	fwrite_big((void *) &(P.N),sizeof(int),1,dat);
	fprintf(stderr,"## %i\n",i++);
	fwrite_big((void *) &(P.NH),sizeof(int),1,dat);
	fprintf(stderr,"## %i\n",i++);
	fwrite_big((void *) &(P.NC),sizeof(int),1,dat);
	fprintf(stderr,"## %i\n",i++);
	fwrite_big((void *) &(P.NCP),sizeof(int),1,dat);
	fprintf(stderr,"## %i\n",i++);
	fwrite_big((void *) &(P.ncw),sizeof(int),1,dat);
	fprintf(stderr,"## %i\n",i++);
	fwrite_big((void *) &(P.nch),sizeof(int),1,dat);
	fprintf(stderr,"## %i\n",i++);
	fwrite_big((void *) &(P.seed1),sizeof(long),1,dat);
	fprintf(stderr,"## %i\n",i++);
	fwrite_big((void *) &(P.seed2),sizeof(long),1,dat);
	fprintf(stderr,"## %i\n",i++);
	fwrite_big((void *) &(P.TimeStep),sizeof(double),1,dat);
	fprintf(stderr,"## %i\n",i++);
	fwrite_big((void *) &(P.SnapshotSaveTime),sizeof(double),1,dat);
	fprintf(stderr,"## %i\n",i++);
	fwrite_big((void *) &(State.CellMemberArray),sizeof(int *),1,dat);
	fprintf(stderr,"## %i\n",i++);
	fwrite_big((void *) &(State.CellSuscMemberArray),sizeof(int *),1,dat);
	fprintf(stderr,"## %i\n",i++);
	
	zfwrite_big((void *) Hosts,sizeof(person),(size_t) P.N,dat);

	fprintf(stderr,"## %i\n",i++);
	zfwrite_big((void *) Households,sizeof(household),(size_t) P.NH,dat);
	fprintf(stderr,"## %i\n",i++);
	zfwrite_big((void *) Cells,sizeof(cell),(size_t) P.NC,dat);
	fprintf(stderr,"## %i\n",i++);
	zfwrite_big((void *) Mcells,sizeof(microcell),(size_t) P.NMC,dat);
	fprintf(stderr,"## %i\n",i++);

	zfwrite_big((void *) State.CellMemberArray,sizeof(int),(size_t) P.N,dat);
	fprintf(stderr,"## %i\n",i++);
	zfwrite_big((void *) State.CellSuscMemberArray,sizeof(int),(size_t) P.N,dat);
	fprintf(stderr,"## %i\n",i++);

	fclose(dat);
}

void UpdateProbs(int DoPlace)
{
	int j,k,m,same_country,mc,mc_ind,q,p,ci;
	float t;

	if(!DoPlace)
		{
#pragma omp parallel for private(j,k,m,t) schedule(static,500)
		for(j=0;j<P.NCP;j++)
			{
			CellLookup[j]->tot_prob=0;
			CellLookup[j]->S0=CellLookup[j]->S+CellLookup[j]->L+CellLookup[j]->I;
			if((P.DoDeath)&&(!P.DoSIS)) 
				{
				CellLookup[j]->S0+=CellLookup[j]->n/5;
				if((CellLookup[j]->n<100)||(CellLookup[j]->S0>CellLookup[j]->n)) CellLookup[j]->S0=CellLookup[j]->n;
				}
			}
		}
	else
		{
#pragma omp parallel for private(j,k,m,t) schedule(static,500) //added j to private
		for(j=0;j<P.NCP;j++)
			{
			CellLookup[j]->S0=CellLookup[j]->S;
			CellLookup[j]->tot_prob=0;
			}
		}
#pragma omp parallel for private(j,k,m,t,q,p,same_country,mc_ind,mc,ci) schedule(static,500)
	for(j=0;j<P.NCP;j++)
		{

			if((P.DoRoadNetwork)&&(P.DoRoadPopEffect)&&((CellLookup[j]->road_connection)||(CellLookup[0]->road_connection)))
			{
				CellLookup[j]->cum_trans[0]=((float) (CellLookup[0]->S0))*CellLookup[j]->road_access*CellLookup[0]->road_access*CellLookup[j]->max_trans[0]; //scaling distance according to accessibility by road - was accessTo*accessFrom CellLookup[j]->road_access*CellLookup[0]->road_access
				t=((float) CellLookup[0]->n)*CellLookup[j]->road_access*CellLookup[0]->road_access*CellLookup[j]->max_trans[0];
			}
			else if((P.DoCapitalCityEffect)&&(P.DoCapitalCityPopEffect)&&(P.DoAdUnits)&&((CellLookup[j]->capital_city)^(CellLookup[0]->capital_city))) //if we're doing the capital city effect and either cell l or cell m (but not both) are in capital cities
			{
				same_country=0;
				if(CellLookup[j]->capital_city)
				{
					//if cell j contains a capital city district, check to see whether any of cell 0 is in the same country
					ci=CellLookup[0]-Cells;
					mc=(ci/P.nch)*P.NMCL*P.nmch+(ci%P.nch)*P.NMCL;
					for(q=0;q<P.NMCL;q++)
					{
						for(p=0;p<P.NMCL;p++)
						{
							//get index of microcell
							mc_ind=mc+p+q*P.nmch;
							if((AdUnits[Mcells[mc_ind].adunit].id/P.CountryDivisor)==(CellLookup[j]->capital_city/P.CountryDivisor))
							{
								same_country=1;//P.RoadAccessDistance;
							}
						}
					}
				}
				else if (CellLookup[0]->capital_city)
				{
					//if cell 0 contains a capital city district, check to see whether any of cell j is in the same country
					ci=CellLookup[j]-Cells;
					mc=(ci/P.nch)*P.NMCL*P.nmch+(ci%P.nch)*P.NMCL;
					for(q=0;q<P.NMCL;q++)
					{
						for(p=0;p<P.NMCL;p++)
						{
							//get index of microcell
							mc_ind=mc+p+q*P.nmch;
							if((AdUnits[Mcells[mc_ind].adunit].id/P.CountryDivisor)==(CellLookup[0]->capital_city/P.CountryDivisor))
							{
								same_country=1;//P.RoadAccessDistance;
							}
						}
					}
				}
				if(same_country)
				{
					CellLookup[j]->cum_trans[0]=((float) (CellLookup[0]->S0))*P.CapitalCityPopEffect*CellLookup[j]->max_trans[0];
					t=((float) CellLookup[0]->n)*P.CapitalCityPopEffect*CellLookup[j]->max_trans[0];
				}
				else
				{
					CellLookup[j]->cum_trans[0]=((float) (CellLookup[0]->S0))*CellLookup[j]->max_trans[0];
					t=((float) CellLookup[0]->n)*CellLookup[j]->max_trans[0];
				}
			}
			else
			{
				CellLookup[j]->cum_trans[0]=((float) (CellLookup[0]->S0))*CellLookup[j]->max_trans[0];
				t=((float) CellLookup[0]->n)*CellLookup[j]->max_trans[0];
			}


		//CellLookup[j]->cum_trans[0]=((float) (CellLookup[0]->S0))*CellLookup[j]->max_trans[0];
		//t=((float) CellLookup[0]->n)*CellLookup[j]->max_trans[0];
		for(m=1;m<P.NCP;m++)
			{

				if((P.DoRoadNetwork)&&(P.DoRoadPopEffect)&&((CellLookup[j]->road_connection)||(CellLookup[m]->road_connection)))
				{
					CellLookup[j]->cum_trans[m]=CellLookup[j]->cum_trans[m-1]+((float) (CellLookup[m]->S0))*CellLookup[j]->road_access*CellLookup[m]->road_access*CellLookup[j]->max_trans[m]; //scaling distance according to accessibility by road - was accessTo*accessFrom
					t+=((float) CellLookup[m]->n)*CellLookup[j]->road_access*CellLookup[m]->road_access*CellLookup[j]->max_trans[m];
				}
				else if((P.DoCapitalCityEffect)&&(P.DoCapitalCityPopEffect)&&(P.DoAdUnits)&&((CellLookup[j]->capital_city)^(CellLookup[m]->capital_city))) //if we're doing the capital city effect and either cell l or cell m (but not both) are in capital cities
				{
					same_country=0;
					if(CellLookup[j]->capital_city)
					{
						//if cell j contains a capital city district, check to see whether any of cell m is in the same country
						ci=CellLookup[m]-Cells;
						mc=(ci/P.nch)*P.NMCL*P.nmch+(ci%P.nch)*P.NMCL;
						for(q=0;q<P.NMCL;q++)
						{
							for(p=0;p<P.NMCL;p++)
							{
								//get index of microcell
								mc_ind=mc+p+q*P.nmch;
								if((AdUnits[Mcells[mc_ind].adunit].id/P.CountryDivisor)==(CellLookup[j]->capital_city/P.CountryDivisor))
								{
									same_country=1;//P.RoadAccessDistance;
								}
							}
						}
					}
					else if (CellLookup[m]->capital_city)
					{
						//if cell m contains a capital city district, check to see whether any of cell j is in the same country
						ci=CellLookup[j]-Cells;
						mc=(ci/P.nch)*P.NMCL*P.nmch+(ci%P.nch)*P.NMCL;
						for(q=0;q<P.NMCL;q++)
						{
							for(p=0;p<P.NMCL;p++)
							{
								//get index of microcell
								mc_ind=mc+p+q*P.nmch;
								if((AdUnits[Mcells[mc_ind].adunit].id/P.CountryDivisor)==(CellLookup[m]->capital_city/P.CountryDivisor))
								{
									same_country=1;//P.RoadAccessDistance;
								}
							}
						}
					}
					if(same_country)
					{
						CellLookup[j]->cum_trans[m]=CellLookup[j]->cum_trans[m-1]+((float) (CellLookup[m]->S0))*P.CapitalCityPopEffect*CellLookup[j]->max_trans[m];
						t+=((float) CellLookup[m]->n)*P.CapitalCityPopEffect*CellLookup[j]->max_trans[m];
					}
					else
					{
						CellLookup[j]->cum_trans[m]=CellLookup[j]->cum_trans[m-1]+((float) (CellLookup[m]->S0))*CellLookup[j]->max_trans[m];
						t+=((float) CellLookup[m]->n)*CellLookup[j]->max_trans[m];
					}
				}
				else
				{
					CellLookup[j]->cum_trans[m]=CellLookup[j]->cum_trans[m-1]+((float) (CellLookup[m]->S0))*CellLookup[j]->max_trans[m];
					t+=((float) CellLookup[m]->n)*CellLookup[j]->max_trans[m];
				}


			//CellLookup[j]->cum_trans[m]=CellLookup[j]->cum_trans[m-1]+((float) (CellLookup[m]->S0))*CellLookup[j]->max_trans[m];
			//t+=((float) CellLookup[m]->n)*CellLookup[j]->max_trans[m];
			}
		CellLookup[j]->tot_prob=CellLookup[j]->cum_trans[P.NCP-1];
		for(m=0;m<P.NCP;m++)
			CellLookup[j]->cum_trans[m]/=CellLookup[j]->tot_prob;
		CellLookup[j]->tot_prob/=t;
		for(k=m=0;k<=1024;k++)
			{
			while(CellLookup[j]->cum_trans[m]*1024<((float) k)) m++;
			CellLookup[j]->InvCDF[k]=m;
			}
		}
}

void TravelReturnSweep(double t)
{
	int i,j,k,l,n,nr,ner,tn;

	if(floor(1+t+P.TimeStep)!=floor(1+t))
		{
		nr=ner=0;
		k=(int) floor(t);
		l=1+k%MAX_TRAVEL_TIME;
#pragma omp parallel for private(i,j,k,n,tn) reduction(+:nr,ner) schedule(static,1)
		for(tn=0;tn<P.NumThreads;tn++)
			{
			for(j=tn;j<P.Nplace[HOTEL_PLACE_TYPE];j+=P.NumThreads)
				{
				n=Places[HOTEL_PLACE_TYPE][j].n;
				for(k=n-1;k>=0;k--)
					{
					i=Places[HOTEL_PLACE_TYPE][j].members[k];
					if(Hosts[i].Travelling==l)
						{
						n--;
/*						if((n<0)||(Places[HOTEL_PLACE_TYPE][j].members[n]<0)||(Places[HOTEL_PLACE_TYPE][j].members[n]>=P.N))
							{fprintf(stderr,"### %i %i %i %i\n",j,k,n,Places[HOTEL_PLACE_TYPE][j].members[n]);ner++;}
						else if((k<0)||(k>n))
							{fprintf(stderr,"@ %i %i %i %i\n",j,k,n,Places[HOTEL_PLACE_TYPE][j].members[n]);ner++;}
						else
*/
						if(k!=n)
							{Places[HOTEL_PLACE_TYPE][j].members[k]=Places[HOTEL_PLACE_TYPE][j].members[n];}
						nr++;
						if(Hosts[i].PlaceLinks[HOTEL_PLACE_TYPE]!=j)
							{ner++;fprintf(stderr,"(%i %i) ",j,Hosts[i].PlaceLinks[HOTEL_PLACE_TYPE]);}
						Hosts[i].PlaceLinks[HOTEL_PLACE_TYPE]=-1;
						Hosts[i].Travelling=0;
						}
					}
				Places[HOTEL_PLACE_TYPE][j].n=n;
				}
			}
		fprintf(stderr," d=%i e=%i>",nr,ner);
		}
}

void TravelDepartSweep(double t)
{
	int c,i,i2,j,k,l,d,d2,m,n,f,f2,f3,mps,nld,nad,nsk,tn,bm,hp;
	double s,s2,nl;
	cell *ct;

	if(floor(1+t-P.TimeStep)!=floor(1+t))
		{
		bm=((P.DoBlanketMoveRestr)&&(t>=P.MoveRestrTimeStart)&&(t<P.MoveRestrTimeStart+P.MoveRestrDuration));
		mps=2*((int) P.PlaceTypeMeanSize[HOTEL_PLACE_TYPE])-P.NumThreads-1;
		k=(int) floor(t);
		d=k%MAX_TRAVEL_TIME;
		nad=nld=nsk=0;
#pragma omp parallel for private(i,i2,j,k,l,d2,m,n,s,f,f2,f3,tn,hp) reduction(+:nad,nsk) schedule(static,1)
		for(tn=0;tn<P.NumThreads;tn++)
			for(i=tn;i<P.Nairports;i+=P.NumThreads)
				if((Airports[i].total_traffic>0)&&(Airports[i].num_mcell>0))
					{
					s=Airports[i].total_traffic;
					if((t>P.AirportCloseTimeStart)&&(t<P.AirportCloseTimeStart+P.AirportCloseTimeStartBase))
						s*=P.AirportCloseEffectiveness;
					n=(s>0)?((int) ignpoi_mt(ranf_seed, (double) s,tn)):0;
					f3=0;
					j=0;
					while(j<n)
						{
						s=ranf_mt(ranf_seed, tn);
						l=Airports[i].Inv_DestMcells[(int) floor(s*1024)];
						while(Airports[i].DestMcells[l].prob<s) l++;
						l=Airports[i].DestMcells[l].id;
						k=(int) (ranf_mt(ranf_seed, tn)*((double) Mcells[l].n));
						i2=Mcells[l].members[k];
						if((abs(Hosts[i2].inf)<2)&&(Hosts[i2].inf!=-2))
							{
							d2=HOST_AGE_GROUP(i2);
							if((P.RelativeTravelRate[d2]==1)||(ranf_mt(ranf_seed, tn)<P.RelativeTravelRate[d2]))
								{
								f2=1;
#pragma omp critical
									{
									if(Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE]==-1)
										{
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE]=-2;
										f2=0;
										}
									}
								if(!f2)
									{
									s=ranf_mt(ranf_seed, tn);
									l=Airports[i].Inv_prop_traffic[(int) floor(s*128)];
									while(Airports[i].prop_traffic[l]<s) l++;
									k=Airports[i].conn_airports[l];
									f2=0;
									if(bm)
										{
										if(dist2_raw(Airports[i].loc_x,Airports[i].loc_y,Airports[k].loc_x,Airports[k].loc_y)>P.MoveRestrRadius2)
											{
											if(ranf_mt(ranf_seed, tn)>P.MoveRestrEffect) 
												{
												f2=1;
												nsk++;
												j++;
#pragma omp critical
												Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE]=-1;
												}
											}
										}
									if(!f2)
										{
										f=1;
										f2=0;
										do
											{
											s=ranf_mt(ranf_seed, tn);
											m=Airports[k].Inv_DestPlaces[(int) floor(s*1024)];
											while(Airports[k].DestPlaces[m].prob<s) m++;
											l=Airports[k].DestPlaces[m].id;
#pragma omp critical
												{
												if((hp=Places[HOTEL_PLACE_TYPE][l].n)<mps)
													{
													f=0;
													Places[HOTEL_PLACE_TYPE][l].n++;
													}
												}
											if(!f)
												{
												f3=0;
												Places[HOTEL_PLACE_TYPE][l].members[hp]=i2;
												d2=(d+P.InvJourneyDurationDistrib[(int) (ranf_mt(ranf_seed, tn)*1024.0)])%MAX_TRAVEL_TIME;
												Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE]=l;
												Hosts[i2].Travelling=1+d2;
												nad++;
												j++;
												}
											f2++;
											}
										while((f)&&(f2<300));
										if(f)
											{
#pragma omp critical
											Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE]=-1;
											if(++f3>100)
												{j++;nsk++;}
											}
										}
									}
								}
							}
						else
							j++;
						}
					}
				fprintf(stderr,"<ar=%i as=%i",nad,nsk);
		nl=((double) P.PlaceTypeMeanSize[HOTEL_PLACE_TYPE])*P.HotelPropLocal/P.MeanLocalJourneyTime;
		nsk=0;
#pragma omp parallel for private(c,i,i2,j,l,d2,m,n,s,s2,ct,f,f2,f3,tn,hp) reduction(+:nld,nsk) schedule(static,1)
		for(tn=0;tn<P.NumThreads;tn++)
			for(i=tn;i<P.Nplace[HOTEL_PLACE_TYPE];i+=P.NumThreads)
				{
				c=((int) (Places[HOTEL_PLACE_TYPE][i].loc_x/P.cwidth))*P.nch+((int) (Places[HOTEL_PLACE_TYPE][i].loc_y/P.cheight));
				n=(int) ignpoi_mt(ranf_seed, nl*Cells[c].tot_prob, tn);
				if(Places[HOTEL_PLACE_TYPE][i].n+n>mps) 
					{
					nsk+=(Places[HOTEL_PLACE_TYPE][i].n+n-mps);
					n=mps-Places[HOTEL_PLACE_TYPE][i].n;
					}
				for(j=0;j<n;j++)
					{
					do
						{
						f=0;
						s=ranf_mt(ranf_seed, tn);
						l=Cells[c].InvCDF[(int) floor(s*1024)];
						while(Cells[c].cum_trans[l]<s) l++;
						ct=CellLookup[l];
						m=(int) (ranf_mt(ranf_seed, tn)*((double) ct->S0));
						if(m<(ct->S+ct->L))
							{
							i2=ct->susceptible[m];
							d2=HOST_AGE_GROUP(i2);
							f3=0;
							if((Hosts[i2].Travelling==0)&&((P.RelativeTravelRate[d2]==1)||(ranf_mt(ranf_seed, tn)<P.RelativeTravelRate[d2])))
								{
#pragma omp critical
									{if(Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE]==-1) {Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE]=-2;f3=1;}}
								}
							if(f3)
								{
								s2=dist2_raw(Households[Hosts[i2].hh].loc_x,Households[Hosts[i2].hh].loc_y,Places[HOTEL_PLACE_TYPE][i].loc_x,Places[HOTEL_PLACE_TYPE][i].loc_y);
								f2=1;
								if((bm)&&(s2>P.MoveRestrRadius2))
									{
									if(ranf_mt(ranf_seed, tn)>=P.MoveRestrEffect)
										{
#pragma omp critical
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE]=-1;
										nsk++;
										f2=0;
										}
									}
								if(f2)
									{
									s=numKernel(s2)/Cells[c].max_trans[l];
									if(ranf_mt(ranf_seed, tn)>=s)
										{
#pragma omp critical
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE]=-1;
										f=1;
										}
									else
										{
										d2=(d+P.InvLocalJourneyDurationDistrib[(int) (ranf_mt(ranf_seed, tn)*1024.0)])%MAX_TRAVEL_TIME;
										hp=Places[HOTEL_PLACE_TYPE][i].n;
										Places[HOTEL_PLACE_TYPE][i].n++;
										Places[HOTEL_PLACE_TYPE][i].members[hp]=i2;
										Hosts[i2].Travelling=1+d2;
										nld++;
#pragma omp critical
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE]=i;
										}
									}
								}
							else
								f=1;
							}
						else
							nsk++;
						}
					while(f);
					}
				}
		fprintf(stderr," l=%i ls=%i ",nld,nsk);
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

double CalcHouseInf(int j, unsigned short int ts)
{
	return (HOST_ISOLATED(j)?P.CaseIsolationHouseEffectiveness:1.0)
			*(HOST_QUARANTINED(j)?P.HQuarantineHouseEffect:1.0)
			*P.HouseholdDenomLookup[Households[Hosts[j].hh].nhr-1]*CalcPersonInf(j,ts);
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
	return (HOST_ISOLATED(j)?P.CaseIsolationEffectiveness:1.0)
			*(HOST_QUARANTINED(j)?P.HQuarantineSpatialEffect:1.0)
			*((Hosts[j].inf==-2)?P.SymptSpatialContactRate:1.0)
			*P.RelativeSpatialContact[HOST_AGE_GROUP(j)]
			*CalcPersonInf(j,ts); 		/*	*Hosts[j].spatial_norm */
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

double CalcHouseSusc(int ai, unsigned short int ts,int infector,int tn)
{
	return CalcPersonSusc(ai,ts,infector,tn)
		*((Mcells[Hosts[ai].mcell].socdist==2)?((Hosts[ai].esocdist_comply)?P.ESocDistHouseholdEffect:P.SocDistHouseholdEffect):1.0);
}

double CalcPlaceSusc(int ai,int k, unsigned short int ts,int infector,int tn)
{
	return CalcPersonSusc(ai,ts,infector,tn)*(HOST_QUARANTINED(ai)?P.HQuarantinePlaceEffect[k]:1.0)
		*((Mcells[Hosts[ai].mcell].socdist==2)?((Hosts[ai].esocdist_comply)?P.ESocDistPlaceEffect[k]:P.SocDistPlaceEffect[k]):1.0);
}

double CalcSpatialSusc(int ai, unsigned short int ts,int infector,int tn)
{
	return CalcPersonSusc(ai,ts,infector,tn)*(HOST_QUARANTINED(ai)?P.HQuarantineSpatialEffect:1.0)
		*((Mcells[Hosts[ai].mcell].socdist==2)?((Hosts[ai].esocdist_comply)?P.ESocDistSpatialEffect:P.SocDistSpatialEffect):1.0);
}


double CalcPersonSusc(int ai, unsigned short int ts,int infector,int tn)
{
	double suscdrop;
	UpdatePersonDemog(ai,tn);

	if ((Hosts[ai].keyworker) && ((Hosts[ai].vacc_start_time <= 0) || (HOST_TO_BE_VACCED(ai)&&Hosts[ai].revacc)))
	{
		suscdrop = P.HCWVaccSuscDrop;
	}
	else
	{
		suscdrop = (HOST_VACCED(ai) ? (HOST_VACCED_SWITCH(ai) ? P.VaccSuscDrop2 : P.VaccSuscDrop) : 1.0);
	}

	return P.WAIFW_Matrix[HOST_AGE_GROUP(ai)][HOST_AGE_GROUP(infector)]
		*P.AgeSusceptibility[HOST_AGE_GROUP(ai)]*Hosts[ai].susc
		*(HOST_TREATED(ai)?P.EvolResistRelProphSusc[Hosts[infector].resist]:1.0)
		*suscdrop*PrevalenceDepTransmission(ai)
		*(((Hosts[ai].hcw||Hosts[ai].flw)&&P.OutbreakDetected)?P.RelSuscPPE:1.0);
}

void EquilibPersonDemog(int ai,int tn)
{
#ifdef NEW_AGE_MODEL
	int b,l;
	double f;

/* Lazy evaluation of deaths - only suited to constant populations and SIS model.
*/
	if(P.DoDeath)
		{
		if((Hosts[ai].birth_time+(int) Hosts[ai].life_expectancy)<P.ts_age)
			{
			b=Hosts[ai].birth_time;
			l=(int) Hosts[ai].life_expectancy;
			while((b+l)<P.ts_age) 
				{
				b+=l;
				f=P.InvLifeExpecDist[0][(int) ceil(1000.0*ranf_mt(ranf_seed, tn))];
				l=((int) (f*P.TimeStepsPerYear));
				}
			Hosts[ai].birth_time=b;
			Hosts[ai].life_expectancy=(unsigned short int) l;
			}
		}
#endif
}

void UpdatePersonDemog(int ai,int tn)
{
#ifdef NEW_AGE_MODEL
	int b,l,i,j,k;
	double f;

/* Lazy evaluation of deaths - only suited to constant populations and SIS model. This only checks to see if someone has died
if they are susceptible and have been contacted. If they are dead they are reborn and susceptibility is reset.
*/
	if(P.DoDeath)
		{
		if(P.DoSIS)
			{
			if((Hosts[ai].birth_time+(int) Hosts[ai].life_expectancy)<P.ts_age)
				{
				b=Hosts[ai].birth_time;l=Hosts[ai].life_expectancy;
				while((b+l)<P.ts_age) 
					{
					b+=l;
					f=P.InvLifeExpecDist[0][(int) ceil(1000.0*ranf_mt(ranf_seed, tn))];
					l=((int) (f*P.TimeStepsPerYear));
					StateT[tn].cumD++;
					StateT[tn].cumDa[HOST_AGE_GROUP(ai)]++;
					}
				if((b+P.usRoutineImmunisationMinAge<P.usRoutineImmunisationStartTime)||(P.RoutineImmunisationEffectiveCoverage==0))
					f=1;
				else
					f=(ranf_mt(ranf_seed, tn)<P.RoutineImmunisationEffectiveCoverage)?0:1;
#ifdef FRESSCA
				if(P.DoDistributionVaccination)
					{
					//Check if dead ai in vacc queue
					i=Hosts[ai].mcell;
					j=Hosts[ai].vacc_queue_pos;
					if(j>=0)
#pragma omp critical (dvacc)
						{
						Mcells[i].dvacc_queue[j]=-1;
						}
					if(!f)
						{
						// Add ai to vaccination queue
						j=Mcells[i].dvacc_count;
						k=(Mcells[i].ndvacc_queue+j)%Mcells[i].n;
#pragma omp critical (dvacc)
							{
							Mcells[i].dvacc_queue[k]=ai;
							Mcells[i].dvacc_min_vacc_time[k]=b+P.usRoutineImmunisationMinAge;
							Mcells[i].dvacc_expiry_time[k]=b+P.usRoutineImmunisationMaxAge;
							Mcells[i].dvacc_count++;
							Hosts[ai].vacc_queue_pos=k;
							}
						f=1;
						}
					}
#endif
#pragma omp critical (hostbirth)
					{
					Hosts[ai].birth_time=b;
					Hosts[ai].life_expectancy=l;
					Hosts[ai].infector=-1;
					Hosts[ai].susc=f;
					}
				}
			}
//		else
//			UpdateSuscPersonDemogSIR(ai,tn);
		}
#endif
}

void UpdatePersonDemogSIR(int ai,int tn)
{
#ifdef NEW_AGE_MODEL
	int i,j,k,b,l,f;
	double coverage;

/* Update demog for SIR model. This checks to see if someone has died and resets birth date and moves them back to susceptible
*/
	if(abs(Hosts[ai].inf)>=3)
		{
		if((Hosts[ai].birth_time+(int) Hosts[ai].life_expectancy)<P.ts_age)
			{
			b=Hosts[ai].birth_time;l=Hosts[ai].life_expectancy;
			do 
				{
				StateT[tn].cumD++;
				StateT[tn].cumDa[l/((int) (P.TimeStepsPerYear*AGE_GROUP_WIDTH))]++;
				b+=l;
				f=P.InvLifeExpecDist[0][(int) ceil(1000.0*ranf_mt(ranf_seed, tn))];
				l=((int) (f*P.TimeStepsPerYear));
				}
			while((b+l)<P.ts_age);
			Hosts[ai].birth_time=b;
			Hosts[ai].life_expectancy=l;
			Hosts[ai].infector=-1;
#ifdef FRESSCA
			if(P.DoDistributionVaccination)
				{
				i=Hosts[ai].mcell;
				if(b+P.usRoutineImmunisationMinAge>=P.usRoutineImmunisationStartTime)
					{
					coverage=((P.DistribNetCountry==-1)||(Mcells[i].country==P.DistribNetCountry))?P.RoutineImmunisationEffectiveCoverage:P.RoutineCoverageNonDistrCountry;
					if(coverage==0)
						f=1;
					else if (coverage==1)
						f=0;
					else
						f=(ranf_mt(ranf_seed, tn)<coverage)?0:1;
					}
				else
					f=1;
				DoReborn(ai);
				//Check if dead ai in vacc queue
				j=Hosts[ai].vacc_queue_pos;
				if(j>=0)
					{
					Mcells[i].dvacc_queue[j]=-1;
					Hosts[ai].vacc_queue_pos=-1;
					}
				if(!f)
					{
					// Add ai to vaccination queue
					j=Mcells[i].dvacc_count;
					k=(Mcells[i].ndvacc_queue+j)%Mcells[i].n;
					Mcells[i].dvacc_queue[k]=ai;
					Mcells[i].dvacc_min_vacc_time[k]=b+P.usRoutineImmunisationMinAge;
					Mcells[i].dvacc_expiry_time[k]=b+P.usRoutineImmunisationMaxAge;
					Mcells[i].dvacc_count++;
					Hosts[ai].vacc_queue_pos=k;
					}
				}
			else
#endif
				{
				if((b+P.usRoutineImmunisationMinAge<P.usRoutineImmunisationStartTime)||(P.RoutineImmunisationEffectiveCoverage==0))
					f=1;
				else
					f=(ranf_mt(ranf_seed, tn)<P.RoutineImmunisationEffectiveCoverage)?0:1;
				if(f) DoReborn(ai);
				}
			}
		}
#endif
}
void UpdateSuscPersonDemogSIR(int ai,int tn)
{
#ifdef NEW_AGE_MODEL
	int i,j,k,b,l,f;
	double coverage;

/* Update demog for SIR model. This checks to see if someone has died and resets birth date and moves them back to susceptible.
*/
	if(Hosts[ai].inf==0)
		{
		if((Hosts[ai].birth_time+Hosts[ai].life_expectancy)<P.ts_age)
			{
			b=Hosts[ai].birth_time;l=Hosts[ai].life_expectancy;
			do
				{
				StateT[tn].cumD++;
				StateT[tn].cumDa[l/((int) (P.TimeStepsPerYear*AGE_GROUP_WIDTH))]++;
				b+=l;
				f=P.InvLifeExpecDist[0][(int) ceil(1000.0*ranf_mt(ranf_seed, tn))];
				l=((int) (f*P.TimeStepsPerYear));
				}
			while((b+l)<P.ts_age);
			Hosts[ai].birth_time=b;
			Hosts[ai].life_expectancy=l;
			Hosts[ai].infector=-1;
#ifdef FRESSCA
			if(P.DoDistributionVaccination)
				{
				//Check if dead ai in vacc queue
				i=Hosts[ai].mcell;
				if(b+P.usRoutineImmunisationMinAge>=P.usRoutineImmunisationStartTime)
					{
					coverage=((P.DistribNetCountry==-1)||(Mcells[i].country==P.DistribNetCountry))?P.RoutineImmunisationEffectiveCoverage:P.RoutineCoverageNonDistrCountry;
					if(coverage==0)
						f=1;
					else if (coverage==1)
						f=0;
					else
						f=(ranf_mt(ranf_seed, tn)<coverage)?0:1;
					}
				else
					f=1;
				j=Hosts[ai].vacc_queue_pos;
				if(j>=0)
					{
					Hosts[ai].vacc_queue_pos=-1;
					Mcells[i].dvacc_queue[j]=-1;
					}
				if(!f)
					{
				// Add ai to vaccination queue
					j=Mcells[i].dvacc_count;
					k=(Mcells[i].ndvacc_queue+j)%Mcells[i].n;
					Mcells[i].dvacc_queue[k]=ai;
					Mcells[i].dvacc_min_vacc_time[k]=(unsigned short int) (b+((int) P.usRoutineImmunisationMinAge));
					Mcells[i].dvacc_expiry_time[k]=(unsigned short int) (b+((int) P.usRoutineImmunisationMaxAge));
					Mcells[i].dvacc_count++;
					Hosts[ai].vacc_queue_pos=k;
					}
				}
			else
#endif
				{
				if((b+P.usRoutineImmunisationMinAge<P.usRoutineImmunisationStartTime)||(P.RoutineImmunisationEffectiveCoverage==0))
					f=1;
				else
					f=(ranf_mt(ranf_seed, tn)<P.RoutineImmunisationEffectiveCoverage)?0:1;
				if(!f) DoImmune(ai);
				}
			}
		}
#endif
}

#ifdef FRESSCA
void UpdateVaccStatus(int ai,int tn)
{
	int i,j,k,b,l,f;

	if(Hosts[ai].inf==0)
		{
		b=Hosts[ai].birth_time;
		l=b+P.usRoutineImmunisationMinAge;
		j=P.usRoutineImmunisationStartTime;
		if((l<j)||(P.RoutineImmunisationEffectiveCoverage==0))
			f=1;
		else
			f=(ranf_mt(ranf_seed, tn)<P.RoutineImmunisationEffectiveCoverage)?0:1;
		if(!f)
			{
		// Add ai to vaccination queue
			i=Hosts[ai].mcell;
			j=Mcells[i].dvacc_count;
			k=(Mcells[i].ndvacc_queue+j)%Mcells[i].n;
			Mcells[i].dvacc_queue[k]=ai;
			Mcells[i].dvacc_min_vacc_time[k]=(unsigned short int) l;
			l=b+P.usRoutineImmunisationMaxAge;
			Mcells[i].dvacc_expiry_time[k]=(unsigned short int) l;
			Mcells[i].dvacc_count++;
			Hosts[ai].vacc_queue_pos=k;
			}
		}
}
#endif

/* function: PrevalenceDepTransmission(int ai)
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
double PrevalenceDepTransmission(int ai)
{
	if(P.DoPrevDependTrans)
	{
		if(P.DoPrevDepTransTotalCases)
		{
			if(((double)(Cells[Hosts[ai].pcell].D+Cells[Hosts[ai].pcell].R))/(double)Cells[Hosts[ai].pcell].n>=P.PrevDependTransThresh)
			{
				return P.PrevDependRelativeSusc;
			}
			else
			{
				return 1.0;
			}
		}
		else if(P.DoPrevDepTransCurrInf)
		{
			if((double)Cells[Hosts[ai].pcell].I/(double)Cells[Hosts[ai].pcell].n>=P.PrevDependTransThresh)
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

void InfectSweepOld(double t, int run) //added int as argument to InfectSweepOld to record run number: ggilani - 15/10/14
{
	int i,j,k,l,m,n,i2,b,i3,f,f2,tn,cq,bm;
	double seasonality,sbeta,hbeta,s,s2,s3,s4,fp;
	cell *c,*ct;
	microcell *mi,*mt,*mp;
	unsigned short int ts; 
	
	if(!P.DoSeasonality)
		seasonality=1.0;
	else
		seasonality=P.Seasonality[((int) t)%DAYS_PER_YEAR];
	ts=(unsigned short int) (P.TimeStepsPerDay*t);
	fp=P.TimeStep/(1-P.FalsePositiveRate);
	sbeta=seasonality*fp*P.LocalBeta;
	hbeta=(P.DoHouseholds)?(seasonality*fp*P.HouseholdTrans):0;
	bm=((P.DoBlanketMoveRestr)&&(t>=P.MoveRestrTimeStart)&&(t<P.MoveRestrTimeStart+P.MoveRestrDuration));
#pragma omp parallel for private(j,k,l,m,n,i2,b,i3,f,f2,s,s2,s3,s4,c,ct,tn,mi,mt,mp,cq) schedule(static,1)
	for(tn=0;tn<P.NumThreads;tn++)
		for(b=tn;b<P.NCP;b+=P.NumThreads)
			{	
			c=CellLookup[b];
			if(hbeta>0)
				{
				for(j=0;j<c->I;j++)
					if((Households[Hosts[c->infected[j]].hh].nh>1)&&(!(Hosts[c->infected[j]].nc_plus_hh_disabled & HH_DISABLED))&&(!Hosts[c->infected[j]].Travelling))
						{
						l=Households[Hosts[c->infected[j]].hh].FirstPerson;
						m=l+Households[Hosts[c->infected[j]].hh].nh;
						s3=hbeta*CalcHouseInf(c->infected[j],ts);
						f=0;
						for(i3=l;(i3<m)&&(!f);i3++)
							for(i2=0;(i2<P.PlaceTypeNum)&&(!f);i2++)
								if(Hosts[i3].PlaceLinks[i2]>=0)
									f=PLACE_CLOSED(i2,Hosts[i3].PlaceLinks[i2]);
						if(f) {s3*=P.PlaceCloseHouseholdRelContact;}/* NumPCD++;}*/
						for(i3=l;i3<m;i3++)
							{
							if((Hosts[i3].inf==0)&& (!(Hosts[i3].nc_plus_hh_disabled & HH_DISABLED))&& (!Hosts[i3].Travelling))
								{
								s=s3*CalcHouseSusc(i3,ts,c->infected[j],tn);
								if(ranf_mt(ranf_seed, tn)<s)
									{
									cq=Hosts[i3].pcell%P.NumThreads;
									if((Hosts[i3].infector==-1)&&(StateT[tn].n_queue[cq]<P.InfQueuePeakLength))
										{
										if((P.FalsePositiveRate>0)&&(ranf_mt(ranf_seed, tn)<P.FalsePositiveRate))
											Hosts[i3].infector=Hosts[i3].infect_type=-1;
										else
											{
											Hosts[i3].infector=c->infected[j];
											Hosts[i3].infect_type=1+INFECT_TYPE_MASK*(1+Hosts[c->infected[j]].infect_type/INFECT_TYPE_MASK);
											}
										StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++]=i3;
										}
									}
								}
							}
						}
				}
			if(P.DoPlaces)
				{
				for(j=0;j<c->I;j++)
				  if(!HOST_ABSENT(c->infected[j]))
					{
					mi=Mcells+Hosts[c->infected[j]].mcell;
					for(k=0;k<P.PlaceTypeNum;k++)
						{
						l=Hosts[c->infected[j]].PlaceLinks[k];
						if((l>=0)&&(!PLACE_CLOSED(k,l)))
							{
							s3=fp*seasonality*CalcPlaceInf(c->infected[j],k,ts);
							mp=Mcells+Places[k][l].mcell;
							if(bm)
								{
								if((dist2_raw(Households[Hosts[c->infected[j]].hh].loc_x,Households[Hosts[c->infected[j]].hh].loc_y,
									Places[k][l].loc_x,Places[k][l].loc_y)>P.MoveRestrRadius2))
										s3*=P.MoveRestrEffect;
								}
							else if((mi->moverest!=mp->moverest)&&((mi->moverest==2)||(mp->moverest==2)))
								s3*=P.MoveRestrEffect;

							if((k!=HOTEL_PLACE_TYPE)&&(!Hosts[c->infected[j]].Travelling))
								{
								i2=(Hosts[c->infected[j]].PlaceGroupLinks[k]);
								s4=s3;
								//s4=s3*(1-P.PlaceTypePropBetweenGroupLinks[k]*P.PlaceTypeGroupSizeParam1[k]/((double) Places[k][l].n));
								if(s4<0) 
								{
									fprintf(stderr,"@@@ %lg\n",s4);
									exit(1);
								}
								else if(s4>=1)
									n=Places[k][l].group_size[i2];
								else
									n=(int) ignbin_mt(ranf_seed, (long) Places[k][l].group_size[i2],s4,tn);
								if(n>0) SampleWithoutReplacement(tn,n,Places[k][l].group_size[i2]);
								for(m=0;m<n;m++)
									{
									i3=Places[k][l].members[Places[k][l].group_start[i2]+SamplingQueue[tn][m]];
									if((Hosts[i3].inf==0)&&(!HOST_ABSENT(i3)))
										{
										mt=Mcells+Hosts[i3].mcell;
										ct=Cells+Hosts[i3].pcell;
										s=CalcPlaceSusc(i3,k,ts,c->infected[j],tn);
										if(bm)
											{
											if((dist2_raw(Households[Hosts[i3].hh].loc_x,Households[Hosts[i3].hh].loc_y,
												Places[k][l].loc_x,Places[k][l].loc_y)>P.MoveRestrRadius2))
												s*=P.MoveRestrEffect;
											}
										else if((mt->moverest!=mp->moverest)&&((mt->moverest==2)||(mp->moverest==2)))
											s*=P.MoveRestrEffect;
										if((s==1)||(ranf_mt(ranf_seed, tn)<s))
											{
											cq=Hosts[i3].pcell%P.NumThreads;
											if((Hosts[i3].infector==-1)&&(StateT[tn].n_queue[cq]<P.InfQueuePeakLength))
												{
												if((P.FalsePositiveRate>0)&&(ranf_mt(ranf_seed, tn)<P.FalsePositiveRate))
													Hosts[i3].infector=Hosts[i3].infect_type=-1;
												else
													{
													Hosts[i3].infector=c->infected[j];
													Hosts[i3].infect_type=2+k+INFECT_TYPE_MASK*(1+Hosts[c->infected[j]].infect_type/INFECT_TYPE_MASK);
													}
												StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++]=i3;
												}
											}
										}
									}
								}
							if((k==HOTEL_PLACE_TYPE)||(!Hosts[c->infected[j]].Travelling))
								{
								s3*=P.PlaceTypePropBetweenGroupLinks[k]*P.PlaceTypeGroupSizeParam1[k]/((double) Places[k][l].n);
								if(s3<0) {
									fprintf(stderr,"@@@ %lg\n",s3);
									exit(1);
								}
								else if(s3>=1)
									n=Places[k][l].n;
								else
									n=(int) ignbin_mt(ranf_seed, (long) Places[k][l].n,s3,tn);
								if(n>0) SampleWithoutReplacement(tn,n,Places[k][l].n);		
								for(m=0;m<n;m++)
									{
									i3=Places[k][l].members[SamplingQueue[tn][m]];
									if((Hosts[i3].inf==0)&&(!HOST_ABSENT(i3)))
										{
										mt=Mcells+Hosts[i3].mcell;
										ct=Cells+Hosts[i3].pcell;
										s=CalcPlaceSusc(i3,k,ts,c->infected[j],tn); 
										if(bm)
											{
											if((dist2_raw(Households[Hosts[i3].hh].loc_x,Households[Hosts[i3].hh].loc_y,
												Places[k][l].loc_x,Places[k][l].loc_y)>P.MoveRestrRadius2))
												s*=P.MoveRestrEffect;
											}
										else if((mt->moverest!=mp->moverest)&&((mt->moverest==2)||(mp->moverest==2)))
											s*=P.MoveRestrEffect;
										if((s==1)||(ranf_mt(ranf_seed, tn)<s))
											{
											cq=Hosts[i3].pcell%P.NumThreads;
											if((Hosts[i3].infector==-1)&&(StateT[tn].n_queue[cq]<P.InfQueuePeakLength))
												{
												if((P.FalsePositiveRate>0)&&(ranf_mt(ranf_seed, tn)<P.FalsePositiveRate))
													Hosts[i3].infector=Hosts[i3].infect_type=-1;
												else
													{
													Hosts[i3].infector=c->infected[j];
													Hosts[i3].infect_type=2+k+NUM_PLACE_TYPES+INFECT_TYPE_MASK*(1+Hosts[c->infected[j]].infect_type/INFECT_TYPE_MASK);
													}
												StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++]=i3;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			if(sbeta>0)
				{
				i2=c->I;
				s3=0;
				for(j=0;j<i2;j++)
					{
					if(Hosts[c->infected[j]].Travelling)
						{s2=0;f=0;}
					else
						{
						s2=CalcSpatialInf(c->infected[j],ts);
						f=0;
						if(P.DoPlaces)
						  for(i3=0;(i3<P.PlaceTypeNum)&&(!f);i3++)
							if(Hosts[c->infected[j]].PlaceLinks[i3]>=0)
								f=PLACE_CLOSED(i3,Hosts[c->infected[j]].PlaceLinks[i3]);
						}
					if(f)
						{
						s2*=P.PlaceCloseSpatialRelContact;
						/* NumPCD++; */
						s3+=s2;
						StateT[tn].cell_inf[j]=-s3;
						}
					else
						{
						s3+=s2;
						StateT[tn].cell_inf[j]=s3;
						}
					}
				if(s3>0)
					{
					n=(int) ignpoi_mt(ranf_seed, s3*sbeta*((double) c->tot_prob),tn);
					if(n>0)
						{
						for(j=0;j<i2-1;j++) StateT[tn].cell_inf[j]/=s3;
						StateT[tn].cell_inf[i2-1]=(StateT[tn].cell_inf[i2-1]<0)?-1:1;
						}
					for(k=0;k<n;k++)
						{
						if(i2==1) 
							j=0;
						else
							{
							s=ranf_mt(ranf_seed, tn);
							j=m=i2/2;
							f=1;
							do
								{
								if(m>1) m/=2;
								if((j>0)&&(fabs(StateT[tn].cell_inf[j-1])>=s))
									{
									j-=m;
									if(j==0)
										f=0;
//									else
//										f=(fabs(StateT[tn].cell_inf[j-1])>=s);
									}
								else if((j<i2-1)&&(fabs(StateT[tn].cell_inf[j])<s))
									{
									j+=m;
									if(j==i2-1)
										f=0;
//									else
//										f=(fabs(StateT[tn].cell_inf[j])<s);
									}
								else
									{
/*									fprintf(stderr,"##%i %i %i %i: %lg %lg %lg\n",i2,j,m,f,s,fabs(StateT[tn].cell_inf[j-1]),fabs(StateT[tn].cell_inf[j]));
*/									f=0;
									}
								}
							while(f);
/*							if(fabs(StateT[tn].cell_inf[j])<s) fprintf(stderr,"##%i %i %i %i: %lg %lg %lg\n",i2,j,m,f,s,fabs(StateT[tn].cell_inf[j-1]),fabs(StateT[tn].cell_inf[j]));
							if((j>0)&&(fabs(StateT[tn].cell_inf[j-1])>=s)) fprintf(stderr,"@@%i %i %i %i: %lg %lg %lg\n",i2,j,m,f,s,fabs(StateT[tn].cell_inf[j-1]),fabs(StateT[tn].cell_inf[j]));
*/							}
						f=(StateT[tn].cell_inf[j]<0);
						do
							{
							s=ranf_mt(ranf_seed, tn);
							l=c->InvCDF[(int) floor(s*1024)];
							while(c->cum_trans[l]<s) l++;
							ct=CellLookup[l];
							m=(int) (ranf_mt(ranf_seed, tn)*((double) ct->S0));
							i3=ct->susceptible[m];
							s2=dist2(Hosts+i3,Hosts+c->infected[j]);
							s=numKernel(s2)/c->max_trans[l];
							f2=0;
							if((ranf_mt(ranf_seed, tn)>=s)||(abs(Hosts[i3].inf)==5))
								{f2=1;}
							else if(m<ct->S)
								{
								if((!Hosts[i3].Travelling)&&((c!=ct)||(Hosts[i3].hh!=Hosts[c->infected[j]].hh)))
									{
									mi=Mcells+Hosts[c->infected[j]].mcell;
									mt=Mcells+Hosts[i3].mcell;
									s=CalcSpatialSusc(i3,ts,c->infected[j],tn); 
									if(bm)
										{
										if((dist2_raw(Households[Hosts[c->infected[j]].hh].loc_x,Households[Hosts[c->infected[j]].hh].loc_y,
											Households[Hosts[i3].hh].loc_x,Households[Hosts[i3].hh].loc_y)>P.MoveRestrRadius2))
											s*=P.MoveRestrEffect;
										}
									else if((mt->moverest!=mi->moverest)&&((mt->moverest==2)||(mi->moverest==2)))
										s*=P.MoveRestrEffect;
									if(!f)
										{
										for(m=f2=0;(m<P.PlaceTypeNum)&&(!f2);m++)
											if(Hosts[i3].PlaceLinks[m]>=0)
												f2=PLACE_CLOSED(m,Hosts[i3].PlaceLinks[m]);
										if(f2) {s*=P.PlaceCloseSpatialRelContact;}/* NumPCD++;} */
										f2=0;
										}
									if((s==1)||(ranf_mt(ranf_seed, tn)<s))
										{
										cq=((int) (ct-Cells))%P.NumThreads;
										if((Hosts[i3].infector==-1)&&(StateT[tn].n_queue[cq]<P.InfQueuePeakLength))
											{
											if((P.FalsePositiveRate>0)&&(ranf_mt(ranf_seed, tn)<P.FalsePositiveRate))
												Hosts[i3].infector=Hosts[i3].infect_type=-1;
											else
												{
												Hosts[i3].infector=c->infected[j];
												Hosts[i3].infect_type=2+2*NUM_PLACE_TYPES+INFECT_TYPE_MASK*(1+Hosts[c->infected[j]].infect_type/INFECT_TYPE_MASK);
												}
											StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++]=i3;
											}
										}
									}
								}
							}
						while(f2);
						}
					}
				}
			}
#pragma omp parallel for private(i,k) schedule(static,1)
		for(j=0;j<P.NumThreads;j++)
			{
			for(k=0;k<P.NumThreads;k++)
				{
				for(i=0;i<StateT[k].n_queue[j];i++)
					{
					if(Hosts[StateT[k].inf_queue[j][i]].infect_type==-1)
						DoFalseCase(StateT[k].inf_queue[j][i],t,ts,j);
					else
						DoInfect(StateT[k].inf_queue[j][i],t,j,run); //added int as argument to DoInfect to record run number: ggilani - 15/10/14
					}
				StateT[k].n_queue[j]=0;
				}
			}
}


void InfectSweep(double t, int run) //added run number as argument in order to record it in event log
{
	int i,j,k,l,m,n,i2,b,i3,f,f2,tn,cq,bm,ci;
	double seasonality,sbeta,hbeta,s,s2,s3,s4,s5,fp,contact_scale;
	cell *c,*ct;
	microcell *mi,*mt,*mp;
	unsigned short int ts; 
	person *si;
	
	if(!P.DoSeasonality)
		seasonality=1.0;
	else
		seasonality=P.Seasonality[((int) t)%DAYS_PER_YEAR];
	ts=(unsigned short int) (P.TimeStepsPerDay*t);
	fp=P.TimeStep/(1-P.FalsePositiveRate);
	sbeta=seasonality*fp*P.LocalBeta;
	hbeta=(P.DoHouseholds)?(seasonality*fp*P.HouseholdTrans):0;
	bm=((P.DoBlanketMoveRestr)&&(t>=P.MoveRestrTimeStart)&&(t<P.MoveRestrTimeStart+P.MoveRestrDuration));
#pragma omp parallel for private(j,k,l,m,n,i2,b,i3,f,f2,s,s2,s3,s4,s5,c,ct,mi,mt,mp,cq,ci,si,contact_scale) schedule(static,1)
	for(tn=0;tn<P.NumThreads;tn++)
		for(b=tn;b<P.NCP;b+=P.NumThreads)
			{	
			c=CellLookup[b];
			s5=0;
			for(j=0;j<c->I;j++)
				{
				ci=c->infected[j];
				si=Hosts+ci;
				
				if(hbeta>0)
					{
					if((Households[si->hh].nh>1)&&(!(si->nc_plus_hh_disabled & HH_DISABLED))&&(!si->Travelling))
						{
						l=Households[si->hh].FirstPerson;
						m=l+Households[si->hh].nh;
						s3=hbeta*CalcHouseInf(ci,ts);
						f=0;
						for(i3=l;(i3<m)&&(!f);i3++)
							for(i2=0;(i2<P.PlaceTypeNum)&&(!f);i2++)
								if(Hosts[i3].PlaceLinks[i2]>=0)
									{
									f=PLACE_CLOSED(i2,Hosts[i3].PlaceLinks[i2]);
									}
						if(f) {s3*=P.PlaceCloseHouseholdRelContact;}/* NumPCD++;}*/
						for(i3=l;i3<m;i3++)
							{
							if((Hosts[i3].inf==0)&& (!(Hosts[i3].nc_plus_hh_disabled & HH_DISABLED))&& (!Hosts[i3].Travelling))
								{
								s=s3*CalcHouseSusc(i3,ts,ci,tn);
								if(ranf_mt(ranf_seed, tn)<s)
									{
									cq=Hosts[i3].pcell%P.NumThreads;
									if((StateT[tn].n_queue[cq]<P.InfQueuePeakLength)) //(Hosts[i3].infector==-1)&&
										{
										if((P.FalsePositiveRate>0)&&(ranf_mt(ranf_seed, tn)<P.FalsePositiveRate))
											Hosts[i3].infector=Hosts[i3].infect_type=-1;
										else
											{
											Hosts[i3].infector=ci;
											Hosts[i3].infect_type=1+INFECT_TYPE_MASK*(1+si->infect_type/INFECT_TYPE_MASK);
											}
										StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++]=i3;
										}
									}
								}
							}
						}
					}
				if(P.DoPlaces)
					{
					if(!HOST_ABSENT(ci))
						{
						mi=Mcells+si->mcell;
						for(k=0;k<P.PlaceTypeNum;k++)
							{
							l=si->PlaceLinks[k];
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
							if((l>=0)&&(!PLACE_CLOSED(k,l)))
								{
								s3=fp*seasonality*contact_scale*CalcPlaceInf(ci,k,ts); //scale contacts so that hospitalised cases split their contacts between network and staff
								mp=Mcells+Places[k][l].mcell;
								if(bm)
									{
									if((dist2_raw(Households[si->hh].loc_x,Households[si->hh].loc_y,
										Places[k][l].loc_x,Places[k][l].loc_y)>P.MoveRestrRadius2))
											s3*=P.MoveRestrEffect;
									}
								else if((mi->moverest!=mp->moverest)&&((mi->moverest==2)||(mp->moverest==2)))
									s3*=P.MoveRestrEffect;

								if((k!=HOTEL_PLACE_TYPE)&&(!si->Travelling))
									{
									i2=(si->PlaceGroupLinks[k]);
									s4=s3;
									//s4=s3*(1-P.PlaceTypePropBetweenGroupLinks[k]*P.PlaceTypeGroupSizeParam1[k]/((double) Places[k][l].n));
									if(s4<0)
									{
										fprintf(stderr,"@@@ %lg\n",s4);
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
											n = (int)ignbin_mt(ranf_seed, (long)Places[k][l].nhcws, s4, tn);
										}
										else
										{
											n = (int)ignbin_mt(ranf_seed, (long)Places[k][l].group_size[i2], s4, tn);
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
									for(m=0;m<n;m++)
										{
										if (k == P.HospPlaceTypeNum)
										{	
											i3 = Places[k][l].members[SamplingQueue[tn][m]];
										}
										else
										{
											i3 = Places[k][l].members[Places[k][l].group_start[i2] + SamplingQueue[tn][m]];
										}
										if((Hosts[i3].inf==0)&&(!HOST_ABSENT(i3)))
											{
											mt=Mcells+Hosts[i3].mcell;
											ct=Cells+Hosts[i3].pcell;
											s=CalcPlaceSusc(i3,k,ts,ci,tn);
											if(bm)
												{
												if((dist2_raw(Households[Hosts[i3].hh].loc_x,Households[Hosts[i3].hh].loc_y,
													Places[k][l].loc_x,Places[k][l].loc_y)>P.MoveRestrRadius2))
													s*=P.MoveRestrEffect;
												}
											else if((mt->moverest!=mp->moverest)&&((mt->moverest==2)||(mp->moverest==2)))
												s*=P.MoveRestrEffect;
											if((s==1)||(ranf_mt(ranf_seed, tn)<s))
												{
												cq=Hosts[i3].pcell%P.NumThreads;
												if((StateT[tn].n_queue[cq]<P.InfQueuePeakLength)) //(Hosts[i3].infector==-1)&&
													{
													if((P.FalsePositiveRate>0)&&(ranf_mt(ranf_seed, tn)<P.FalsePositiveRate))
														Hosts[i3].infector=Hosts[i3].infect_type=-1;
													else
														{
														Hosts[i3].infector=ci;
														Hosts[i3].infect_type=2+k+INFECT_TYPE_MASK*(1+si->infect_type/INFECT_TYPE_MASK);
														}
													StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++]=i3;
													}
												}
											}
										}
									}
								if((k==HOTEL_PLACE_TYPE)||(!si->Travelling))
									{
									s3*=P.PlaceTypePropBetweenGroupLinks[k]*P.PlaceTypeGroupSizeParam1[k]/((double) Places[k][l].n);
									if(s3<0) 
									{
										fprintf(stderr,"@@@ %lg\n",s3);
										exit(1);
									}
									else if(s3>=1)
										n=Places[k][l].n;
									else
										n=(int) ignbin_mt(ranf_seed, (long) Places[k][l].n,s3,tn);
									if(n>0) SampleWithoutReplacement(tn,n,Places[k][l].n);		
									for(m=0;m<n;m++)
										{
										i3=Places[k][l].members[SamplingQueue[tn][m]];
										if((Hosts[i3].inf==0)&&(!HOST_ABSENT(i3)))
											{
											mt=Mcells+Hosts[i3].mcell;
											ct=Cells+Hosts[i3].pcell;
											s=CalcPlaceSusc(i3,k,ts,ci,tn); 
											if(bm)
												{
												if((dist2_raw(Households[Hosts[i3].hh].loc_x,Households[Hosts[i3].hh].loc_y,
													Places[k][l].loc_x,Places[k][l].loc_y)>P.MoveRestrRadius2))
													s*=P.MoveRestrEffect;
												}
											else if((mt->moverest!=mp->moverest)&&((mt->moverest==2)||(mp->moverest==2)))
												s*=P.MoveRestrEffect;
											if((s==1)||(ranf_mt(ranf_seed, tn)<s))
												{
												cq=Hosts[i3].pcell%P.NumThreads;
												if((StateT[tn].n_queue[cq]<P.InfQueuePeakLength))//(Hosts[i3].infector==-1)&&
													{
													if((P.FalsePositiveRate>0)&&(ranf_mt(ranf_seed, tn)<P.FalsePositiveRate))
														Hosts[i3].infector=Hosts[i3].infect_type=-1;
													else
														{
														Hosts[i3].infector=ci;
														Hosts[i3].infect_type=2+k+NUM_PLACE_TYPES+INFECT_TYPE_MASK*(1+si->infect_type/INFECT_TYPE_MASK);
														}
													StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++]=i3;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				if(sbeta>0)
					{
					if(si->Travelling)
						{s2=0;f=0;}
					else
						{
						s2=CalcSpatialInf(ci,ts);
						f=0;
						if(P.DoPlaces)
						  for(i3=0;(i3<P.PlaceTypeNum)&&(!f);i3++)
							if(si->PlaceLinks[i3]>=0)
								{
								f=PLACE_CLOSED(i3,si->PlaceLinks[i3]);
								}
						}
					if(f)
						{
						s2*=P.PlaceCloseSpatialRelContact;
						/* NumPCD++; */
						s5+=s2;
						StateT[tn].cell_inf[j]=-s5;
						}
					else
						{
						s5+=s2;
						StateT[tn].cell_inf[j]=s5;
						}
					}
				}
			if(s5>0)
				{
				n=(int) ignpoi_mt(ranf_seed, s5*sbeta*((double) c->tot_prob),tn);
				i2=c->I;
				if(n>0)
					{
					for(j=0;j<i2-1;j++) StateT[tn].cell_inf[j]/=s5;
					StateT[tn].cell_inf[i2-1]=(StateT[tn].cell_inf[i2-1]<0)?-1:1;
					}
				for(k=0;k<n;k++)
					{
					if(i2==1) 
						j=0;
					else
						{
						s=ranf_mt(ranf_seed, tn);
						j=m=i2/2;
						f=1;
						do
							{
							if(m>1) m/=2;
							if((j>0)&&(fabs(StateT[tn].cell_inf[j-1])>=s))
								{
								j-=m;
								if(j==0)
									f=0;
//									else
//										f=(fabs(StateT[tn].cell_inf[j-1])>=s);
								}
							else if((j<i2-1)&&(fabs(StateT[tn].cell_inf[j])<s))
								{
								j+=m;
								if(j==i2-1)
									f=0;
//									else
//										f=(fabs(StateT[tn].cell_inf[j])<s);
								}
							else
								{
/*									fprintf(stderr,"##%i %i %i %i: %lg %lg %lg\n",i2,j,m,f,s,fabs(StateT[tn].cell_inf[j-1]),fabs(StateT[tn].cell_inf[j]));
*/									f=0;
								}
							}
						while(f);
/*							if(fabs(StateT[tn].cell_inf[j])<s) fprintf(stderr,"##%i %i %i %i: %lg %lg %lg\n",i2,j,m,f,s,fabs(StateT[tn].cell_inf[j-1]),fabs(StateT[tn].cell_inf[j]));
						if((j>0)&&(fabs(StateT[tn].cell_inf[j-1])>=s)) fprintf(stderr,"@@%i %i %i %i: %lg %lg %lg\n",i2,j,m,f,s,fabs(StateT[tn].cell_inf[j-1]),fabs(StateT[tn].cell_inf[j]));
*/						}
					f=(StateT[tn].cell_inf[j]<0);
					ci=c->infected[j];
					si=Hosts+ci;
					do
						{
						s=ranf_mt(ranf_seed, tn);
						l=c->InvCDF[(int) floor(s*1024)];
						while(c->cum_trans[l]<s) l++;
						ct=CellLookup[l];
						m=(int) (ranf_mt(ranf_seed, tn)*((double) ct->S0));
						i3=ct->susceptible[m];
						s2=dist2(Hosts+i3,Hosts+ci);
						s=numKernel(s2)/c->max_trans[l];
						//alter acceptance probability for cross border effect here: testing the cross border effect - ggilani 19/01/15
						if((AdUnits[Mcells[Hosts[ci].mcell].adunit].id/P.CountryDivisor)!=(AdUnits[Mcells[Hosts[i3].mcell].adunit].id/P.CountryDivisor))
						{
							s*=P.PropCrossBorderInf;
						}
						f2=0;
						if((ranf_mt(ranf_seed, tn)>=s)||(abs(Hosts[i3].inf)==5))
							{f2=1;}
						//add cross border effect - relative contact assumed to be less if infector and infectee live in different countries: ggilani 09/12/14
						//else if(((AdUnits[Mcells[Hosts[ci].mcell].adunit].id/P.CountryDivisor)!=(AdUnits[Mcells[Hosts[i3].mcell].adunit].id/P.CountryDivisor))&&(ranf_mt(tn)>=P.PropCrossBorderInf)) //checking to see if they are in the same country
						//{
						//	f2=1;
						//}
						else if(m<ct->S)
							{
							if((!Hosts[i3].Travelling)&&((c!=ct)||(Hosts[i3].hh!=si->hh)))
								{
								mi=Mcells+si->mcell;
								mt=Mcells+Hosts[i3].mcell;
								s=CalcSpatialSusc(i3,ts,ci,tn); 
								if(bm)
									{
									if((dist2_raw(Households[si->hh].loc_x,Households[si->hh].loc_y,
										Households[Hosts[i3].hh].loc_x,Households[Hosts[i3].hh].loc_y)>P.MoveRestrRadius2))
										s*=P.MoveRestrEffect;
									}
								else if((mt->moverest!=mi->moverest)&&((mt->moverest==2)||(mi->moverest==2)))
									s*=P.MoveRestrEffect;
								if(!f)
									{
									for(m=f2=0;(m<P.PlaceTypeNum)&&(!f2);m++)
										if(Hosts[i3].PlaceLinks[m]>=0)
											{
											f2=PLACE_CLOSED(m,Hosts[i3].PlaceLinks[m]);
											}
									if(f2) {s*=P.PlaceCloseSpatialRelContact;}/* NumPCD++;} */
									f2=0;
									}
								if((s==1)||(ranf_mt(ranf_seed, tn)<s))
									{
									cq=((int) (ct-Cells))%P.NumThreads;
									if((Hosts[i3].inf==0)&&(StateT[tn].n_queue[cq]<P.InfQueuePeakLength)) //Hosts[i3].infector==-1
										{
										if((P.FalsePositiveRate>0)&&(ranf_mt(ranf_seed, tn)<P.FalsePositiveRate))
											Hosts[i3].infector=Hosts[i3].infect_type=-1;
										else
											{
											Hosts[i3].infector=ci;
											Hosts[i3].infect_type=2+2*NUM_PLACE_TYPES+INFECT_TYPE_MASK*(1+si->infect_type/INFECT_TYPE_MASK);
											}
										StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++]=i3;
										}
									}
								}
							}
						}
					while(f2);
					}
				}
			}

	
#pragma omp parallel for private(i,k) schedule(static,1)
		for(j=0;j<P.NumThreads;j++)
			{
			for(k=0;k<P.NumThreads;k++)
				{
				for(i=0;i<StateT[k].n_queue[j];i++)
					{
					if(Hosts[StateT[k].inf_queue[j][i]].infect_type==-1)
						DoFalseCase(StateT[k].inf_queue[j][i],t,ts,j);
					else
						DoInfect(StateT[k].inf_queue[j][i],t,j, run);
					}
				StateT[k].n_queue[j]=0;
				}
			}
}

void IncubRecoverySweep(double t,int run)
{
	int i,j,k,l,b,tn,ci,day;
	double currInfSafeFuneral;
	cell *c;
	person *si;
	unsigned short int ts,tc;
	double currPropSafeFunerals,currRelInfFuneral;

	ts=(unsigned short int) (P.TimeStepsPerDay*t);

	//if we're doing hospitalisation by admin unit, set number of people in discharge and admission queues to zero; otherwise if doing hospitalisation by place, only set state hospitalisation queue to zero
	if(P.DoHospitalisation)
	{
		if(P.DoETUByAdUnit)
		{
			for(i=0;i<P.NumAdunits;i++)
			{
				AdUnits[i].nh_queue=0;
				for(j=0;j<P.NumThreads;j++)
				{
					StateT[j].nh_queue[i]=0;
					StateT[j].nhd_queue[i]=0;
				}
			}
		}
		else
		//If we're doing hospitalisation by place, we'll still use the same queue to store individuals who need hospitalisation, but we'll store them all in the same queue, not by admin unit
		{
			for(j=0;j<P.NumThreads;j++)
				{
					StateT[j].nhd_queue[0]=0;
				}
		}
	}

#pragma omp parallel for private(j,k,l,b,c,tn,tc,ci,si,day) schedule(static,1)
	for(tn=0;tn<P.NumThreads;tn++)
	{
		for(b=tn;b<P.NCP;b+=P.NumThreads)
		{
			c=CellLookup[b];
			for(j=((int) c->L-1);j>=0;j--)
				if(ts>=Hosts[c->latent[j]].latent_time) DoIncub(c->latent[j],ts,tn,run);
			StateT[tn].n_queue[0]=0;
			for(j=c->I-1;j>=0;j--)
			{
				ci=c->infected[j];
				si=Hosts+ci;
				tc=si->latent_time+((int) (P.LatentToSymptDelay/P.TimeStep));
				if((P.DoSymptoms)&&(ts==tc))
					DoCase(ci,t,ts,tn);
				if ((ts==si->detect_time)&&(si->detected))
				{
					DoDetectedCase(ci, t, ts, tn);
				}

				//Now considered which of infected have reached time to hospitalisation
				if ((P.DoHospitalisation) && ((ts >= si->hospital_time) && (ts < (si->hospital_time + (int)(P.HospWaitingTime * P.TimeStepsPerDay)))) && (abs(si->inf) != 6) && (!si->hospitalised) &&(!si->etu) && (si->detected))
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
				if(ts>=si->recovery_time)
				{
					if((P.DoHospitalisation)&&(P.DoETUByAdUnit)&&((si->hospitalised)||(si->etu)))
					{
						//if someone has reached their recovery time, regardless of whether it's a death or recovery, they will leave hospital at this point if they are hospitalised
						//mark someone to be discharged; however, if we're doing hospitalisation by place, this gets taken care of in HospitalSweep
						StateT[tn].hd_queue[Mcells[si->mcell].adunit][StateT[tn].nhd_queue[Mcells[si->mcell].adunit]++]=ci;
					}
					if((si->to_die)&&(P.DoFuneralTransmission)&&(abs(si->inf)!=6)) //only go through this if block if using waiting time until leaving infectious class, and if host is alive and infectious
					{
						//set infection status flag to dead and infectious: -6
						si->inf=-6;
						//do some accounting here to make sure we count the death at the right time
						StateT[tn].cumD++;
						StateT[tn].cumDa[HOST_AGE_GROUP(ci)]++;
						if (P.DoAdUnits) StateT[tn].cumD_adunit[Mcells[si->mcell].adunit]++;

						// set recovery time to current recovery time plus length of funeral transmission duration
						si->recovery_time+=(unsigned short int)(P.FuneralTransmissionDuration*P.TimeStepsPerDay);
						
						if((P.DoHospitalisation)&&((si->hospitalised)||(si->etu))) //this is regardless of whether we are considering hospitalisation by admin unit or place
						{
							//if in hospital, they are definitely detected when they die
							StateT[tn].cumDD++;
							if (P.DoAdUnits) StateT[tn].cumDD_adunit[Mcells[si->mcell].adunit]++;
							if ((t >= P.FuneralControlTimeStart)&&(State.cumSDB_adunit[Mcells[si->mcell].adunit]<AdUnits[Mcells[si->mcell].adunit].maxSDB))
							{
								//if someone has died in hospital, we assumed that they will have a safe burial
								si->infectiousMult = (P.RelativeInfectiousnessFuneral * P.RelInfSafeFuneral);
								//if in hospital, they definitely have a safe burial
								StateT[tn].cumSDB++;
								if (P.DoAdUnits) StateT[tn].cumSDB_adunit[Mcells[si->mcell].adunit]++;
							}
						}
						else if(si->detected)
						{
							//if a detected case, they are definitely detected when they die
							StateT[tn].cumDD++;
							if (P.DoAdUnits) StateT[tn].cumDD_adunit[Mcells[si->mcell].adunit]++;
							//alter host's infectiousness, taking into account relative reduction in infectiousness due to safe burial
							if ((t >= P.FuneralControlTimeStart) && (ranf_mt(ranf_seed, tn) <= P.ProportionSafeFuneral)&& (State.cumSDB_adunit[Mcells[si->mcell].adunit] < AdUnits[Mcells[si->mcell].adunit].maxSDB))
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
							si->infectiousMult=P.RelativeInfectiousnessFuneral;
							si->safeBurial = 0;
						}
						//}

					}
					if((!si->to_die)&&(ts>=si->recovery_time))
						StateT[tn].inf_queue[0][StateT[tn].n_queue[0]++]=ci;
					else if((si->to_die)&&(ts>=si->recovery_time))
					{
						if((!HOST_TREATED(ci))||(abs(si->inf)==6)) // if infectious status is 6, then host is already dead and cannot be treated!
							DoDeath(ci,tn,run);
						else if(ranf_mt(ranf_seed, tn)>=P.EvolResistRelTreatDeathDrop[si->resist])
							StateT[tn].inf_queue[0][StateT[tn].n_queue[0]++]=ci;
						else
							DoDeath(ci,tn,run);
					}

				}
				
			}
				
			for(j=0;j<StateT[tn].n_queue[0];j++) DoRecover(StateT[tn].inf_queue[0][j],run,tn);
			StateT[tn].n_queue[0]=0;
			l=(P.EvolInfectMutationRate==0)?0:((int) ignbin_mt(ranf_seed, (long) c->I,P.TimeStep*P.EvolInfectMutationRate,tn));
			for(j=0;j<l;j++)
			{
				k=(int) (((double) c->I)*ranf_mt(ranf_seed, tn));
				if(Hosts[c->infected[k]].base_inf_level<P.EvolInfectMax)
					Hosts[c->infected[k]].base_inf_level+=P.EvolInfectStep;
			}
		}
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
	
	if((P.DoHospitalisation)&(P.IncludeHospitalPlaceType))
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
							if (ranf_mt(ranf_seed, tn)<P.propContactLost) //if this contact is going to be lost to follow up
							{
								Hosts[AdUnits[i].ct_queue[j]].contactTraced_end_time = ts + (unsigned short int) ((P.contactTraceDuration*ranf_mt(ranf_seed, tn)) * P.TimeStepsPerDay);
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

	for(i=0;i<P.NumAdunits;i++)
	{
		AdUnits[i].nct_queue=0;
		for(j=0;j<P.NumThreads;j++)
		{
			StateT[j].nct_queue[i]=0;
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
					DoVacc(State.vacc_queue[i], ts,1);
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
	int i,j;//,nHospBeds,pop,newBeds;
	int nCases, numBedsTotal, numBedsInUse;
	double capacityFlag;
	//pop=P.N;
	
	//first update current mean time from onset to hospitalisation
	for(i=P.CurrIndMeanTimeToHosp;i<P.NMeanTimeToHosp;i++)
	{
		if(t>=P.ChangePointMeanTimeToHosp[i]+P.ETUTimeStart)
		{
			P.HospitalisationTime=P.MeanTimeToHosp[i];
			P.CurrIndMeanTimeToHosp++;
		}
	}

	//first update current mean time from onset to hospitalisation
	if (P.DoContactTracing)
	{
		for (i = P.CurrIndMeanTimeToHospCT; i < P.NMeanTimeToHospCT; i++)
		{
			if (t >= P.ChangePointMeanTimeToHospCT[i]+P.ETUTimeStart)
			{
				P.HospitalisationTime_contactTrace = P.MeanTimeToHospCT[i];
				P.CurrIndMeanTimeToHospCT++;
			}
		}
	}

	if ((P.DoReactETUBeds)&&(t>=P.ETUTimeStart))
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
				else if((AdUnits[i].ETUbedsActive)&&(AdUnits[i].nextTimeToETUBeds<t))
				{
					//calculate current capacity
					numBedsTotal= BedsAvailablePerAdUnit(t, i);
					numBedsInUse= AdUnits[i].currentETUBeds;
					capacityFlag = (double)(numBedsInUse) / (double)(numBedsTotal);
					//if this has changed (more specifically, if it has got bigger because we've crossed the threshold for adding beds again)
					if (capacityFlag>P.CapacityToMoreETUBeds) //we also can't add more beds until we've added the last set
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
				AdUnits[i].totalETUBeds +=AdUnits[i].nextETUBeds;
				State.NumBeds+= AdUnits[i].nextETUBeds;
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
			if ((t >= AdUnits[i].timeETUBedsAvailable)&&(AdUnits[i].ETUbedsActive==0))
			{
				AdUnits[i].totalETUBeds = AdUnits[i].initialETUBeds;
				AdUnits[i].ETUbedsActive = 1;
			}
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
			if ((AdUnits[i].contactTraceThresholdCrossed)&&(AdUnits[i].nextTimeToCT < t))
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

	if  (t >= P.FuneralControlTimeStart)
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
void UpdateVaccination(double t,int n)
{

	int i,weeklyDC;//,nHospBeds,pop,newBeds;
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

	if ((t >= P.TimeToIncVaccRing)&&(P.NVaccRingsActive<3))
	{
		P.NVaccRingsActive++;
	}

	//update vacc dose per day
	//first find total number of detected cases over the past 7 days
	if (P.UpdateVaccDosePerDay)
	{
		weeklyDC = 0;
		for (i = max(n - 6,0); i <= n; i++)
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
			P.MaxVaccDosePerDay = P.MaxVaccGeoDosePerDay;
			P.UpdateVaccDosePerDay = 0;
		}
	}


}


void UpdateCaseDetection(double t)
{
	int i,j;
	
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
	else if ((P.DoUpdateCaseDetectionByCases)&&(P.UpdateCaseDetectionByCasesFlag==0))
	{
		if (State.cumDC>P.CaseThresholdUntilUpdateCaseDetection)
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

void DemogSweep(double t)
{
	int j,b,tn,f;
	cell *c;

	f=0;
#pragma omp parallel for private(j,b,c,tn) reduction(+:f) schedule(static,1)
	for(tn=0;tn<P.NumThreads;tn++)
		for(b=tn;b<P.NCP;b+=P.NumThreads)
			{
			c=CellLookup[b];
			for(j=c->S-1;j>=0;j--)
				UpdateSuscPersonDemogSIR(c->susceptible[j],tn);
			for(j=c->S+c->L+c->I;j<c->n;j++)
				UpdatePersonDemogSIR(c->susceptible[j],tn);
			if(c->S>c->S0) f++;
			}
	if(f>0) {fprintf(stderr,"** UpdateProbs - %i\n",f);UpdateProbs(0);}
}


void SIASweep(unsigned short int ts)
{
#ifdef NEW_AGE_MODEL
#ifdef FRESSCA
	int f,i2,i,j,k,l,m,a,tn;

	//fprintf(stderr,"\nApplying SIA\n");
#pragma omp parallel for private(f,i2,i,j,k,l,m,a,tn) schedule(static,1)
	for(tn=0;tn<P.NumThreads;tn++)
		for(m=tn;m<P.NMCP;m+=P.NumThreads)
			{
			i=McellLookup[m]-Mcells;
			if(((P.SIADoAllCountries)||(Mcells[i].country==P.DistribNetCountry))&&(P.SIAEffectiveCoverage>0))
			  for(i2=0;i2<Mcells[i].n;i2++)
				{
				l=Mcells[i].members[i2];
				a=((int) ts)-Hosts[l].birth_time;
				if((a>=P.usSIAMinAge)&&(a<P.usSIAMaxAge))
					{
					if(P.SIAEffectiveCoverage==1)
						f=1;
					else
						f=(ranf_mt(ranf_seed, tn)<P.SIAEffectiveCoverage);
					if(f)
						{
						j=Hosts[l].vacc_queue_pos;
						if(j>=0)
							{
							Mcells[i].dvacc_queue[j]=-1;
							Hosts[l].vacc_queue_pos=-1;
							}
						// Add l to vaccination queue
						j=Mcells[i].dvacc_count;
						k=(Mcells[i].ndvacc_queue+j)%Mcells[i].n;
						Mcells[i].dvacc_queue[k]=l;
						Mcells[i].dvacc_min_vacc_time[k]=ts;
						Mcells[i].dvacc_expiry_time[k]=ts+P.usSIADuration;
						Mcells[i].dvacc_count++;
						Hosts[l].vacc_queue_pos=k;
						}
					}
				}
			}
#endif
#endif
}


int TreatSweep(double t)
	{
	int i,j,i2,j2,k,l,m,b,rd,prd,mrd,f,f1,f2,f3,f4,tn,bs,ad,ad2,ii,jj,npropvacc;
	int minx, maxx, miny, trig_thresh, nckwp,FirstPerson,LastPerson;
	double maxdist,maxdist2,dosepercase, dosepercell,logpop;
	unsigned short int ts,tstf,tstb,tsvb,tspf,tsmb,tsmf,tssdf,tskwpf;
	int global_trig;
	double r,rad2;

	ts=(unsigned short int) (P.TimeStepsPerDay*t);
	rd=(int) ceil(P.TreatRadius/P.mcwidth);
	prd=(int) ceil(P.PlaceCloseRadius/P.mcwidth);
	mrd=(int) ceil(P.MoveRestrRadius/P.mcwidth);
	f=f1=0;
	if(P.DoGlobalTriggers)
		{
		if(P.DoPerCapitaTriggers)
			global_trig=(int) floor(((double) State.trigDC)*P.GlobalIncThreshPop/((double) P.N));
		else
			global_trig=State.trigDC;
		}
	else
		global_trig=0;
	if((P.DoPlaces)&&(t>=P.TreatTimeStart)&&(t<P.TreatTimeStart+P.TreatPlaceGeogDuration)&&(State.cumT<P.TreatMaxCourses))
		{
		tstf=(unsigned short int) (P.TimeStepsPerDay*(t+P.TreatDelayMean+P.TreatProphCourseLength));
#pragma omp parallel for private(i,j,k,l,m,f) reduction(+:f1) schedule(static,1)
		for(i=0;i<P.NumThreads;i++)
			{
			for(j=0;j<P.PlaceTypeNum;j++)
				{
				for(k=0;k<StateT[i].np_queue[j];k++)
					{
					l=StateT[i].p_queue[j][k];
					if(P.DoPlaceGroupTreat)
						{
						f=StateT[i].pg_queue[j][k];
						for(m=((int) Places[j][l].group_start[f]);m<((int) (Places[j][l].group_start[f]+Places[j][l].group_size[f]));m++)
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
*/							if((!HOST_TO_BE_TREATED(Places[j][l].members[m]))&&((P.TreatPlaceTotalProp[j]==1)||(ranf_mt(ranf_seed, i)<P.TreatPlaceTotalProp[j])))
								DoProph(Places[j][l].members[m],ts,i);
							}
						}
					else
						{
						if((Places[j][l].treat)&&(!PLACE_TREATED(j,l)))
							{
							f1=1;
							Places[j][l].treat_end_time=tstf;
							for(m=0;m<Places[j][l].n;m++)
								if(!HOST_TO_BE_TREATED(Places[j][l].members[m]))
									{
									if((P.TreatPlaceTotalProp[j]==1)||(ranf_mt(ranf_seed, i)<P.TreatPlaceTotalProp[j]))
										DoProph(Places[j][l].members[m],ts,i);
									}
							}
						Places[j][l].treat=0;
						}
					}
				StateT[i].np_queue[j]=0;
				}
			}
		}
	if((P.DoMassVacc)&&(t>=P.VaccTimeStart))
		for(j=0;j<2;j++)
			{
			m=P.VaccMaxCourses;
			if(m>State.n_mvacc) m=State.n_mvacc;
#pragma omp parallel for private(i) schedule(static,1000)
			for(i=State.mvacc_cum;i<m;i++)
				DoVacc(State.mvacc_queue[i],ts,0);
			State.mvacc_cum=m;
			}
	if((t>=P.TreatTimeStart)||(t>=P.VaccTimeStart)||(t>=P.PlaceCloseTimeStart)||(t>=P.MoveRestrTimeStart)||(t>=P.SocDistTimeStart)||(t>=P.KeyWorkerProphTimeStart))
		{
		tstf=(unsigned short int) (P.TimeStepsPerDay*(t+P.TreatProphCourseLength)-1);
		tstb=(unsigned short int) (P.TimeStepsPerDay*(t+P.TreatDelayMean));
		//tsvb= (unsigned short int) (P.TimeStepsPerDay * t)+1; //modified this so that the delay gets taken care of in the vaccine queue and ring and geo vaccinations align timewise
		tsvb=(unsigned short int) (P.TimeStepsPerDay*(t+P.VaccDelayMean));
		tspf=(unsigned short int) ceil(P.TimeStepsPerDay*(t+P.PlaceCloseDelayMean+P.PlaceCloseDuration));
		tsmf=(unsigned short int) ceil(P.TimeStepsPerDay*(t+P.MoveRestrDuration));
		tsmb=(unsigned short int) floor(P.TimeStepsPerDay*(t+P.MoveDelayMean));
		tssdf=(unsigned short int) ceil(P.TimeStepsPerDay*(t+P.SocDistDuration));
		tskwpf=(unsigned short int) ceil(P.TimeStepsPerDay*(t+P.KeyWorkerProphRenewalDuration));
		nckwp=(int) ceil(P.KeyWorkerProphDuration/P.TreatProphCourseLength);
#pragma omp parallel for private(tn,i,i2,j,j2,k,l,m,b,bs,minx,maxx,miny,f2,f3,f4,trig_thresh,r,ad,ad2,maxdist,maxdist2,dosepercase,logpop,ii,jj,npropvacc,rad2,FirstPerson,LastPerson) reduction(+:f) schedule(static,1) //,ii,ll,nvacc,maxvacc,maxdist,maxdist2,dosepercase,logpop
		for(tn=0;tn<P.NumThreads;tn++)
			for(bs=tn;bs<P.NMCP;bs+=P.NumThreads)
				{
				b=(int) (McellLookup[bs]-Mcells);
				ad=(P.DoAdUnits)?AdUnits[Mcells[b].adunit].id:0;
				if((!P.RestrictTreatToTarget)||(Mcells[b].country==P.TargetCountry))
					{
					if((Mcells[b].treat==2)&&(ts>=Mcells[b].treat_end_time))
						{
						f=1;
						Mcells[b].treat=0;
						}
					if((Mcells[b].treat==1)&&(ts>=Mcells[b].treat_start_time))
						{
						f=1;
						Mcells[b].treat=2;
						Mcells[b].treat_trig=0;
						Mcells[b].treat_end_time=tstf;
						for(i=0;i<Mcells[b].n;i++)
							{
							l=Mcells[b].members[i];
							if((!HOST_TO_BE_TREATED(l))&&((P.TreatPropRadial==1)||(ranf_mt(ranf_seed, tn)<P.TreatPropRadial)))
								DoProphNoDelay(l,ts,tn,1);
							}
						}
					
					trig_thresh=(P.DoPerCapitaTriggers)?((int) ceil(((double) (Mcells[b].n*P.TreatCellIncThresh))/P.IncThreshPop)):P.TreatCellIncThresh;
					if((t>=P.TreatTimeStart)&&(Mcells[b].treat==0)&&(((Mcells[b].treat_trig>=trig_thresh)&&(!P.DoGlobalTriggers))||(global_trig>=P.TreatCellIncThresh))&&(P.TreatRadius2>0))
						{
						minx=(b/P.nmch);miny=(b%P.nmch);
						k=b;
						maxx=0;
						i=j=m=f2=0;
						l=f3=1;
						if((!P.TreatByAdminUnit)||(ad>0))
							{
							ad2=ad/P.TreatAdminUnitDivisor;
							do
								{
								if((minx>=0)&&(minx<P.nmcw)&&(miny>=0)&&(miny<P.nmch))
									{
									if(P.TreatByAdminUnit)
										f4=(AdUnits[Mcells[k].adunit].id/P.TreatAdminUnitDivisor==ad2);
									else
										f4=((r=dist2_mm(Mcells+b,Mcells+k))<P.TreatRadius2);
									if(f4)
										{
										f=f2=1;
										if((Mcells[k].n>0)&&(Mcells[k].treat==0)&&((!P.RestrictTreatToTarget)||(Mcells[k].country==P.TargetCountry)))
											{
											Mcells[k].treat_start_time=tstb;
											Mcells[k].treat=1;
											maxx+=Mcells[k].n;
											}
										}
									}
								if(j==0)
									minx=minx+1;
								else if(j==1)
									miny=miny-1;
								else if(j==2)
									minx=minx-1;
								else if(j==3)
									miny=miny+1;
								m=(m+1)%l;
								if(m==0)
									{
									j=(j+1)%4;
									i=(i+1)%2;
									if(i==0) l++;
									if(j==1) {f3=f2;f2=0;}
									}
								k=((minx+P.nmcw)%P.nmcw)*P.nmch+(miny+P.nmch)%P.nmch;
								}
							while((f3)&&(maxx<P.TreatMaxCoursesPerCase));
							}
						}
					//if Mcells.vacc==2, it is being vaccinated as a trigger cell, if Mcells.vacc==1, it is being vaccinated as a peripheral cell, but can still trigger vaccination, if Mcells.vacc==0, it
					//is currently not scheduled to be vaccinated
					
					if ((P.LimitGeoVaccDosesPerCase)&&(P.StopVaccinationPostThreshold))
					{
						npropvacc=Mcells[b].n*(0.9*P.VaccProp);
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
					trig_thresh=(P.DoPerCapitaTriggers)?((int) ceil(((double) (Mcells[b].n*P.VaccCellIncThresh))/P.IncThreshPop)):P.VaccCellIncThresh;
					if (P.LimitGeoVaccDosesPerCase)
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
						if ((!P.DoMassVacc) && (P.VaccRadius2 > 0) && (t >= P.VaccTimeStart) && (Mcells[b].vacc==0) && (((Mcells[b].vacc_trig >= trig_thresh) && (!P.DoGlobalTriggers) && (State.cumV < P.VaccMaxCourses)) || (global_trig >= P.VaccCellIncThresh))) //changed from VaccTimeStart to VaccTimeStarGeo
						{
							minx = (b / P.nmch); miny = (b % P.nmch);
							k = b;
							//number of vaccine doses to distribute geographically in this cell
							//nvacc=maxvacc=P.VaccDosesPerCasePerCell * Mcells[b].vacc_trig;
							//max distance between cells;
							maxx=maxdist = maxdist2 = dosepercase = dosepercell = 0;
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
					if((Mcells[b].placeclose==2)&&(ts>=Mcells[b].place_end_time))
						{
						f=1;
						Mcells[b].placeclose=P.DoPlaceCloseOnceOnly;
						}
					trig_thresh=(P.DoPerCapitaTriggers)?((int) ceil(((double) (Mcells[b].n*P.PlaceCloseCellIncThresh))/P.IncThreshPop)):P.PlaceCloseCellIncThresh;
					if((P.DoPlaces)&&(t>=P.PlaceCloseTimeStart)&&(Mcells[b].placeclose==0))
						{
						if(((P.PlaceCloseByAdminUnit)&&(AdUnits[Mcells[b].adunit].place_close_trig<USHRT_MAX-1)&&(((double) AdUnits[Mcells[b].adunit].place_close_trig)/((double) AdUnits[Mcells[b].adunit].NP)>P.PlaceCloseAdunitPropThresh))
						||((!P.PlaceCloseByAdminUnit)&&(((Mcells[b].place_trig>=trig_thresh)&&(!P.DoGlobalTriggers))||(global_trig>=P.PlaceCloseCellIncThresh))))
							{
//							if(P.PlaceCloseByAdminUnit) AdUnits[Mcells[b].adunit].place_close_trig=USHRT_MAX-1; // This means schools only close once
							AdUnits[Mcells[b].adunit].place_close_trig=0; //Schools can close multiple times
							minx=(b/P.nmch);miny=(b%P.nmch);
							k=b;
							i=j=m=f2=0;
							l=f3=1;
							if((!P.PlaceCloseByAdminUnit)||(ad>0))
								{
								ad2=ad/P.PlaceCloseAdminUnitDivisor;
								do
									{
									if((minx>=0)&&(minx<P.nmcw)&&(miny>=0)&&(miny<P.nmch))
										{
										if(P.PlaceCloseByAdminUnit)
											f4=(AdUnits[Mcells[k].adunit].id/P.PlaceCloseAdminUnitDivisor==ad2);
										else
											f4=((r=dist2_mm(Mcells+b,Mcells+k))<P.PlaceCloseRadius2);
										if(f4)
											{
											f2=f=1;
											if((Mcells[k].n>0)&&(Mcells[k].placeclose==0)&&((!P.RestrictTreatToTarget)||(Mcells[k].country==P.TargetCountry)))
												{
												Mcells[k].place_end_time=tspf;
												Mcells[k].place_trig=0;
												Mcells[k].placeclose=2;
												for(j2=0;j2<P.PlaceTypeNum;j2++)
													if(j2!=HOTEL_PLACE_TYPE)
														for(i2=0;i2<Mcells[k].np[j2];i2++)
															DoPlaceClose(j2,Mcells[k].places[j2][i2],ts,tn,1);
												}
											}
										}
									if(j==0)
										minx=minx+1;
									else if(j==1)
										miny=miny-1;
									else if(j==2)
										minx=minx-1;
									else if(j==3)
										miny=miny+1;
									m=(m+1)%l;
									if(m==0)
										{
										j=(j+1)%4;
										i=(i+1)%2;
										if(i==0) l++;
										if(j==1) {f3=f2;f2=0;}
										}
									k=((minx+P.nmcw)%P.nmcw)*P.nmch+(miny+P.nmch)%P.nmch;
									if((k<0)||(k>P.nmcw*P.nmch)) fprintf(stderr,"### %i %i %i\n ",k,minx,miny);
									}
								while(f3);
								}
							}
						}
					if((Mcells[b].moverest==2)&&(ts>=Mcells[b].move_end_time))
						{
						f=1;
						Mcells[b].moverest=0;
						}
					if((Mcells[b].moverest==1)&&(ts>=Mcells[b].move_start_time))
						{
						f=1;
						Mcells[b].moverest=2;
						Mcells[b].move_trig=0;
						Mcells[b].move_end_time=tsmf;
						}
					trig_thresh=(P.DoPerCapitaTriggers)?((int) ceil(((double) (Mcells[b].n*P.MoveRestrCellIncThresh))/P.IncThreshPop)):P.MoveRestrCellIncThresh;
					if((t>=P.MoveRestrTimeStart)&&(Mcells[b].moverest==0)&&(((Mcells[b].move_trig>=trig_thresh)&&(!P.DoGlobalTriggers))||(global_trig>=P.MoveRestrCellIncThresh)))
						{
						minx=(b/P.nmch);miny=(b%P.nmch);
						k=b;
						i=j=m=f2=0;
						l=f3=1;
						if((!P.MoveRestrByAdminUnit)||(ad>0))
							{
							ad2=ad/P.MoveRestrAdminUnitDivisor;
							do
								{
								if((minx>=0)&&(minx<P.nmcw)&&(miny>=0)&&(miny<P.nmch))
									{
									if(P.MoveRestrByAdminUnit)
										f4=(AdUnits[Mcells[k].adunit].id/P.MoveRestrAdminUnitDivisor==ad2);
									else
										f4=((r=dist2_mm(Mcells+b,Mcells+k))<P.MoveRestrRadius2);
									if(f4)
										{
										f=f2=1;
										if((Mcells[k].n>0)&&(Mcells[k].moverest==0)&&((!P.RestrictTreatToTarget)||(Mcells[k].country==P.TargetCountry)))
											{
											Mcells[k].move_start_time=tsmb;
											Mcells[k].moverest=1;
											}
										}
									}
								if(j==0)
									minx=minx+1;
								else if(j==1)
									miny=miny-1;
								else if(j==2)
									minx=minx-1;
								else if(j==3)
									miny=miny+1;
								m=(m+1)%l;
								if(m==0)
									{
									j=(j+1)%4;
									i=(i+1)%2;
									if(i==0) l++;
									if(j==1) {f3=f2;f2=0;}
									}
								k=((minx+P.nmcw)%P.nmcw)*P.nmch+(miny+P.nmch)%P.nmch;
								}
							while(f3);
							}
						}					

					if((Mcells[b].socdist==2)&&(ts>=Mcells[b].socdist_end_time))
						{
						f=1;
						Mcells[b].socdist=0;
						}
					trig_thresh=(P.DoPerCapitaTriggers)?((int) ceil(((double) (Mcells[b].n*P.SocDistCellIncThresh))/P.IncThreshPop)):P.SocDistCellIncThresh;
					if((t>=P.SocDistTimeStart)&&(Mcells[b].socdist==0)&&(((Mcells[b].socdist_trig>=trig_thresh)&&(!P.DoGlobalTriggers))||(global_trig>=P.SocDistCellIncThresh)))
						{
						minx=(b/P.nmch);miny=(b%P.nmch);
						k=b;
						i=j=m=f2=0;
						l=f3=1;
						do
							{
							if((minx>=0)&&(minx<P.nmcw)&&(miny>=0)&&(miny<P.nmch))
								if(dist2_mm(Mcells+b,Mcells+k)<P.SocDistRadius2)
									{
									f=f2=1;
									if((Mcells[k].n>0)&&(Mcells[k].socdist==0)&&((!P.RestrictTreatToTarget)||(Mcells[k].country==P.TargetCountry)))
										{
										Mcells[k].socdist=2;
										Mcells[k].socdist_trig=0;
										Mcells[k].socdist_end_time=tssdf;
										}
									}
								if(j==0)
									minx=minx+1;
								else if(j==1)
									miny=miny-1;
								else if(j==2)
									minx=minx-1;
								else if(j==3)
									miny=miny+1;
								m=(m+1)%l;
								if(m==0)
									{
									j=(j+1)%4;
									i=(i+1)%2;
									if(i==0) l++;
									if(j==1) {f3=f2;f2=0;}
									}
								k=((minx+P.nmcw)%P.nmcw)*P.nmch+(miny+P.nmch)%P.nmch;
							}
						while(f3);
						}					

					if((Mcells[b].keyworkerproph==2)&&(ts>=Mcells[b].keyworkerproph_end_time))
						{
						f=1;
						Mcells[b].keyworkerproph=0;
						}
					trig_thresh=(P.DoPerCapitaTriggers)?((int) ceil(((double) (Mcells[b].n*P.KeyWorkerProphCellIncThresh))/P.IncThreshPop)):P.KeyWorkerProphCellIncThresh;
					if((P.DoPlaces)&&(t>=P.KeyWorkerProphTimeStart)&&(Mcells[b].keyworkerproph==0)&&(((Mcells[b].keyworkerproph_trig>=trig_thresh)&&(!P.DoGlobalTriggers))||(global_trig>=P.KeyWorkerProphCellIncThresh)))
						{
						minx=(b/P.nmch);miny=(b%P.nmch);
						k=b;
						i=j=m=f2=0;
						l=f3=1;
						do
							{
							if((minx>=0)&&(minx<P.nmcw)&&(miny>=0)&&(miny<P.nmch))
								if(dist2_mm(Mcells+b,Mcells+k)<P.KeyWorkerProphRadius2)
									{
									f=f2=1;
									if((Mcells[k].n>0)&&(Mcells[k].keyworkerproph==0)&&((!P.RestrictTreatToTarget)||(Mcells[k].country==P.TargetCountry)))
										{
										Mcells[k].keyworkerproph=2;
										Mcells[k].keyworkerproph_trig=0;
										Mcells[k].keyworkerproph_end_time=tskwpf;
										for(i2=0;i2<Mcells[k].n;i2++)
											{
											j2=Mcells[k].members[i2];
											if((Hosts[j2].keyworker)&&(!HOST_TO_BE_TREATED(j2)))
												DoProphNoDelay(j2,ts,tn,nckwp);
											}
										}
									}
								if(j==0)
									minx=minx+1;
								else if(j==1)
									miny=miny-1;
								else if(j==2)
									minx=minx-1;
								else if(j==3)
									miny=miny+1;
								m=(m+1)%l;
								if(m==0)
									{
									j=(j+1)%4;
									i=(i+1)%2;
									if(i==0) l++;
									if(j==1) {f3=f2;f2=0;}
									}
								k=((minx+P.nmcw)%P.nmcw)*P.nmch+(miny+P.nmch)%P.nmch;
							}
						while(f3);
						}					

					}
				}
			for(i=0;i<P.NumThreads;i++)
				{
				State.cumT+=StateT[i].cumT;
				State.cumTP+=StateT[i].cumTP;
				State.cumUT+=StateT[i].cumUT;
				//State.cumV+=StateT[i].cumV;
				StateT[i].cumT=StateT[i].cumUT=StateT[i].cumTP=0; //StateT[i].cumV=0;
				}
		}
	f+=f1;


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
	


	return (f>0);
}


#ifdef FRESSCA
void DistrVaccSweep(double t)n
{
	unsigned short int ts;
	microcell * m2vacc;
	int i,j,k,f,n,imc,ids,nq,nv,nvs,nvl,nvm,nqc,nw,offset,tn;
	int NumDistrCenters,DistIndx,DistrStock,nvaccAvail,start,end,vacc_count,leftover;
	GEO_DISTR_CENTER *DC;

	ts=(unsigned short int) (t*P.TimeStepsPerDay);
	nq=nv=nvs=nvl=nvm=nqc=nw=0;
  //Distribution Network Vaccination Loop
  // this will essentially be a Mass Vaccination strategy for now.
	if((P.DoDistributionVaccination)&&(ts>=P.usRoutineImmunisationStartTime))
		{
		NumDistrCenters = VaccineDistributionNetwork->get_num_centers(HEALTH);
    // loop through the distribution centers //This is where we will parallelize, should be enough.
#pragma omp parallel for private(m2vacc,ids,i,j,k,f,n,imc,DistIndx,DistrStock,nvaccAvail,start,end,vacc_count,leftover,DC,offset) reduction(+:nvm,nq,nv,nvs,nvl,nqc,nw) schedule(static,1)
		for(tn=0;tn<P.NumThreads;tn++)
		  for(ids=tn;ids<=NumDistrCenters;ids+=P.NumThreads)
			{
			if(ids<NumDistrCenters)
				{
				DC = VaccineDistributionNetwork->get_center_from_level_index(HEALTH,ids);
				DistrStock = DC->get_stock();
				nvs+=DistrStock;
				k=DC->get_ndemanders();
       //Lets Shuffle that list of microcells, we don't want to be democrats :)
				//FYShuffle_mt(MCellDistrIndex[ids],DC->get_ndemanders(),0);
			//fprintf(stderr,"\nStock available at Store %d = %d %d",ids,DistrStock,DC->get_ndemanders());
				}
			else
				k=noCenterMcellNum;
       //Let's loop over the Mcells for each center
			offset=(int) (ranf_mt(ranf_seed, tn)*((double) k)); //randomise the mcell to go first
			for(imc=0;imc<k;imc++)
				{
				m2vacc= &Mcells[MCellDistrIndex[ids][(imc+offset)%k]];
				n=m2vacc->n;
				nq+=m2vacc->dvacc_count;
				if(ids<NumDistrCenters)
					nvl+=(nvaccAvail = DC->Consume_Avail(m2vacc->dvacc_count));
				else
					nvaccAvail=1000000000;
	 //setup our place in the queue
				start = m2vacc->ndvacc_queue;
				end = start+m2vacc->dvacc_count;
	//Vaccinate
                vacc_count=0;
				int q_count=0;
				int q_missed=0;
				f=1;
				for(i=start;(i<end)&&(f);i++)
					{
					q_count++;
					j=m2vacc->dvacc_queue[i%n];
					if(j<0)
						{
						if(end<start+m2vacc->dvacc_count) end++;
						m2vacc->ndvacc_queue=(m2vacc->ndvacc_queue+1)%n;
						m2vacc->dvacc_count--;
						}
					else if(ts>m2vacc->dvacc_expiry_time[i%n])
						{
						if(end<start+m2vacc->dvacc_count) end++;
						Hosts[j].vacc_queue_pos=-1;
						m2vacc->ndvacc_queue=(m2vacc->ndvacc_queue+1)%n;
						m2vacc->dvacc_count--;
						}
					else if(ts<m2vacc->dvacc_min_vacc_time[i%n])
						{f=0;q_count--;}
					else if(nvaccAvail>0)
						{
						if(!DoDistribVacc(j,ts))
							{
							vacc_count++;
							nvaccAvail--;
							}
						else
							{
							if(end<start+m2vacc->dvacc_count) end++;
							}
						Hosts[j].vacc_queue_pos=-1;
						m2vacc->ndvacc_queue=(m2vacc->ndvacc_queue+1)%n;
						m2vacc->dvacc_count--;
						}
					else
						q_missed++;
					}
				if(ids<NumDistrCenters)
					{
					nv+=vacc_count;
					nvm+=q_missed;
					nqc+=q_count;
					//calculate wastage of doses in phial not used by the end of the day
					if(vacc_count>0)
						{
						double divs=((double) P.UpdatesPerDemogUpdate)/(P.VaccPhialLifetime*((double) (P.TimeStepsPerDay)));
						double vacc_per_div=((double) vacc_count)/divs;
						i=vacc_count;
						j=0;
						while(i>0)
							{
							n=(int) ignpoi_mt(ranf_seed, vacc_per_div, tn);
							if(n>0)
								{
								i-=n;
								if(i>0) j+=n%P.VaccDosesPerPhial;
								}
							}
						j+=(i+n)%P.VaccDosesPerPhial;
						if(j>nvaccAvail) j=nvaccAvail;
						nvaccAvail-=j;
						nw+=j;
						}
					if(nvaccAvail > 0) DC->GiveBack(nvaccAvail);
					}
				}
			}	
		}
	State.cumT+=nw; //use this to store wastage
	j=0;
	for(i=0;i<VaccineDistributionNetwork->num_nodes();i++) j+=VaccineDistributionNetwork->centers()[i]->get_cumulWastage();
	State.cumT+=j-State.dvacc_wastage;
	State.dvacc_wastage=j;
	State.cumTP+=nvm; //use this to store shortfall
	State.cumV+=nv;
	fprintf(stderr,"\n%i (%i) in queue. %i vaccinated. %i doses in stock. %i doses distributed. %i shortfall. %i waste. \n",nqc,nq,nv,nvs,nvl,nvm,nw);
}
#endif

void RecordSample(double t,int n)
{
	int i,j,k,S,L,I,R,D,cumC,cumTC,cumI,cumR,cumD,cumDC,cumDD,cumSDB,cumFC,cumFI; //added cumulative funeral infections: ggilani 24/10/14
	int cumETU; //add number of ETU cases, cumulative ETU cases: ggilani 28/10/14
	int cumH; //added cumulative number of hospitalisations - ggilani 01/07/24
	int cumCT; //added cumulative number of contact traced: ggilani 15/06/17
	int cumCC; //added cumulative number of cases who are contacts: ggilani 28/05/2019
	int cumHQ,cumAC,cumAH,cumAA,cumACS,cumAPC,cumAPA,cumAPCS,numPC,trigDC;
	int cumC_country[MAX_COUNTRIES]; //add cumulative cases per country
	int cumC_adunit[MAX_ADUNITS]; //added cumulative cases per adunit
	int cumDC_adunit[MAX_ADUNITS]; //added cumulative detected cases per adunit
	int cumDD_adunit[MAX_ADUNITS]; //added cumulative detected deaths per adunit <-this and cumulative detected recoveries will help calculate current active detected cases per adunit to use when determine bed capacity increases
	int cumDR_adunit[MAX_ADUNITS]; //added cumulative detected recoveries per adunit - ggilani 13/09/23
	int cumSDB_adunit[MAX_ADUNITS];
	int cumV_adunit[MAX_ADUNITS];
	int cumVG_adunit[MAX_ADUNITS];
	int cumD_adunit[MAX_ADUNITS];
	cell *ct;
	unsigned short int ts;
	int nCase; //added to calculate triggers for individual admin unit funeral controls - ggilani 07/03/2017
	float nsy,*temp2,*temp3,*temp7;

	nsy=DAYS_PER_YEAR/P.SampleStep;
	ts=(unsigned short int) (P.TimeStepsPerDay*t);
	S=L=I=R=D=cumI=cumC=cumDC=cumDD=cumSDB=cumTC=cumFC=cumFI=cumHQ=cumAC=cumAA=cumAH=cumACS=cumAPC=cumAPA=cumAPCS=cumD=cumETU=cumH=cumCT=cumCC=0;
	for(i=0;i<MAX_COUNTRIES;i++) cumC_country[i]=0;
	for(i=0;i<MAX_ADUNITS;i++) cumC_adunit[i]=cumDC_adunit[i]=cumD_adunit[i]=cumDD_adunit[i]=cumDR_adunit[i]=cumSDB_adunit[i]=cumV_adunit[i]=cumVG_adunit[i]=0; //added cumulative detected cases, detected deaths, detected recoveries
#pragma omp parallel for private(i,ct) schedule(static,10000) reduction(+:S,L,I,R,D,cumTC) //added i to private
	for(i=0;i<P.NCP;i++)
		{
		ct=CellLookup[i];
		S+=(int) ct->S;
		L+=(int) ct->L;
		I+=(int) ct->I;
		R+=(int) ct->R;
		D+=(int) ct->D;
		cumTC+=(int) ct->cumTC;
		}
	cumR=R;
	//cumD=D;
	State.sumRad2=0;
	for(j=0;j<P.NumThreads;j++) 
		{
		cumI+=StateT[j].cumI;
		cumC+=StateT[j].cumC;
		cumDC+=StateT[j].cumDC;
		cumFC+=StateT[j].cumFC;
		cumFI+=StateT[j].cumFI; //added cumulative Funeral infections
		cumETU+=StateT[j].cumETU; //added cumulative hospitalisation
		cumH += StateT[j].cumH;
		cumCT+=StateT[j].cumCT; //added contact tracing
		cumCC+=StateT[j].cumCC; //added cases who are contacts
		State.sumRad2+=StateT[j].sumRad2;
		cumHQ+=StateT[j].cumHQ;
		cumAC+=StateT[j].cumAC;
		cumAA+=StateT[j].cumAA;
		cumAPC+=StateT[j].cumAPC;
		cumAPA+=StateT[j].cumAPA;
		cumAPCS+=StateT[j].cumAPCS;
		cumAH+=StateT[j].cumAH;
		cumACS+=StateT[j].cumACS;
		cumD+=StateT[j].cumD;
		cumDD += StateT[j].cumDD;
		cumSDB += StateT[j].cumSDB;
		//add up cumulative country counts: ggilani - 12/11/14
		for(i=0;i<MAX_COUNTRIES;i++) cumC_country[i]+=StateT[j].cumC_country[i];
		//add up cumulative adunit counts: ggilani - 12/11/14
		for(i=0;i<P.NumAdunits;i++) cumC_adunit[i]+=StateT[j].cumC_adunit[i];
		for(i=0;i<P.NumAdunits;i++) cumDC_adunit[i]+=StateT[j].cumDC_adunit[i];
		for (i = 0; i < P.NumAdunits; i++) cumD_adunit[i] += StateT[j].cumD_adunit[i];
		for (i = 0; i < P.NumAdunits; i++) cumSDB_adunit[i] += StateT[j].cumSDB_adunit[i];
		for (i = 0; i < P.NumAdunits; i++) cumV_adunit[i] += StateT[j].cumV_adunit[i];
		for (i = 0; i < P.NumAdunits; i++) cumVG_adunit[i] += StateT[j].cumVG_adunit[i]; //added geographically targeted vaccination by admin_unit
		for (i = 0; i < P.NumAdunits; i++) cumDD_adunit[i] += StateT[j].cumDD_adunit[i]; //detected deaths
		for (i = 0; i < P.NumAdunits; i++) cumDR_adunit[i] += StateT[j].cumDR_adunit[i]; //detected recoveries
		if(State.maxRad2<StateT[j].maxRad2) State.maxRad2=StateT[j].maxRad2;
		}
	for(j=0;j<P.NumThreads;j++)
		StateT[j].maxRad2=State.maxRad2;
	TimeSeries[n].t=t;
	TimeSeries[n].S=(double) S;
	TimeSeries[n].L=(double) L;
	TimeSeries[n].I=(double) I;
	TimeSeries[n].R=(double) R;
	TimeSeries[n].D=(double) D;
	TimeSeries[n].incI=(double) (cumI-State.cumI);
	TimeSeries[n].incC=(double) (cumC-State.cumC);
	TimeSeries[n].incFC=(double) (cumFC-State.cumFC);
	TimeSeries[n].incFI=(double) (cumFI-State.cumFI); //added funeral infections
	TimeSeries[n].incETU=(double) (cumETU-State.cumETU); //added incidence of hospitalisation
	TimeSeries[n].incH = (double)(cumH - State.cumH);
	TimeSeries[n].incCT=(double) (cumCT-State.cumCT); // added contact tracing
	TimeSeries[n].incCC = (double)(cumCC-State.cumCC); // added cases who are contacts
	TimeSeries[n].incDC=(double) (cumDC-State.cumDC); //added incidence of detected cases
	TimeSeries[n].incD = (double)(cumD - State.cumD); //added incidence of detected deaths
	TimeSeries[n].incDD = (double)(cumDD - State.cumDD); //added incidence of detected deaths
	TimeSeries[n].incSDB = (double)(cumSDB - State.cumSDB); //added incidence of safe burials
	TimeSeries[n].incTC=(double) (cumTC-State.cumTC);
	TimeSeries[n].incR=(double) (cumR-State.cumR);
	TimeSeries[n].incD=(double) (cumD-State.cumD);
	TimeSeries[n].incHQ=(double) (cumHQ-State.cumHQ);
	TimeSeries[n].incAC=(double) (cumAC-State.cumAC);
	TimeSeries[n].incAH=(double) (cumAH-State.cumAH);
	TimeSeries[n].incAA=(double) (cumAA-State.cumAA);
	TimeSeries[n].incACS=(double) (cumACS-State.cumACS);
	TimeSeries[n].incAPC=(double) (cumAPC-State.cumAPC);
	TimeSeries[n].incAPA=(double) (cumAPA-State.cumAPA);
	TimeSeries[n].incAPCS=(double) (cumAPCS-State.cumAPCS);
	TimeSeries[n].cumT=State.cumT;
	TimeSeries[n].cumUT=State.cumUT;
	TimeSeries[n].cumTP=State.cumTP;
	TimeSeries[n].cumV=State.cumV;
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
	
	TimeSeries[n].cumDC=cumDC;
	TimeSeries[n].cumDD = cumDD;
	TimeSeries[n].cumSDB = cumSDB;
	//incidence per country
	for(i=0;i<MAX_COUNTRIES;i++) TimeSeries[n].incC_country[i]=(double)(cumC_country[i]-State.cumC_country[i]);
	trigDC=cumDC;
	if(n>=P.TriggersSamplingInterval) trigDC-=TimeSeries[n-P.TriggersSamplingInterval].cumDC;
	if(trigDC>State.trigDC) State.trigDC=trigDC;
	State.S=S;
	State.L=L;
	State.I=I;
	State.R=R;
	State.D=D;
	State.cumI=cumI;
	State.cumDC=cumDC;
	State.cumDD = cumDD;
	State.cumSDB = cumSDB;
	State.cumTC=cumTC;
	State.cumFC=cumFC;
	State.cumFI=cumFI; //added cumulative funeral infections
	State.cumETU=cumETU; //added cumulative hospitalisation
	State.cumH = cumH;
	State.cumCT=cumCT; //added cumulative contact tracing
	State.cumCC=cumCC; //added cumulative cases who are contacts
	State.cumC=cumC;
	State.cumR=cumR;
	State.cumD=cumD;
	State.cumHQ=cumHQ;
	State.cumAC=cumAC;
	State.cumAH=cumAH;
	State.cumAA=cumAA;
	State.cumACS=cumACS;
	State.cumAPC=cumAPC;
	State.cumAPA=cumAPA;
	State.cumAPCS=cumAPCS;
	//update cumulative cases per country
	for(i=0;i<MAX_COUNTRIES;i++) State.cumC_country[i]=cumC_country[i];
	//update overall state variable for cumulative cases per adunit
	for(i=0;i<P.NumAdunits;i++) State.cumC_adunit[i]+=cumC_adunit[i]; //here we increment rather than set equal to as the cumulative admin unit counts per thread are reset to zero after recording a sample
	for(i=0;i<P.NumAdunits;i++) State.cumDC_adunit[i]+=cumDC_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumD_adunit[i] += cumD_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumDD_adunit[i] += cumDD_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumDR_adunit[i] += cumDR_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumDD_adunit[i] += cumDD_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumSDB_adunit[i] += cumSDB_adunit[i];
	for (i = 0; i < P.NumAdunits; i++) State.cumV_adunit[i] += cumV_adunit[i]; //added vaccination: ggilani 14/01/25
	for (i = 0; i < P.NumAdunits; i++) State.cumVG_adunit[i] += cumVG_adunit[i]; //added geographically targeted vaccination: ggilani 22/08/25
	TimeSeries[n].rmsRad=(State.cumI>0) ? sqrt(State.sumRad2/((double)State.cumI)):0;
	TimeSeries[n].maxRad=sqrt(State.maxRad2);
	TimeSeries[n].extinct=((((P.SmallEpidemicCases>=0)&&(State.R<=P.SmallEpidemicCases))||(P.SmallEpidemicCases<0))&&(State.I+State.L==0))?1:0;
	TimeSeries[n].detected = P.OutbreakDetected;
	for(i=0;i<NUM_AGE_GROUPS;i++)
		{
		TimeSeries[n].incCa[i]=TimeSeries[n].incIa[i]=TimeSeries[n].incDa[i]= TimeSeries[n].incDCa[i] = TimeSeries[n].incETUa[i] = TimeSeries[n].incHa[i] = TimeSeries[n].incVa[i] = 0;
		for(j=0;j<P.NumThreads;j++)
			{
			TimeSeries[n].incCa[i]+=(double) StateT[j].cumCa[i];
			TimeSeries[n].incIa[i]+=(double) StateT[j].cumIa[i];
			TimeSeries[n].incDa[i]+=(double) StateT[j].cumDa[i];
			TimeSeries[n].incDCa[i] += (double)StateT[j].cumDCa[i];
			TimeSeries[n].incETUa[i] += (double)StateT[j].cumETUa[i];
			TimeSeries[n].incHa[i] += (double)StateT[j].cumHa[i];
			//TimeSeries[n].incVa[i] += (double)StateT[j].cumVa[i];
			}
		TimeSeries[n].incVa[i] += State.cumVa[i]; //changed vacc by age to here, because it is incremented outside of the threads: ggilani 22/02/22
		}
	for(i=0;i<P.EvolResistNumTypes;i++)
		{
		TimeSeries[n].incC_resist[i]=TimeSeries[n].incI_resist[i]=TimeSeries[n].cumT_resist[i]=0;
		for(j=0;j<P.NumThreads;j++)
			{
			TimeSeries[n].incC_resist[i]+=(double) StateT[j].cumC_resist[i];
			TimeSeries[n].incI_resist[i]+=(double) StateT[j].cumI_resist[i];
			TimeSeries[n].cumT_resist[i]+=(double) StateT[j].cumT_resist[i];
			StateT[j].cumC_resist[i]=StateT[j].cumI_resist[i]=0;
			}
		}
	for(i=0;i<2;i++)
		{
		TimeSeries[n].incC_keyworker[i]=TimeSeries[n].incI_keyworker[i]=TimeSeries[n].cumT_keyworker[i]=TimeSeries[n].incD_keyworker[i]=0;
		for(j=0;j<P.NumThreads;j++)
			{
			TimeSeries[n].incC_keyworker[i]+=(double) StateT[j].cumC_keyworker[i];
			TimeSeries[n].incI_keyworker[i]+=(double) StateT[j].cumI_keyworker[i];
			TimeSeries[n].cumT_keyworker[i]+=(double) StateT[j].cumT_keyworker[i];
			TimeSeries[n].incD_keyworker[i] += (double)StateT[j].cumD_keyworker[i];
			StateT[j].cumC_keyworker[i]=StateT[j].cumI_keyworker[i]=StateT[j].cumD_keyworker[i]=0;
			}
		}

	for(i=0;i<INFECT_TYPE_MASK;i++)
		{
		TimeSeries[n].incItype[i]=0;
		for(j=0;j<P.NumThreads;j++)
			{
			TimeSeries[n].incItype[i]+=(double) StateT[j].cumItype[i];
			StateT[j].cumItype[i]=0;
			}
		}
	if(P.DoAdUnits)
		for(i=0;i<=P.NumAdunits;i++)
		{
			TimeSeries[n].incI_adunit[i]=TimeSeries[n].incC_adunit[i]= TimeSeries[n].incDC_adunit[i]= TimeSeries[n].incD_adunit[i]= TimeSeries[n].incDD_adunit[i]= TimeSeries[n].incDR_adunit[i]= TimeSeries[n].incSDB_adunit[i]= TimeSeries[n].incV_adunit[i]= TimeSeries[n].incVG_adunit[i]=TimeSeries[n].incETU_adunit[i]=TimeSeries[n].incH_adunit[i]=TimeSeries[n].cumT_adunit[i]=TimeSeries[n].incCT_adunit[i]=TimeSeries[n].incCC_adunit[i]=TimeSeries[n].capETU_adunit[i]= TimeSeries[n].ETU_adunit[i]=0; //added detected cases: ggilani 03/02/15
			for(j=0;j<P.NumThreads;j++)
			{
				TimeSeries[n].incI_adunit[i]+=(double) StateT[j].cumI_adunit[i];
				TimeSeries[n].incC_adunit[i]+=(double) StateT[j].cumC_adunit[i];
				TimeSeries[n].incDC_adunit[i]+=(double) StateT[j].cumDC_adunit[i]; //added detected cases: ggilani 03/02/15
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
				TimeSeries[n].incETU_adunit[i]+=(double) StateT[j].cumETU_adunit[i]; //added hospitalisation
				TimeSeries[n].incH_adunit[i] += (double)StateT[j].cumH_adunit[i];
				TimeSeries[n].incCT_adunit[i]+=(double) StateT[j].cumCT_adunit[i]; //added contact tracing: ggilani 15/06/17
				TimeSeries[n].incCC_adunit[i] += (double)StateT[j].cumCC_adunit[i]; //added cases who are contacts: ggilani 28/05/2019
				TimeSeries[n].cumT_adunit[i]+=(double) StateT[j].cumT_adunit[i];
				StateT[j].cumI_adunit[i]=StateT[j].cumC_adunit[i]=StateT[j].cumH_adunit[i]=StateT[j].cumDC_adunit[i]= StateT[j].cumD_adunit[i]=StateT[j].cumDD_adunit[i]= StateT[j].cumDR_adunit[i]=StateT[j].cumSDB_adunit[i]=StateT[j].cumCT_adunit[i]=StateT[j].cumV_adunit[i]= StateT[j].cumVG_adunit[i]=StateT[j].cumETU_adunit[i]= StateT[j].cumCC_adunit[i]=0; //added hospitalisation, detected cases, contact tracing: ggilani 03/02/15, 15/06/17
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
		
		if(P.DoContactTracing) //added this to print out total number of beds in use at each time point: ggilani 31/10/14
		{
			for(i=0;i<P.NumAdunits;i++)
			{
				TimeSeries[n].CT_adunit[i]=(double) AdUnits[i].nct;
			}
		}
	if((P.DoPlaces)&&(t>=P.PlaceCloseTimeStart))
		{
		for(i=0;i<NUM_PLACE_TYPES;i++)
			{
			numPC=0;
			for(j=0;j<P.Nplace[i];j++)
				if(PLACE_CLOSED(i,j)) numPC++;
			State.NumPlacesClosed[i]=numPC;
			TimeSeries[n].PropPlacesClosed[i]=((double) numPC)/((double) P.Nplace[i]);
			}
		}
	//if (State.cumC > P.NumUndetectedInfPreOutbreakAlert)
	//{
	//	P.ControlPropCasesId = P.PostAlertControlPropCasesId;
	//}
	if(State.cumDC>=P.PreControlClusterIdCaseThreshold)
	{
		if (P.OutbreakDetected == 0)
		{
			P.OutbreakDetected = 1; //mart outbreak as detected
			P.PropHospSeek *= P.RelChangeHospSeekPostOutbreak;
		}
		if(P.DoGlobalTriggers)
			{
			if(P.DoPerCapitaTriggers)
				D=(int) floor(((double) State.trigDC)*P.GlobalIncThreshPop/((double) P.N));
			else
				D=State.trigDC;
			if(D>=P.TreatCellIncThresh) 
				{
				if(P.TreatTimeStart>=1e10) P.TreatTimeStart=t+P.TreatTimeStartBase;
				if(P.CaseIsolationTimeStart>=1e10) P.CaseIsolationTimeStart=t+P.CaseIsolationTimeStartBase;
				if(P.HQuarantineTimeStart>=1e10) P.HQuarantineTimeStart=t+P.HQuarantineTimeStartBase;
				}
			if(D>=P.VaccCellIncThresh)
				{
				if (P.VaccTimeStart >= 1e10)
					{
					P.VaccTimeStart = t + P.VaccTimeStartBase;
					}
				if (P.VaccNewCoursesStartTime >= 1e10)
				{
					P.VaccNewCoursesStartTime = t + P.VaccNewCoursesStartTimeBase;
				}
				}
			if(D>=P.SocDistCellIncThresh)
				{
				if(P.SocDistTimeStart>=1e10) P.SocDistTimeStart=t+P.SocDistTimeStartBase;
				}
			if(D>=P.PlaceCloseCellIncThresh)
				{
				if(P.PlaceCloseTimeStart>=1e10) P.PlaceCloseTimeStart=t+P.PlaceCloseTimeStartBase;
				if((P.PlaceCloseTimeStart2>=1e10)&&(t>=P.PlaceCloseDuration+P.PlaceCloseTimeStart))
					{
					P.PlaceCloseTimeStart=t+P.PlaceCloseTimeStartBase2-P.PlaceCloseTimeStartBase;
					P.PlaceCloseDuration=P.PlaceCloseDuration2;
					}
				}
			if(D>=P.MoveRestrCellIncThresh)
				{
				if(P.MoveRestrTimeStart>=1e10) P.MoveRestrTimeStart=t+P.MoveRestrTimeStartBase;
				}
			if(D>=P.KeyWorkerProphCellIncThresh)
				{
				if(P.KeyWorkerProphTimeStart>=1e10) P.KeyWorkerProphTimeStart=t+P.KeyWorkerProphTimeStartBase;
				}
			/*if (D >= P.GeoVaccCellIncThresh)
			{
				if (P.GeoVaccTimeStart >= 1e10) P.GeoVaccTimeStart = t + P.GeoVaccTimeStartBase;
			}*/
			//added this to enable funeral control start time when not doing funeral controls on an admin unit basis
			if((D>=P.FuneralControlCellIncThresh))
				{
					if(P.FuneralControlTimeStart>=1e10)
					{
						P.FuneralControlTimeStart=t+P.FuneralControlTimeStartBase;
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
			if(P.TreatTimeStart>=1e10) P.TreatTimeStart=t+P.TreatTimeStartBase;
			if(P.CaseIsolationTimeStart>=1e10) P.CaseIsolationTimeStart=t+P.CaseIsolationTimeStartBase;
			if(P.HQuarantineTimeStart>=1e10) P.HQuarantineTimeStart=t+P.HQuarantineTimeStartBase;
			if (P.VaccTimeStart >= 1e10)
			{
				P.VaccTimeStart = t + P.VaccTimeStartBase;
				//fprintf(stderr, "t=%lg, P.VaccTimeStart=%lg, P.VaccTimeStartBase=%lg\n", t, P.VaccTimeStart, P.VaccTimeStartBase);
			}
			if (P.VaccNewCoursesStartTime >= 1e10)
			{
				P.VaccNewCoursesStartTime = t + P.VaccNewCoursesStartTimeBase;
			}
			if(P.SocDistTimeStart>=1e10) P.SocDistTimeStart=t+P.SocDistTimeStartBase;
			if(P.PlaceCloseTimeStart>=1e10) P.PlaceCloseTimeStart=t+P.PlaceCloseTimeStartBase;
			if(P.MoveRestrTimeStart>=1e10) P.MoveRestrTimeStart=t+P.MoveRestrTimeStartBase;
			if(P.KeyWorkerProphTimeStart>=1e10) P.KeyWorkerProphTimeStart=t+P.KeyWorkerProphTimeStartBase;
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
	
		if(P.AirportCloseTimeStart>=1e10) P.AirportCloseTimeStart=t+P.AirportCloseTimeStartBase;
	}
	if((P.PlaceCloseIndepThresh>0)&&(((double) State.cumDC)>=P.PlaceCloseIndepThresh))
		{
		if(P.PlaceCloseTimeStart>=1e10) P.PlaceCloseTimeStart=t+P.PlaceCloseTimeStartBase;
		}
	if(P.OutputBitmap>=1)
		{
		TSMean=TSMeanNE;TSVar=TSVarNE;
		CaptureBitmap(n,0);
		OutputBitmap(t,0);
		}

}

void RecordInfTypes(void)
{
	int i,j,k,l,lc,lc2,b,c,n,nf,i2;
	double *res,*res_av,*res_var,t,s;

	for(n=0;n<P.NumSamples;n++)
		{
		for(i=0;i<INFECT_TYPE_MASK;i++) TimeSeries[n].Rtype[i]=0;
		for(i=0;i<NUM_AGE_GROUPS;i++) TimeSeries[n].Rage[i]=0;
		TimeSeries[n].Rdenom=0;
		}
	for(i=0;i<INFECT_TYPE_MASK;i++) inftype[i]=0;
	for(i=0;i<MAX_COUNTRIES;i++) infcountry[i]=0;
	for(i=0;i<MAX_SEC_REC;i++) 
		for(j=0;j<MAX_GEN_REC;j++)
			indivR0[i][j]=0;
	for(i=0;i<=MAX_HOUSEHOLD_SIZE;i++)
		for(j=0;j<=MAX_HOUSEHOLD_SIZE;j++)
			inf_household[i][j]=case_household[i][j]=0;
	for(b=0;b<P.NC;b++)
		if((Cells[b].S!=Cells[b].n)||(Cells[b].R>0))
			for(c=0;c<Cells[b].n;c++)
				Hosts[Cells[b].members[c]].listpos=0;
//	for(b=0;b<P.NC;b++)
//		if((Cells[b].S!=Cells[b].n)||(Cells[b].R>0))
			{
			j=k=l=lc=lc2=0;t=1e10;
//			for(c=0;c<Cells[b].n;c++)
			for(i=0;i<P.N;i++)
				{
//				i=Cells[b].members[c];
				if(j==0) j=k=Households[Hosts[i].hh].nh;
				if((Hosts[i].inf!=0)&&(Hosts[i].inf!=4))
					{
					if(Hosts[i].latent_time*P.TimeStep<=P.SampleTime)
						TimeSeries[(int) (Hosts[i].latent_time*P.TimeStep/P.SampleStep)].Rdenom++;
					infcountry[Mcells[Hosts[i].mcell].country]++;
					if(abs(Hosts[i].inf)<3)
						l=-1;
					else if(l>=0)
						l++;
					if((l>=0)&&((Hosts[i].inf==-3)||(Hosts[i].inf==-5)))
						{
						lc2++;
						if(Hosts[i].latent_time*P.TimeStep<=t) // This convoluted logic is to pick up households where the index is symptomatic
							{lc=1;t=Hosts[i].latent_time*P.TimeStep;}
						}
					else if((l>0)&&(Hosts[i].latent_time*P.TimeStep<t))
						{lc=0;t=Hosts[i].latent_time*P.TimeStep;}
					i2=Hosts[i].infector;
					if(i2>=0)
						{
						Hosts[i2].listpos++;
						if(Hosts[i2].latent_time*P.TimeStep<=P.SampleTime)
							{
							TimeSeries[(int) (Hosts[i2].latent_time*P.TimeStep/P.SampleStep)].Rtype[Hosts[i].infect_type%INFECT_TYPE_MASK]++;
							TimeSeries[(int) (Hosts[i2].latent_time*P.TimeStep/P.SampleStep)].Rage[HOST_AGE_GROUP(i)]++;
							}
						}
					}
				inftype[Hosts[i].infect_type%INFECT_TYPE_MASK]++;
				j--;
				if(j==0) 
					{
					if(l<0) l=0;
					inf_household[k][l]++;
					case_household[k][lc2]++; //now recording total symptomatic cases, rather than infections conditional on symptomatic index
					l=lc=lc2=0;t=1e10;
					}
				}
			}
		for(b=0;b<P.NC;b++)
			if((Cells[b].S!=Cells[b].n)||(Cells[b].R>0))
				for(c=0;c<Cells[b].n;c++)
					{
					i=Cells[b].members[c];
					if((abs(Hosts[i].inf)==3)||(abs(Hosts[i].inf)==5))
						{
						l=Hosts[i].infect_type/INFECT_TYPE_MASK;
						if((l<MAX_GEN_REC)&&(Hosts[i].listpos<MAX_SEC_REC)) indivR0[Hosts[i].listpos][l]++;
						}
					}
/* 	if(!TimeSeries[P.NumSamples-1].extinct) */
		{
		for(i=0;i<INFECT_TYPE_MASK;i++) inftype_av[i]+=inftype[i];
		for(i=0;i<MAX_COUNTRIES;i++)
			{
			infcountry_av[i]+=infcountry[i];
			if(infcountry[i]>0) infcountry_num[i]++;
			}
		for(i=0;i<MAX_SEC_REC;i++)
			for(j=0;j<MAX_GEN_REC;j++)
				indivR0_av[i][j]+=indivR0[i][j];
		for(i=0;i<=MAX_HOUSEHOLD_SIZE;i++)
			for(j=0;j<=MAX_HOUSEHOLD_SIZE;j++)
				{
				inf_household_av[i][j]+=inf_household[i][j];
				case_household_av[i][j]+=case_household[i][j];
				}
		}
	for(n=0;n<P.NumSamples;n++)
		{
		s=0;
		if(TimeSeries[n].Rdenom==0) TimeSeries[n].Rdenom=1e-10;
		for(i=0;i<NUM_AGE_GROUPS;i++)
			TimeSeries[n].Rage[i]/=TimeSeries[n].Rdenom;
		for(i=0;i<INFECT_TYPE_MASK;i++)
			s+=(TimeSeries[n].Rtype[i]/=TimeSeries[n].Rdenom);
		TimeSeries[n].Rdenom=s;
		}
	nf=(sizeof(results)-3*sizeof(float *))/sizeof(double);
	if(!P.DoAdUnits) nf-=MAX_ADUNITS;
	fprintf(stderr,"extinct=%i (%i)\n",TimeSeries[P.NumSamples-1].extinct,P.NumSamples-1);
	if(TimeSeries[P.NumSamples-1].extinct)
		{TSMean=TSMeanE;TSVar=TSVarE;P.NRactE++;}
	else
		{TSMean=TSMeanNE;TSVar=TSVarNE;P.NRactNE++;}
	s=0;
	for(n=0;n<P.NumSamples;n++)
		{
		if(s<TimeSeries[n].incC) {s=TimeSeries[n].incC;t=P.SampleStep*((double) n);}
		res=(double *) &TimeSeries[n];
		res_av=(double *) &TSMean[n];
		res_var=(double *) &TSVar[n];
		for(i=0;i<nf;i++)
			{
			res_av[i]+=res[i];
			res_var[i]+=res[i]*res[i];
			}
		if(TSMean[n].cumTmax<TimeSeries[n].cumT) TSMean[n].cumTmax=TimeSeries[n].cumT;
		if(TSMean[n].cumVmax<TimeSeries[n].cumV) TSMean[n].cumVmax=TimeSeries[n].cumV;
		}
	PeakHeightSum+=s;
	PeakHeightSS+=s*s;
	PeakTimeSum+=t;
	PeakTimeSS+=t*t;
}

void DoReborn(int ai)
// This transfers a person straight from immune to susceptible. Used to update demography for SIR model.
{	
	person *a;
	int j,c;
	double x,y;

	a=Hosts+ai;
	c=a->pcell;
	if(abs(a->inf)>=3)
		{
		a->inf=0;
		a->infector=-1;
		if(a->listpos>Cells[c].S+Cells[c].L+Cells[c].I)
			{
			Cells[c].susceptible[a->listpos]=Cells[c].susceptible[Cells[c].S+Cells[c].L+Cells[c].I];
			Hosts[Cells[c].susceptible[a->listpos]].listpos=a->listpos;
			}
		if(Cells[c].I>0)
			{
			Cells[c].susceptible[Cells[c].S+Cells[c].L+Cells[c].I]=Cells[c].susceptible[Cells[c].S+Cells[c].L];
			Hosts[Cells[c].susceptible[Cells[c].S+Cells[c].L+Cells[c].I]].listpos=Cells[c].S+Cells[c].L+Cells[c].I;
			}
		if(Cells[c].L>0)
			{
			Cells[c].susceptible[Cells[c].S+Cells[c].L]=Cells[c].susceptible[Cells[c].S];
			Hosts[Cells[c].susceptible[Cells[c].S+Cells[c].L]].listpos=Cells[c].S+Cells[c].L;
			}
		Cells[c].susceptible[Cells[c].S]=ai;
		a->listpos=Cells[c].S;
		Cells[c].latent++;
		Cells[c].infected++;
		Cells[c].S++;
		Cells[c].R--;
		if(P.OutputBitmap)
			{
			x=((int) (Households[a->hh].loc_x*P.scalex))-P.bminx;
			y=((int) (Households[a->hh].loc_y*P.scaley))-P.bminy;
			if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
				{
				j=y*bmh->width+x;
				if((j<bmh->imagesize)&&(j>=0))
					{
#pragma omp atomic
					bmi3[j]--;
					}
				}
			}
		}
}

void DoImmune(int ai)
// This transfers a person straight from susceptible to immune. Used to start a run with a partially immune population.
{	
	person *a;
	int j,c;
	double x,y;

	a=Hosts+ai;
	if(a->inf==0)
		{
		c=a->pcell;
		a->inf=4;
		Cells[c].S--;
		if(a->listpos<Cells[c].S)
			{
			Cells[c].susceptible[a->listpos]=Cells[c].susceptible[Cells[c].S];
			Hosts[Cells[c].susceptible[a->listpos]].listpos=a->listpos;
			}
		if(Cells[c].L>0)
			{
			Cells[c].susceptible[Cells[c].S]=Cells[c].susceptible[Cells[c].S+Cells[c].L];
			Hosts[Cells[c].susceptible[Cells[c].S]].listpos=Cells[c].S;
			}
		if(Cells[c].I>0)
			{
			Cells[c].susceptible[Cells[c].S+Cells[c].L]=Cells[c].susceptible[Cells[c].S+Cells[c].L+Cells[c].I];
			Hosts[Cells[c].susceptible[Cells[c].S+Cells[c].L]].listpos=Cells[c].S+Cells[c].L;
			}
		if(a->listpos<Cells[c].S+Cells[c].L+Cells[c].I)
			{
			Cells[c].susceptible[Cells[c].S+Cells[c].L+Cells[c].I]=ai;
			a->listpos=Cells[c].S+Cells[c].L+Cells[c].I;
			}
		Cells[c].latent--;
		Cells[c].infected--;
		Cells[c].R++;
		if(P.OutputBitmap)
			{
			x=((int) (Households[a->hh].loc_x*P.scalex))-P.bminx;
			y=((int) (Households[a->hh].loc_y*P.scaley))-P.bminy;
			if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
				{
				j=y*bmh->width+x;
				if((j<bmh->imagesize)&&(j>=0))
					{
#pragma omp atomic
					bmi3[j]++;
					}
				}
			}
		}
}


void DoInfect(int ai,double t,int tn, int run) //added int as argument to DoInfect to record run number: ggilani - 15/10/14
{

	int i,j;
	unsigned short int ts;
	double q,x,y;
	person *a;

	a=Hosts+ai;
	
	if(a->inf==0)
		{
		ts=(unsigned short int) (P.TimeStepsPerDay*t);
		a->inf=1;
		a->infection_time=(unsigned short int) ts;
		if(a->infector>=0) 
			{
			a->resist=Hosts[a->infector].resist;
			a->base_inf_level=Hosts[a->infector].base_inf_level;
			if(abs(Hosts[a->infector].inf)==6) //if host has an infector (i.e. they are not a seed case), check to see if infection status flag of host is 6, i.e. that it is a funeral contact, ggilani: 24/10/14
				{
					StateT[tn].cumFI++;
				}
			}
		StateT[tn].cumI++;
		StateT[tn].cumItype[a->infect_type%INFECT_TYPE_MASK]++;
		StateT[tn].cumIa[HOST_AGE_GROUP(ai)]++;
		x=(Households[a->hh].loc_x-P.LocationInitialInfection[0][0]);
		y=(Households[a->hh].loc_y-P.LocationInitialInfection[0][1]);
		q=x*x+y*y;
		StateT[tn].sumRad2+=q;
		if(q>StateT[tn].maxRad2) StateT[tn].maxRad2=q;
			{
			Cells[a->pcell].S--;
			if(Cells[a->pcell].S>0)
				{
				Cells[a->pcell].susceptible[a->listpos]=Cells[a->pcell].susceptible[Cells[a->pcell].S];
				Hosts[Cells[a->pcell].susceptible[a->listpos]].listpos=a->listpos;
				}
			Cells[a->pcell].latent--;
			Cells[a->pcell].L++;
			a->listpos=Cells[a->pcell].S;
			Cells[a->pcell].latent[0]=ai;
			}
		if((HOST_TREATED(ai))&&(a->resist<(MAX_NUM_RESIST_TYPES-1))&&(ranf_mt(ranf_seed, tn)<P.EvolResistProphMutationRate)) a->resist++;
		StateT[tn].cumI_resist[a->resist]++;
		StateT[tn].cumI_keyworker[a->keyworker]++;
		if(P.DoLatent)
			{
			i=(int) floor((q=ranf_mt(ranf_seed, tn)*CDF_RES));
			q-=((double) i);
			a->latent_time=(unsigned short int) floor(0.5+(t-P.LatentPeriod*log(q*P.latent_icdf[i+1]+(1.0-q)*P.latent_icdf[i]))*P.TimeStepsPerDay);
			}
		else
			a->latent_time=(unsigned short int) (t*P.TimeStepsPerDay);
		if(P.DoAdUnits)
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
		if(P.OutputBitmap)
			{
				if((P.OutputBitmapDetected==0)||((P.OutputBitmapDetected==1)&&(Hosts[ai].detected==1)))
				{
					x=((int) (Households[a->hh].loc_x*P.scalex))-P.bminx;
					y=((int) (Households[a->hh].loc_y*P.scaley))-P.bminy;
					if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
					{
						j=y*bmh->width+x;
						if((j<bmh->imagesize)&&(j>=0))
						{
#pragma omp atomic
							bmi2[j]++;
						}
					}
				}
			}
		//added this to record event if flag is set to 1 : ggilani - 10/10/2014
		if(P.DoRecordInfEvents)
		{
			if(*nEvents<P.MaxInfEvents)
			{
				RecordEvent(t,ai,run,0,tn); //added int as argument to RecordEvent to record run number: ggilani - 15/10/14
			}
		}
		if((t>0)&&(P.DoOneGen))
			{
			DoIncub(ai,ts,tn,run);
			DoCase(ai,t,ts,tn);			
			DoRecover(ai,run,tn);
			}
		}
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

	bi=Hosts[ai].infector;

	//Save information to event
#pragma omp critical (inf_event)
	{
		InfEventLog[*nEvents].run=run;
		InfEventLog[*nEvents].type=type;
		InfEventLog[*nEvents].t=t;
		InfEventLog[*nEvents].infectee_ind=ai;
		InfEventLog[*nEvents].infectee_adunit=Mcells[Hosts[ai].mcell].adunit;
		InfEventLog[*nEvents].infectee_x=Households[Hosts[ai].hh].loc_x+P.SpatialBoundingBox[0];
		InfEventLog[*nEvents].infectee_y=Households[Hosts[ai].hh].loc_y+P.SpatialBoundingBox[1];
		InfEventLog[*nEvents].listpos = Hosts[ai].listpos;
		InfEventLog[*nEvents].infectee_cell = Hosts[ai].pcell;
		InfEventLog[*nEvents].thread=tn;
		if(type==0) //infection event - record time of onset of infector and infector
		{
			InfEventLog[*nEvents].infector_ind=bi;
			if(bi<0)
			{
			InfEventLog[*nEvents].t_infector=-1;
			InfEventLog[*nEvents].infector_cell=-1;
			}
			else
			{
			InfEventLog[*nEvents].t_infector=(int)(Hosts[bi].infection_time/P.TimeStepsPerDay);
			InfEventLog[*nEvents].infector_cell = Hosts[bi].pcell;
			}
		}
		else if(type==1) //onset event - record infectee's onset time
		{
			InfEventLog[*nEvents].t_infector=(int)(Hosts[ai].infection_time/P.TimeStepsPerDay);
		}
		else if((type==2)||(type==3)) //recovery or death event - record infectee's onset time
		{
			InfEventLog[*nEvents].t_infector=(int)(Hosts[ai].latent_time/P.TimeStepsPerDay);
		}

		//increment the index of the infection event
		(*nEvents)++;
	}
	
}

void DoIncub(int ai,unsigned short int ts,int tn, int run)
	{	
	int i;
	person *a;
	double q,ti,cfr;
	int age,day; 

	age=HOST_AGE_GROUP(ai);
	if(age>=NUM_AGE_GROUPS) age=NUM_AGE_GROUPS-1;

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

	a=Hosts+ai;
	if(a->inf==1)
		{
		if(P.InfectiousnessSD==0)
			a->infectiousness=P.AgeInfectiousness[age];
		else
			//a->infectiousness=P.AgeInfectiousness[age]*gen_beta_mt(ranf_seed, P.InfectiousnessBetaA,P.InfectiousnessBetaB,tn);
			a->infectiousness=P.AgeInfectiousness[age]*gen_gamma_mt(ranf_seed, P.InfectiousnessGamA,P.InfectiousnessGamR,tn);
		q=P.ProportionSymptomatic[age]
			*(HOST_TREATED(ai)?P.EvolResistRelTreatSympDrop[Hosts[ai].resist]:1)
			*(HOST_VACCED(ai)?(1-P.VaccSympDrop):1);
		if(ranf_mt(ranf_seed, tn)<q) 
			{
			a->inf=-1;
			a->infectiousness=-P.SymptInfectiousness*a->infectiousness;
			}
		else
		{
			a->inf = 2;
		}
		// commenting this out for now to redo cfr/mortality taking into account the impact of vaccination
		
		if(P.DoInfectiousnessProfile)
			a->recovery_time=a->latent_time+(unsigned short int) (P.InfectiousPeriod*P.TimeStepsPerDay);
		else
			{
			i=(int) floor(q=ranf_mt(ranf_seed, tn)*CDF_RES);
			q-=((double) i);
			ti=-P.InfectiousPeriod*log(q*P.infectious_icdf[i+1]+(1.0-q)*P.infectious_icdf[i]);
			a->recovery_time=a->latent_time+(unsigned short int) floor(0.5+(ti*P.TimeStepsPerDay));
			}

		if (P.DoMortality) //added DoMortality to allow different approaches to assigning mortality: ggilani - 22/10/2014
		{
			if (P.DoEventMortality)
			{
				day = (int)ceil(P.TimeStep * (a->recovery_time - a->latent_time));
				if (ranf_mt(ranf_seed, tn) <= P.RecoveryProb[day])
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
				
				if (ranf_mt(ranf_seed, tn) < cfr)
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
		if(Cells[a->pcell].L>0)
			{
			Cells[a->pcell].susceptible[a->listpos]=Cells[a->pcell].latent[Cells[a->pcell].L];
			Hosts[Cells[a->pcell].susceptible[a->listpos]].listpos=a->listpos;
			}
		Cells[a->pcell].infected--;
		a->listpos=Cells[a->pcell].S+Cells[a->pcell].L;
		Cells[a->pcell].I++;
		Cells[a->pcell].infected[0]=ai;
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

void DoDetectedCase(int ai,double t,unsigned short int ts,int tn)
{	

	int j,k,f,j1,j2,age;
	int i, ii, i1, i2, i3, k2, k3, l, vaccflag;
	int currentRing, cnt, nActiveCases;
	person *a;

	a=Hosts+ai;
	age = HOST_AGE_GROUP(ai);
	if((!P.RestrictTreatToTarget)||(Mcells[a->mcell].country==P.TargetCountry))
	{
		//only increase triggers if outbreak has been declared
		//if (P.OutbreakDetected)
		//{
			if (Mcells[a->mcell].treat_trig < USHRT_MAX - 1) Mcells[a->mcell].treat_trig++;
			if ((!P.OnlyDoGeoVaccWhenNoRing)|| (Hosts[ai].vacc_accept > P.ProbEstablishRing))
			{
				if (Mcells[a->mcell].vacc_trig < USHRT_MAX - 1) Mcells[a->mcell].vacc_trig++;
			}
			if (Mcells[a->mcell].move_trig < USHRT_MAX - 1) Mcells[a->mcell].move_trig++;
			if (Mcells[a->mcell].socdist_trig < USHRT_MAX - 1) Mcells[a->mcell].socdist_trig++;
			if (Mcells[a->mcell].keyworkerproph_trig < USHRT_MAX - 1) Mcells[a->mcell].keyworkerproph_trig++;
		//}
#ifndef ABSENTEEISM_PLACE_CLOSURE
#ifdef PLACE_CLOSE_ROUND_HOUSEHOLD
		if(Mcells[a->mcell].place_trig<USHRT_MAX-1) Mcells[a->mcell].place_trig++;
#endif
		if(t>=P.PlaceCloseTimeStart)
			for(j=0;j<P.PlaceTypeNum;j++)
				if((j!=HOTEL_PLACE_TYPE)&&(a->PlaceLinks[j]>=0))
				{
					if((!P.RestrictTreatToTarget)||(Places[j][a->PlaceLinks[j]].country==P.TargetCountry))
					{
						DoPlaceClose(j,a->PlaceLinks[j],ts,tn,0);
#ifndef PLACE_CLOSE_ROUND_HOUSEHOLD
						if(Mcells[Places[j][a->PlaceLinks[j]].mcell].place_trig<USHRT_MAX-1)
						{
#pragma omp critical (place_trig)
							Mcells[Places[j][a->PlaceLinks[j]].mcell].place_trig++;
						}
#endif
					}
				}
#endif
		if(t>=P.TreatTimeStart)
		{
			if((P.DoHouseholds)&&(Households[Hosts[ai].hh].stockpile==1))
			{
				if((P.PrivateTreatPropCases==1)||(ranf_mt(ranf_seed, tn)<P.PrivateTreatPropCases))
				{
					DoPrivateTreatCase(ai,ts,tn);
					if((t<P.TreatTimeStart+P.TreatHouseholdsDuration)&&((P.PrivateTreatPropCaseHouseholds==1)||(ranf_mt(ranf_seed, tn)<P.PrivateTreatPropCaseHouseholds)))
					{
						j1=Households[Hosts[ai].hh].FirstPerson;j2=j1+Households[Hosts[ai].hh].nh;
						for(j=j1;j<j2;j++)
							if(!HOST_TO_BE_TREATED(j)) DoPrivateProph(j,ts,tn);
					}
				}
			}
			else if((P.TreatPropCases==1)||(ranf_mt(ranf_seed, tn)<P.TreatPropCases))
			{
				DoTreatCase(ai,ts,tn);
				if(P.DoHouseholds)
				{
					if((t<P.TreatTimeStart+P.TreatHouseholdsDuration)&&((P.TreatPropCaseHouseholds==1)||(ranf_mt(ranf_seed, tn)<P.TreatPropCaseHouseholds)))
					{
						j1=Households[Hosts[ai].hh].FirstPerson;j2=j1+Households[Hosts[ai].hh].nh;
						for(j=j1;j<j2;j++)
							if(!HOST_TO_BE_TREATED(j)) DoProph(j,ts,tn);
					}
				}
				if(P.DoPlaces)
				{
					if(t<P.TreatTimeStart+P.TreatPlaceGeogDuration)
						for(j=0;j<P.PlaceTypeNum;j++)
							if(a->PlaceLinks[j]>=0)
							{
								if((!P.RestrictTreatToTarget)||(Places[j][a->PlaceLinks[j]].country==P.TargetCountry))
								{
									if(P.DoPlaceGroupTreat)
									{
										if((P.TreatPlaceProbCaseId[j]==1)||(ranf_mt(ranf_seed, tn)<P.TreatPlaceProbCaseId[j]))
										{
											StateT[tn].p_queue[j][StateT[tn].np_queue[j]]=a->PlaceLinks[j];
											StateT[tn].pg_queue[j][StateT[tn].np_queue[j]++]=a->PlaceGroupLinks[j];
										}
									}
									else
									{
										f=0;
#pragma omp critical (starttreat)
										if(!Places[j][a->PlaceLinks[j]].treat) f=Places[j][a->PlaceLinks[j]].treat=1;
										if(f) 
										{
											if((P.TreatPlaceProbCaseId[j]==1)||(ranf_mt(ranf_seed, tn)<P.TreatPlaceProbCaseId[j]))
												StateT[tn].p_queue[j][StateT[tn].np_queue[j]++]=a->PlaceLinks[j];
											else
												Places[j][a->PlaceLinks[j]].treat=0;
										}
									}
								}
							}
				}
			}
		}
		if(P.DoHouseholds)
		{
			if((!P.DoMassVacc)&&(t>=P.VaccTimeStart)&&(State.cumV<P.VaccMaxCourses))
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
			if((t>=P.HQuarantineTimeStart)&&(t<P.HQuarantineTimeStart+P.HQuarantinePolicyDuration))
			{
				j1=Households[Hosts[ai].hh].FirstPerson;j2=j1+Households[Hosts[ai].hh].nh;
				if(!HOST_TO_BE_QUARANTINED(j1))
				{
					Hosts[j1].quar_start_time=ts+((unsigned short int) (P.TimeStepsPerDay*P.HQuarantineHouseDelay));
					k=(ranf_mt(ranf_seed, tn)<P.HQuarantinePropHouseCompliant)?1:0;
					if(k) StateT[tn].cumHQ++;

					Hosts[j1].quar_comply=((k==0)?0:((ranf_mt(ranf_seed, tn)<P.HQuarantinePropIndivCompliant)?1:0));
					if((Hosts[j1].quar_comply)&&(!HOST_ABSENT(j1)))
					{
						if(HOST_AGE_YEAR(j1)>=P.CaseAbsentChildAgeCutoff)
						{
							if(Hosts[j1].PlaceLinks[NUM_PLACE_TYPES_NOAIR-1]>=0) StateT[tn].cumAH++;
						}
						else
							StateT[tn].cumACS++;
					}
					for(j=j1+1;j<j2;j++)
					{
						Hosts[j].quar_start_time=Hosts[j1].quar_start_time;
						Hosts[j].quar_comply=((k==0)?0:((ranf_mt(ranf_seed, tn)<P.HQuarantinePropIndivCompliant)?1:0));
						if((Hosts[j].quar_comply)&&(!HOST_ABSENT(j)))
						{
							if(HOST_AGE_YEAR(j)>=P.CaseAbsentChildAgeCutoff)
							{
								if(Hosts[j].PlaceLinks[NUM_PLACE_TYPES_NOAIR-1]>=0) StateT[tn].cumAH++;
							}
							else
								StateT[tn].cumACS++;
						}
					}
				}
			}
		}
		if((t>=P.CaseIsolationTimeStart)&&(t<P.CaseIsolationTimeStart+P.CaseIsolationPolicyDuration))
		{
			if((P.CaseIsolationProp==1)||(ranf_mt(ranf_seed, tn)<P.CaseIsolationProp))
			{
				Hosts[ai].isolation_start_time=ts;
				if(HOST_ABSENT(ai))
				{
					if(a->absent_stop_time<ts+P.usCaseAbsenteeismDelay+P.usCaseIsolationDuration)
						a->absent_stop_time=ts+P.usCaseAbsenteeismDelay+P.usCaseIsolationDuration;
				}
				else if(P.DoRealSymptWithdrawal)
/* This calculates adult absenteeism from work due to care of isolated children.  */
				{
					Hosts[ai].absent_start_time=ts+P.usCaseIsolationDelay;
					Hosts[ai].absent_stop_time=ts+P.usCaseIsolationDelay+P.usCaseIsolationDuration;
					if(P.DoPlaces)
					{
						if((!HOST_QUARANTINED(ai))&&(Hosts[ai].PlaceLinks[NUM_PLACE_TYPES_NOAIR-1]>=0)&&(HOST_AGE_YEAR(ai)>=P.CaseAbsentChildAgeCutoff))
							StateT[tn].cumAC++;
					}
					if((P.DoHouseholds)&&(P.DoPlaces)&&(HOST_AGE_YEAR(ai)<P.CaseAbsentChildAgeCutoff))
					{
						if(!HOST_QUARANTINED(ai)) StateT[tn].cumACS++;
						if((P.CaseAbsentChildPropAdultCarers==1)||(ranf_mt(ranf_seed, tn)<P.CaseAbsentChildPropAdultCarers))
						{
							j1=Households[Hosts[ai].hh].FirstPerson;j2=j1+Households[Hosts[ai].hh].nh;
							f=0;
							for(j=j1;(j<j2)&&(!f);j++)
								f=((abs(Hosts[j].inf)!=5)&&(HOST_AGE_YEAR(j)>=P.CaseAbsentChildAgeCutoff)&&(Hosts[j].PlaceLinks[NUM_PLACE_TYPES_NOAIR-1]<0));
							if(!f)
							{
								for(j=j1;(j<j2)&&(!f);j++)
									f=((abs(Hosts[j].inf)!=5)&&(HOST_AGE_YEAR(j)>=P.CaseAbsentChildAgeCutoff)&&((HOST_ABSENT(j))||(HOST_QUARANTINED(j))));
								if(!f)
								{
									f=0;k=-1;
									for(j=j1;(j<j2)&(!f);j++)
										if((HOST_AGE_YEAR(j)>=P.CaseAbsentChildAgeCutoff)&&(abs(Hosts[j].inf)!=5)) {k=j;f=1;}
									if(f)
									{
										Hosts[k].absent_start_time=ts+P.usCaseIsolationDelay;
										Hosts[k].absent_stop_time=ts+P.usCaseIsolationDelay+P.usCaseIsolationDuration;
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
	if ((Hosts[ai].detected) && (Hosts[ai].vacc_accept<P.ProbEstablishRing) && ((P.DoRingVaccination&&(t>=P.RingVaccTimeStart)) || (P.DoContactTracing&&(t>=P.ContactTracingTimeStart))))
	{
		nActiveCases = State.cumDC_adunit[Mcells[Hosts[ai].mcell].adunit];// -(State.cumDD_adunit[Mcells[Hosts[ai].mcell].adunit] + State.cumDR_adunit[Mcells[Hosts[ai].mcell].adunit]);
		//first check to see if we've reached the threshold to start contact tracing in the admin unit of the initial case
		if ((t>=P.ContactTracingTimeStart)&&(nActiveCases >= AdUnits[Mcells[Hosts[ai].mcell].adunit].contactTraceCaseThreshold) && (AdUnits[Mcells[Hosts[ai].mcell].adunit].contactTraceThresholdCrossed == 0)) //changed to detected cases
		{
			//mark that we've now reached the threshold and mark the day on which it occurs
			AdUnits[Mcells[Hosts[ai].mcell].adunit].contactTraceThresholdCrossed = 1;
			AdUnits[Mcells[Hosts[ai].mcell].adunit].contactTraceStartDay = (int)t;
		}

		if (((t >= P.RingVaccTimeStart) && (State.cumV < P.VaccMaxCourses)) || AdUnits[Mcells[Hosts[ai].mcell].adunit].contactTraceThresholdCrossed == 1) //modified this to include criteria for having to be past contact tracing threshold
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
					if ((Hosts[ai].PlaceLinks[j] >= 0) && (j!=P.HospPlaceTypeNum))
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
					if ((P.propContactTraced == 0) || (ranf_mt(ranf_seed, tn) < P.propContactTraced))
					{
						StateT[tn].ct_queue[Mcells[Hosts[StateT[tn].ringvacclist[i]].mcell].adunit][StateT[tn].nct_queue[Mcells[Hosts[StateT[tn].ringvacclist[i]].mcell].adunit]++] = StateT[tn].ringvacclist[i];
					}
					//}
				}
			}

			if ((State.cumV < P.VaccMaxCourses) && (t >= P.RingVaccTimeStart))
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
						if (Mcells[Places[P.HospPlaceTypeNum][i].mcell].adunit==Mcells[Hosts[ai].mcell].adunit)
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
					if((Hosts[ai].PlaceLinks[j] >= 0) && (j != P.HospPlaceTypeNum))
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
							if ((Hosts[i].PlaceLinks[j] >= 0) &&(j!=P.HospPlaceTypeNum))
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
						if ((Hosts[ai].PlaceLinks[j] > 0)&&(j!=P.HospPlaceTypeNum))
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
							if((Hosts[ai].PlaceLinks[j] >= 0) && (j != P.HospPlaceTypeNum))
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

void DoCase(int ai,double t,unsigned short int ts,int tn)
{	
	int j,k,f,j1,j2,m,q;
	person *a;
	int age;
	//int *RingVaccList,*RingVaccNonHouseholdList; //added this for my ring vaccination code - ggilani 15/02/17
	//int *RingVaccRingList; //added this to keep track of which ring each contact belongs too - ggilani 29/05/2019
	int currentRing; //to keep track of the current ring
	int nVacc,casePlaceType,i,i1,i2,cnt,k2,k3,l,nAlreadyVacc,ii,vaccflag;
	double propVacc,adjPropToVacc;

	casePlaceType = 0; //added this to initialise casePlaceType memory - ggilani 21/10/19
	currentRing = 0; //to initialise - ggilani 29/10/19

	age=HOST_AGE_GROUP(ai);
	if(age>=NUM_AGE_GROUPS) age=NUM_AGE_GROUPS-1;
	a=Hosts+ai;
	if(a->inf==-1)
	{
		a->inf=-2;
		if(HOST_ABSENT(ai))
		{
			if(a->absent_stop_time<ts+P.usCaseAbsenteeismDelay+P.usCaseAbsenteeismDuration)
				a->absent_stop_time=ts+P.usCaseAbsenteeismDelay+P.usCaseAbsenteeismDuration;
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
						if ((P.SymptPlaceTypeWithdrawalProp[j] == 1) || (ranf_mt(ranf_seed, tn) < P.SymptPlaceTypeWithdrawalProp[j]))
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
								if ((P.CaseAbsentChildPropAdultCarers == 1) || (ranf_mt(ranf_seed, tn) < P.CaseAbsentChildPropAdultCarers))
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
		
		if(HOST_TREATED(ai)) Cells[Hosts[ai].pcell].cumTC++;
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
				Hosts[ai].detect_time = Hosts[ai].latent_time+((int)(P.LatentToSymptDelay / P.TimeStep)); //if detected immediately, set detect_time to be the same time
			}

			if ((P.DoHospitalisation))//&&(t>=P.ETUTimeStart))
			{
				if ((P.OutbreakDetected)&&(ranf_mt(ranf_seed, tn)<P.PropHospSeek))
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

		if(P.DoAdUnits) StateT[tn].cumC_adunit[Mcells[a->mcell].adunit]++;
	}
}

void DoFalseCase(int ai,double t,unsigned short int ts,int tn)
{	
	int j,k,f,j1,j2;
	person *a;

	/* Arguably adult absenteeism to take care of sick kids could be included here, but then output absenteeism would not be 'excess' absenteeism */
	a=Hosts+ai;
	if((P.ControlPropCasesId==1)||(ranf_mt(ranf_seed, tn)<P.ControlPropCasesId))
		{
		if((!P.DoEarlyCaseDiagnosis)||(State.cumDC>=P.PreControlClusterIdCaseThreshold)) StateT[tn].cumDC++;
		DoDetectedCase(ai,t,ts,tn);
		}
	StateT[tn].cumFC++;
}

void DoRecover(int ai, int run, int tn)
{	
	int i,j,x,y;
	person *a;


	a=Hosts+ai;
	if(abs(a->inf)==2)
		{
		i=a->listpos;
		Cells[a->pcell].I--;
		if(P.DoSIS)
			{
			if(Cells[a->pcell].I>0)
				{
				Cells[a->pcell].susceptible[i]=Cells[a->pcell].infected[0];
				Hosts[Cells[a->pcell].susceptible[i]].listpos=i;
				if(Cells[a->pcell].L>0)
					{
					Cells[a->pcell].latent[Cells[a->pcell].L]=Cells[a->pcell].latent[0];
					Hosts[Cells[a->pcell].latent[Cells[a->pcell].L]].listpos=Cells[a->pcell].S+Cells[a->pcell].L;
					}
				}
			else if(Cells[a->pcell].L>0)
				{
				Cells[a->pcell].susceptible[i]=Cells[a->pcell].latent[0];
				Hosts[Cells[a->pcell].susceptible[i]].listpos=i;
				}
			Cells[a->pcell].susceptible[Cells[a->pcell].S]=ai;
			a->listpos=Cells[a->pcell].S;
			Cells[a->pcell].S++;
			Cells[a->pcell].latent++;
			Cells[a->pcell].infected++;
			a->susc*=P.SuscReductionFactorPerInfection;
			a->inf=0;
			a->infector=-1;
			if(P.OutputBitmap)
			{
				if((P.OutputBitmapDetected==0)||((P.OutputBitmapDetected==1)&&(Hosts[ai].detected==1)))
				{

					x=((int) (Households[a->hh].loc_x*P.scalex))-P.bminx;
					y=((int) (Households[a->hh].loc_y*P.scaley))-P.bminy;
					if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
					{
						j=y*bmh->width+x;
						if((j<bmh->imagesize)&&(j>=0))
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
			if(i<Cells[a->pcell].S+Cells[a->pcell].L+Cells[a->pcell].I)
				{
				Cells[a->pcell].susceptible[i]=Cells[a->pcell].infected[Cells[a->pcell].I];
				Hosts[Cells[a->pcell].susceptible[i]].listpos=i;
				}
			a->inf=3*a->inf/abs(a->inf);
			a->listpos=Cells[a->pcell].S+Cells[a->pcell].L+Cells[a->pcell].I;
			Cells[a->pcell].susceptible[a->listpos]=ai;
			if(P.OutputBitmap)
			{
				if((P.OutputBitmapDetected==0)||((P.OutputBitmapDetected==1)&&(Hosts[ai].detected==1)))
				{
					x=((int) (Households[a->hh].loc_x*P.scalex))-P.bminx;
					y=((int) (Households[a->hh].loc_y*P.scaley))-P.bminy;
					if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
					{
						j=y*bmh->width+x;
						if((j<bmh->imagesize)&&(j>=0))
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

void DoDeath(int ai,int tn, int run)
{	
	int j,x,y;
	person *a;

	a=Hosts+ai;
	if((abs(a->inf)==2)||(abs(a->inf)==6)) //added infection status of 6 as well, as this also needs to deal with 'dead' infectious individuals: ggilani 25/10/14
		{
		a->inf=5*a->inf/abs(a->inf);
		Cells[a->pcell].D++;
		Cells[a->pcell].I--;
		if(Cells[a->pcell].I>0)
			{
			Cells[a->pcell].susceptible[a->listpos]=Cells[a->pcell].infected[Cells[a->pcell].I];
			Hosts[Cells[a->pcell].susceptible[a->listpos]].listpos=a->listpos;
			}
		a->listpos=Cells[a->pcell].S+Cells[a->pcell].L+Cells[a->pcell].I;
		Cells[a->pcell].susceptible[a->listpos]=ai;
/*		a->listpos=-1; */
		if(abs(a->inf) != 6) //if if is someone who is infectious after death, we'll do this bit of accounting in IncubRecoverySweep to make sure their death is recorded on the right day, even though they will stay in the infectious list
		{
			StateT[tn].cumD++;
			StateT[tn].cumDa[HOST_AGE_GROUP(ai)]++;
			if (P.DoAdUnits)
			{
				StateT[tn].cumD_adunit[Mcells[a->mcell].adunit]++;
				if(a->detected) StateT[tn].cumDD_adunit[Mcells[a->mcell].adunit]++;
			}
			StateT[tn].cumD_keyworker[a->keyworker]++;
		}
		if(P.OutputBitmap)
		{
			if((P.OutputBitmapDetected==0)||((P.OutputBitmapDetected==1)&&(Hosts[ai].detected==1)))
			{
				x=((int) (Households[a->hh].loc_x*P.scalex))-P.bminx;
				y=((int) (Households[a->hh].loc_y*P.scaley))-P.bminy;
				if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
				{
					j=y*bmh->width+x;
					if((j<bmh->imagesize)&&(j>=0))
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

void DoTreatCase(int ai,unsigned short int ts,int tn)
{
	int j;
	double x,y;

	if(State.cumT<P.TreatMaxCourses)
		{
		if((!HOST_TREATED(ai))&&(Hosts[ai].resist<(MAX_NUM_RESIST_TYPES-1))&&(ranf_mt(ranf_seed, tn)<P.EvolResistTreatMutationRate)) Hosts[ai].resist++;
#ifdef NO_TREAT_PROPH_CASES
		if(!HOST_TO_BE_TREATED(ai))
#endif
			{
			Hosts[ai].treat_start_time=ts+((unsigned short int) (P.TimeStepsPerDay*P.TreatDelayMean));
			Hosts[ai].treat_stop_time=ts+((unsigned short int) (P.TimeStepsPerDay*(P.TreatDelayMean+P.TreatCaseCourseLength)));
			StateT[tn].cumT++;
			if((abs(Hosts[ai].inf)>0)&&(Hosts[ai].inf!=5)) Cells[Hosts[ai].pcell].cumTC++;
			StateT[tn].cumT_resist[Hosts[ai].resist]++;
			StateT[tn].cumT_keyworker[Hosts[ai].keyworker]++;
			if((++Hosts[ai].num_treats)<2) StateT[tn].cumUT++;
			Cells[Hosts[ai].pcell].tot_treat++;
			if(P.DoAdUnits) StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit]++;
			if(P.OutputBitmap)
				{
				x=((int) (Households[Hosts[ai].hh].loc_x*P.scalex))-P.bminx;
				y=((int) (Households[Hosts[ai].hh].loc_y*P.scaley))-P.bminy;
				if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
					{
					j=y*bmh->width+x;
					if((j<bmh->imagesize)&&(j>=0))
						{
#pragma omp atomic
						bmi4[j]++;
						}
					}
				}
			}
		}
}

void DoProph(int ai,unsigned short int ts,int tn)
{
	int j;
	double x,y;

	if(State.cumT<P.TreatMaxCourses)
		{
		if((!HOST_TREATED(ai))&&(Hosts[ai].resist<(MAX_NUM_RESIST_TYPES-1))&&((abs(Hosts[ai].inf)==1)||(abs(Hosts[ai].inf)==2))&&(ranf_mt(tn)<P.EvolResistProphMutationRate)) Hosts[ai].resist++;
		Hosts[ai].treat_start_time=ts+((unsigned short int) (P.TimeStepsPerDay*P.TreatDelayMean));
		Hosts[ai].treat_stop_time=ts+((unsigned short int) (P.TimeStepsPerDay*(P.TreatDelayMean+P.TreatProphCourseLength)));
		StateT[tn].cumT++;
//		StateT[tn].cumT_resist[Hosts[ai].resist]++;
		StateT[tn].cumT_keyworker[Hosts[ai].keyworker]++;
		if((++Hosts[ai].num_treats)<2) StateT[tn].cumUT++;
		if(P.DoAdUnits)	StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit]++;
#pragma omp critical (tot_treat)
		Cells[Hosts[ai].pcell].tot_treat++;
		if(P.OutputBitmap)
			{
			x=((int) (Households[Hosts[ai].hh].loc_x*P.scalex))-P.bminx;
			y=((int) (Households[Hosts[ai].hh].loc_y*P.scaley))-P.bminy;
			if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
				{
				j=y*bmh->width+x;
				if((j<bmh->imagesize)&&(j>=0))
					{
#pragma omp atomic
					bmi4[j]++;
					}
				}
			}
		}
}

void DoPrivateTreatCase(int ai,unsigned short int ts,int tn)
{
	int j;
	double x,y;

	if((!HOST_TREATED(ai))&&(Hosts[ai].resist<(MAX_NUM_RESIST_TYPES-1))&&(ranf_mt(tn)<P.EvolResistTreatMutationRate)) Hosts[ai].resist++;
#ifdef NO_TREAT_PROPH_CASES
	if(!HOST_TO_BE_TREATED(ai))
#endif
		{
		Households[Hosts[ai].hh].stockpile=-1;
		Hosts[ai].treat_start_time=ts+((unsigned short int) (P.TimeStepsPerDay*P.PrivateTreatDelayMean));
		Hosts[ai].treat_stop_time=ts+((unsigned short int) (P.TimeStepsPerDay*(P.PrivateTreatDelayMean+P.TreatCaseCourseLength)));
		StateT[tn].cumTP++;
		if((abs(Hosts[ai].inf)>0)&&(Hosts[ai].inf!=5))Cells[Hosts[ai].pcell].cumTC++;
		StateT[tn].cumT_resist[Hosts[ai].resist]++;
		StateT[tn].cumT_keyworker[Hosts[ai].keyworker]++;
		if((++Hosts[ai].num_treats)<2) StateT[tn].cumUT++;
		Cells[Hosts[ai].pcell].tot_treat++;
		if(P.DoAdUnits) StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit]++;
		if(P.OutputBitmap)
			{
			x=((int) (Households[Hosts[ai].hh].loc_x*P.scalex))-P.bminx;
			y=((int) (Households[Hosts[ai].hh].loc_y*P.scaley))-P.bminy;
			if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
				{
				j=y*bmh->width+x;
				if((j<bmh->imagesize)&&(j>=0))
					{
#pragma omp atomic
					bmi4[j]++;
					}
				}
			}
		}
}

void DoPrivateProph(int ai,unsigned short int ts,int tn)
{
	int j;
	double x,y;

	if((!HOST_TREATED(ai))&&(Hosts[ai].resist<(MAX_NUM_RESIST_TYPES-1))&&((abs(Hosts[ai].inf)==1)||(abs(Hosts[ai].inf)==2))&&(ranf_mt(tn)<P.EvolResistProphMutationRate)) Hosts[ai].resist++;
	Hosts[ai].treat_start_time=ts+((unsigned short int) (P.TimeStepsPerDay*P.PrivateTreatDelayMean));
	Hosts[ai].treat_stop_time=ts+((unsigned short int) (P.TimeStepsPerDay*(P.PrivateTreatDelayMean+P.TreatProphCourseLength)));
	StateT[tn].cumTP++;
	//StateT[tn].cumT_resist[Hosts[ai].resist]++;
	StateT[tn].cumT_keyworker[Hosts[ai].keyworker]++;
	if((++Hosts[ai].num_treats)<2) StateT[tn].cumUT++;
	if(P.DoAdUnits)	StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit]++;
#pragma omp critical (tot_treat)
	Cells[Hosts[ai].pcell].tot_treat++;
	if(P.OutputBitmap)
		{
		x=((int) (Households[Hosts[ai].hh].loc_x*P.scalex))-P.bminx;
		y=((int) (Households[Hosts[ai].hh].loc_y*P.scaley))-P.bminy;
		if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
			{
			j=y*bmh->width+x;
			if((j<bmh->imagesize)&&(j>=0))
				{
#pragma omp atomic
				bmi4[j]++;
				}
			}
		}
}
void DoProphNoDelay(int ai,unsigned short int ts,int tn,int nc)
{
	int j;
	double x,y;

	if(State.cumT<P.TreatMaxCourses)
		{
		if((!HOST_TREATED(ai))&&(Hosts[ai].resist<(MAX_NUM_RESIST_TYPES-1))&&((abs(Hosts[ai].inf)==1)||(abs(Hosts[ai].inf)==2))&&(ranf_mt(tn)<P.EvolResistProphMutationRate)) Hosts[ai].resist++;
		Hosts[ai].treat_start_time=ts;
		Hosts[ai].treat_stop_time=ts+((unsigned short int) (P.TimeStepsPerDay*P.TreatProphCourseLength*nc));
		StateT[tn].cumT+=nc;
		//StateT[tn].cumT_resist[Hosts[ai].resist]+=nc;
		StateT[tn].cumT_keyworker[Hosts[ai].keyworker]+=nc;
		if((++Hosts[ai].num_treats)<2) StateT[tn].cumUT++;
		if(P.DoAdUnits) StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit]+=nc;
#pragma omp critical (tot_treat)
		Cells[Hosts[ai].pcell].tot_treat++;
		if(P.OutputBitmap)
			{
			x=((int) (Households[Hosts[ai].hh].loc_x*P.scalex))-P.bminx;
			y=((int) (Households[Hosts[ai].hh].loc_y*P.scaley))-P.bminy;
			if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
				{
				j=y*bmh->width+x;
				if((j<bmh->imagesize)&&(j>=0))
					{
#pragma omp atomic
					bmi4[j]++;
					}
				}
			}
		}
}


void DoPlaceClose(int i,int j,unsigned short int ts,int tn,int DoAnyway)
{
	int k,ai,j1,j2,l,f,m,f2;
	unsigned short trig;
#ifdef ABSENTEEISM_PLACE_CLOSURE
	unsigned short int t_old,t_new;
#endif


	f2=0;
/*	if((j<0)||(j>=P.Nplace[i]))
		fprintf(stderr,"** %i %i *\n",i,j);
	else
*/
#ifdef ABSENTEEISM_PLACE_CLOSURE
	t_new=ts/P.TimeStepsPerDay;
#endif
	trig=0;
#pragma omp critical (closeplace)
	  {
	  if(!PLACE_TO_BE_CLOSED(i,j))
		{
		if((!DoAnyway)&&(Places[i][j].control_trig<USHRT_MAX-2))
			{
#ifdef ABSENTEEISM_PLACE_CLOSURE
			t_old=Places[i][j].AbsentLastUpdateTime;
			if(t_new>=t_old+MAX_ABSENT_TIME)
				for(l=0;l<MAX_ABSENT_TIME;l++) Places[i][j].Absent[l]=0;
			else
				for(l=t_old;l<t_new;l++) Places[i][j].Absent[l%MAX_ABSENT_TIME]=0;
			for(l=t_new;l<t_new+P.usCaseAbsenteeismDuration/P.TimeStepsPerDay;l++) Places[i][j].Absent[l%MAX_ABSENT_TIME]++;
			trig=Places[i][j].Absent[t_new%MAX_ABSENT_TIME];
			Places[i][j].AbsentLastUpdateTime=t_new;
			if((P.PlaceCloseByAdminUnit)&&(P.PlaceCloseAdunitPlaceTypes[i]>0)
				&&(((double) trig)/((double) Places[i][j].n)>P.PlaceCloseCasePropThresh))
				{
				//fprintf(stderr,"** %i %i %i %i %lg ## ",i,j,(int) Places[i][j].control_trig, (int) Places[i][j].n,P.PlaceCloseCasePropThresh);
				k=Mcells[Places[i][j].mcell].adunit;
				if(AdUnits[k].place_close_trig<USHRT_MAX-1) AdUnits[k].place_close_trig++;
				}
#else
			trig=(++Places[i][j].control_trig);
			if((P.PlaceCloseByAdminUnit)&&(P.PlaceCloseAdunitPlaceTypes[i]>0)
				&&(((double) Places[i][j].control_trig)/((double) Places[i][j].n)>P.PlaceCloseCasePropThresh))
				{
				//fprintf(stderr,"** %i %i %i %i %lg ## ",i,j,(int) Places[i][j].control_trig, (int) Places[i][j].n,P.PlaceCloseCasePropThresh);
				k=Mcells[Places[i][j].mcell].adunit;
				if(AdUnits[k].place_close_trig<USHRT_MAX-1) AdUnits[k].place_close_trig++;
				}
#endif
			}
		if(Places[i][j].control_trig<USHRT_MAX-1)
			{
			if(P.PlaceCloseFracIncTrig>0)
				k=(((double) trig)/((double) Places[i][j].n)>P.PlaceCloseFracIncTrig);
			else
				k=(((int) trig)>=P.PlaceCloseIncTrig);
			if(((!P.PlaceCloseByAdminUnit)&&(k))||(DoAnyway))
				{
				if(P.DoPlaceCloseOnceOnly)
					Places[i][j].control_trig=USHRT_MAX-1;  // Places only close once
				else
					Places[i][j].control_trig=0;
				if((P.PlaceCloseEffect[i]==0)||(ranf_mt(tn)>=P.PlaceCloseEffect[i]))
					{
					Places[i][j].close_start_time=ts+((unsigned short int) (P.TimeStepsPerDay*P.PlaceCloseDelayMean));
					Places[i][j].close_end_time=ts+((unsigned short int) (P.TimeStepsPerDay*(P.PlaceCloseDelayMean+P.PlaceCloseDuration)));
					f2=1;
					}
				else
					Places[i][j].close_start_time=Places[i][j].close_end_time=ts+((unsigned short int) (P.TimeStepsPerDay*(P.PlaceCloseDelayMean+P.PlaceCloseDuration)));
				}
			}
		}
	  }
	if(f2)
		{
		if(P.DoRealSymptWithdrawal)
			for(k=0;k<Places[i][j].n;k++)
					{
					ai=Places[i][j].members[k];
					Hosts[ai].absent_start_time=Places[i][j].close_start_time;
					Hosts[ai].absent_stop_time=Places[i][j].close_end_time;
					if((HOST_AGE_YEAR(ai)>=P.CaseAbsentChildAgeCutoff)&&(Hosts[ai].PlaceLinks[NUM_PLACE_TYPES_NOAIR-1]>=0)) StateT[tn].cumAPC++;
					if((HOST_AGE_YEAR(ai)<P.CaseAbsentChildAgeCutoff)&&(!HOST_ABSENT(ai))&&(!HOST_QUARANTINED(ai)))
						{
						StateT[tn].cumAPCS++;
						if((P.CaseAbsentChildPropAdultCarers==1)||(ranf_mt(tn)<P.CaseAbsentChildPropAdultCarers))
							{
							j1=Households[Hosts[ai].hh].FirstPerson;j2=j1+Households[Hosts[ai].hh].nh;
							if((j1<0)||(j2>=P.N)) fprintf(stderr,"++ %i %i %i (%i %i %i)##  ",ai,j1,j2,i,j,k);
							f=0;
							for(l=j1;(l<j2)&&(!f);l++)
								f=((abs(Hosts[l].inf)!=5)&&(HOST_AGE_YEAR(l)>=P.CaseAbsentChildAgeCutoff)&&
								((Hosts[l].PlaceLinks[NUM_PLACE_TYPES_NOAIR-1]<0)||(HOST_ABSENT(l))||(HOST_QUARANTINED(l)))); 
							if(!f)
								{
								for(l=j1;(l<j2)&&(!f);l++)
									if((HOST_AGE_YEAR(l)>=P.CaseAbsentChildAgeCutoff)&&(abs(Hosts[l].inf)!=5)) {m=l;f=1;}
								if(f)
									{
									Hosts[m].absent_start_time=Places[i][j].close_start_time;
									Hosts[m].absent_stop_time=Places[i][j].close_end_time;
									StateT[tn].cumAPA++;
									}
								}
							}
						}
					}
		}
}


int DoVacc(int ai,int ts,int ringflag)
{
	int j;
	double x,y;
	int age;

	age = HOST_AGE_GROUP(ai);

	if((!P.DoDistributionVaccination)&&(State.cumV>=P.VaccMaxCourses))
		return 2;
	else if((HOST_TO_BE_VACCED(ai)&&(Hosts[ai].vacc_start_time>0))||(Hosts[ai].inf<-1)||(Hosts[ai].inf>=5))
		return 1;
	else
		{
		if (Hosts[ai].keyworker && Hosts[ai].vacc_start_time <= 0)
		{
			Hosts[ai].revacc = 1;
		}
		Hosts[ai].vacc_start_time=ts+((int) (P.TimeStepsPerDay*P.VaccDelayMean));
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
		else if (ringflag==0) //for geographic vaccination
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
		if(P.OutputBitmap)
			{
			x=((int) (Households[Hosts[ai].hh].loc_x*P.scalex))-P.bminx;
			y=((int) (Households[Hosts[ai].hh].loc_y*P.scaley))-P.bminy;
			if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
				{
				j=y*bmh->width+x;
				if((j<bmh->imagesize)&&(j>=0))
					{
#pragma omp atomic
					bmi4[j]++;
					}
				}
			}
		}
	return 0;
}

#ifdef FRESSCA
int DoDistribVacc(int ai,unsigned short int ts)
{
	int j;
	double x,y;

	if((Hosts[ai].inf<-1)||(Hosts[ai].inf>=5))
		return 1;
	else
		{
		DoImmune(ai);
#pragma omp critical (tot_vacc)
		Cells[Hosts[ai].pcell].tot_vacc++;
		if(P.OutputBitmap)
			{
			x=((int) (Households[Hosts[ai].hh].loc_x*P.scalex))-P.bminx;
			y=((int) (Households[Hosts[ai].hh].loc_y*P.scaley))-P.bminy;
			if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
				{
				j=y*bmh->width+x;
				if((j<bmh->imagesize)&&(j>=0))
					{
#pragma omp atomic
					bmi4[j]++;
					}
				}
			}
		}
	return 0;
}
#endif

void DoVaccNoDelay(int ai,int ts)
	{
	int j;
	double x,y;

	if((State.cumV<P.VaccMaxCourses)&&(!HOST_TO_BE_VACCED(ai))&&(Hosts[ai].inf>=-1)&&(Hosts[ai].inf<5)) //changing cumulative vaccine doses to VG - to separate ring from geographic counts
		{
		Hosts[ai].vacc_start_time=ts;
#pragma omp critical (state_cumV) //changed to VG
		State.cumV++; //changed to VG
#pragma omp critical (state_cumV_daily)
		if(P.VaccDosePerDay>=0)
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
		if(P.OutputBitmap)
			{
			x=((int) (Households[Hosts[ai].hh].loc_x*P.scalex))-P.bminx;
			y=((int) (Households[Hosts[ai].hh].loc_y*P.scaley))-P.bminy;
			if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
				{
				j=y*bmh->width+x;
				if((j<bmh->imagesize)&&(j>=0))
					{
#pragma omp atomic
					bmi4[j]++;
					}
				}
			}
		}
	}

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
	int tn,i,j,k,l,m,p;
	double total_flow,flow;
	int cl_from,cl_to,cl_from_mcl,cl_to_mcl,mcl_from,mcl_to;
	double pop_dens_from[MAX_ADUNITS],pop_dens_to[MAX_ADUNITS];

#pragma omp parallel for private(tn,i,j,k,l,m,p,total_flow,mcl_from,mcl_to,cl_from,cl_to,cl_from_mcl,cl_to_mcl,pop_dens_from,pop_dens_to,flow) schedule(static) //reduction(+:s,t2)
	for(tn=0;tn<P.NumThreads;tn++)
	{
		for(i=tn;i<P.NCP;i+=P.NumThreads)
		{
			//reset pop density matrix to zero
			for(k=0;k<P.NumAdunits;k++)
			{
				pop_dens_from[k]=0.0;
			}

			//find index of cell from which flow travels
			cl_from=CellLookup[i]-Cells;
			cl_from_mcl=(cl_from/P.nch)*P.NMCL*P.nmch+(cl_from%P.nch)*P.NMCL;

			//loop over microcells in these cells to find populations in each admin unit and so flows
			for(k=0;k<P.NMCL;k++)
			{
				for(l=0;l<P.NMCL;l++)
				{
					//get index of microcell
					mcl_from=cl_from_mcl+l+k*P.nmch;
					if(Mcells[mcl_from].n>0)
					{
						//get proportion of each population of cell that exists in each admin unit
						pop_dens_from[Mcells[mcl_from].adunit]+=(((double) Mcells[mcl_from].n)/((double) Cells[cl_from].n));
					}
				}
			}

			for(j=i;j<P.NCP;j++)
			{
				//reset pop density matrix to zero
				for(m=0;m<P.NumAdunits;m++)
				{
					pop_dens_to[m]=0.0;
				}

				//find index of cell which flow travels to
				cl_to=CellLookup[j]-Cells;
				cl_to_mcl=(cl_to/P.nch)*P.NMCL*P.nmch+(cl_to%P.nch)*P.NMCL;
				//calculate distance and kernel between the cells
				//total_flow=Cells[cl_from].max_trans[j]*Cells[cl_from].n*Cells[cl_to].n;
				if(j==0)
				{
					total_flow=Cells[cl_from].cum_trans[j]*Cells[cl_from].n;
				}
				else
				{
					total_flow=(Cells[cl_from].cum_trans[j]-Cells[cl_from].cum_trans[j-1])*Cells[cl_from].n;
				}

				//loop over microcells within destination cell
				for(m=0;m<P.NMCL;m++)
				{
					for(p=0;p<P.NMCL;p++)
					{
						//get index of microcell
						mcl_to=cl_to_mcl+p+m*P.nmch;
						if(Mcells[mcl_to].n>0)
						{
							//get proportion of each population of cell that exists in each admin unit
							pop_dens_to[Mcells[mcl_to].adunit]+=(((double) Mcells[mcl_to].n)/((double) Cells[cl_to].n));
						}
					}
				}

				for(m=0;m<P.NumAdunits;m++)
				{
					for(p=0;p<P.NumAdunits;p++)
					{
						if(m!=p)
						{
							if(AdUnits[m].id/P.CountryDivisor==AdUnits[p].id/P.CountryDivisor)
							{
								flow=total_flow*pop_dens_from[m]*pop_dens_to[p];
							}
							else
							{
								flow=total_flow*pop_dens_from[m]*pop_dens_to[p]*P.PropCrossBorderInf;
							}
							StateT[tn].origin_dest[m][p]+=flow;
							StateT[tn].origin_dest[p][m]+=flow;
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
	for(i=0;i<P.NumAdunits;i++)
	{
		for(j=0;j<P.NumAdunits;j++)
		{
			for(k=0;k<P.NumThreads;k++)
			{
				AdUnits[i].origin_dest[j]+=StateT[k].origin_dest[i][j];
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


double ExpKernel(double r2)
{
	return exp(-sqrt(r2)/P.KernelScale);
}

double PowerKernel(double r2)
{
	double t;

	t=-P.KernelShape*log(sqrt(r2)/P.KernelScale+1);

	return (t<-690)?0:exp(t);
}

double PowerKernelB(double r2)
{
	double t;

	t=0.5*P.KernelShape*log(r2/(P.KernelScale*P.KernelScale));

	return (t>690)?0:(1/(exp(t)+1));
}

double PowerKernelUS(double r2)
	{
	double t;

	t=log(sqrt(r2)/P.KernelScale+1);

	return (t<-690)?0:(exp(-P.KernelShape*t)+P.KernelP3*exp(-P.KernelP4*t))/(1+P.KernelP3);
	}


double GaussianKernel(double r2)
{
	return exp(-r2/(P.KernelScale*P.KernelScale));
}

double StepKernel(double r2)
{
	return (r2>P.KernelScale*P.KernelScale)?0:1;
}

double PowerExpKernel(double r2)
{
	double d,t,t2;

	d=sqrt(r2);
	t=-P.KernelShape*log(d/P.KernelScale+1);

	return (t<-690)?0:exp(t-pow(d/P.KernelP3,P.KernelP4));
}

double numKernel(double r2)
{
	double t,s;

	t=r2/P.KernelDelta;
	if(t>NKR) 
	{
		fprintf(stderr,"** %lg  %lg  %lg**\n",r2,P.KernelDelta,t);
		ERR_CRITICAL("r too large in NumKernel\n");
	}
	s=t*NK_HR;
	if(s<NKR)
		{
		t=s-floor(s);
		t=(1-t)*nKernelHR[(int) s]+t*nKernelHR[(int) (s+1)];
		}
	else
		{
		s=t-floor(t);
		t=(1-s)*nKernel[(int) t]+s*nKernel[(int) (t+1)];
		}
	return t;
}


double dist2UTM(double x1,double y1,double x2,double y2)
	{
	double x,y,cy1,cy2,yt,xi,yi;

	x=fabs(x1-x2)/2;
	y=fabs(y1-y2)/2;
	xi=floor(x);
	yi=floor(y);
	x-=xi;
	y-=yi;
	x=(1-x)*sinx[(int) xi]+x*sinx[((int) xi)+1];
	x=x*x;
	y=(1-y)*sinx[(int) yi]+y*sinx[((int) yi)+1];
	y=y*y;
	yt=fabs(y1+P.SpatialBoundingBox[1]);
	yi=floor(yt);
	cy1=yt-yi;
	cy1=(1-cy1)*cosx[((int) yi)]+cy1*cosx[((int) yi)+1];
	yt=fabs(y2+P.SpatialBoundingBox[1]);
	yi=floor(yt);
	cy2=yt-yi;
	cy2=(1-cy2)*cosx[((int) yi)]+cy2*cosx[((int) yi)+1];
	x=fabs(1000*(y+x*cy1*cy2));
	xi=floor(x);
	x-=xi;
	y=(1-x)*asin2sqx[((int) xi)]+x*asin2sqx[((int) xi)+1];
	return 4*EARTHRADIUS*EARTHRADIUS*y;	
	}

double dist2(person *a, person *b)
{
	double x,y;

	if(P.DoUTM_coords)
		return dist2UTM(Households[a->hh].loc_x,Households[a->hh].loc_y,Households[b->hh].loc_x,Households[b->hh].loc_y);
	else
		{
		x=fabs(Households[a->hh].loc_x-Households[b->hh].loc_x);
		y=fabs(Households[a->hh].loc_y-Households[b->hh].loc_y);
		if(P.DoPeriodicBoundaries)
			{
			if(x>P.width*0.5) x=P.width-x;
			if(y>P.height*0.5) y=P.height-y;
			}
		return x*x+y*y;
		}
}


double dist2_cc(cell *a, cell *b)
{
	double x,y;
	int l,m;

	l=(int) (a-Cells);
	m=(int) (b-Cells);
	if(P.DoUTM_coords)
		return dist2UTM(P.cwidth*fabs((double) (l/P.nch)),P.cheight*fabs((double) (l%P.nch)),
			P.cwidth*fabs((double) (m/P.nch)),P.cheight*fabs((double) (m%P.nch)));
	else
		{
		x=P.cwidth*fabs((double) (l/P.nch-m/P.nch));
		y=P.cheight*fabs((double) (l%P.nch-m%P.nch));
		if(P.DoPeriodicBoundaries)
			{
			if(x>P.width*0.5) x=P.width-x;
			if(y>P.height*0.5) y=P.height-y;
			}
		return x*x+y*y;
		}
}

double dist2_cc_min(cell *a, cell *b)
{
	double x,y;
	int l,m,i,j;
	
	l=(int) (a-Cells);
	m=(int) (b-Cells);
	i=l;j=m;
	if(P.DoUTM_coords)
		{
		if(P.cwidth*((double) abs(m/P.nch-l/P.nch))>PI)
			{
			if(m/P.nch>l/P.nch)
				j+=P.nch;
			else if(m/P.nch<l/P.nch)
				i+=P.nch;
			}
		else
			{
			if(m/P.nch>l/P.nch)
				i+=P.nch;
			else if(m/P.nch<l/P.nch)
				j+=P.nch;
			}
		if(m%P.nch>l%P.nch)
			i++;
		else if(m%P.nch<l%P.nch)
			j++;
		return dist2UTM(P.cwidth*fabs((double) (i/P.nch)),P.cheight*fabs((double) (i%P.nch)),
			P.cwidth*fabs((double) (j/P.nch)),P.cheight*fabs((double) (j%P.nch)));
		}
	else
		{
		if((P.DoPeriodicBoundaries)&&(P.cwidth*((double) abs(m/P.nch-l/P.nch))>P.width*0.5))
			{
			if(m/P.nch>l/P.nch)
				j+=P.nch;
			else if(m/P.nch<l/P.nch)
				i+=P.nch;
			}
		else
			{
			if(m/P.nch>l/P.nch)
				i+=P.nch;
			else if(m/P.nch<l/P.nch)
				j+=P.nch;
			}
		if((P.DoPeriodicBoundaries)&&(P.height*((double) abs(m%P.nch-l%P.nch))>P.height*0.5))
			{
			if(m%P.nch>l%P.nch)
				j++;
			else if(m%P.nch<l%P.nch)
				i++;
			}
		else
			{
			if(m%P.nch>l%P.nch)
				i++;
			else if(m%P.nch<l%P.nch)
				j++;
			}
		x=P.cwidth*fabs((double) (i/P.nch-j/P.nch));
		y=P.cheight*fabs((double) (i%P.nch-j%P.nch));
		if(P.DoPeriodicBoundaries)
			{
			if(x>P.width*0.5) x=P.width-x;
			if(y>P.height*0.5) y=P.height-y;
			}
		return x*x+y*y;
		}
}


double dist2_mm(microcell *a, microcell *b)
	{
	double x,y;
	int l,m;

	l=(int) (a-Mcells);
	m=(int) (b-Mcells);
	if(P.DoUTM_coords)
		return dist2UTM(P.mcwidth*fabs((double) (l/P.nmch)),P.mcheight*fabs((double) (l%P.nmch)),
			P.mcwidth*fabs((double) (m/P.nmch)),P.mcheight*fabs((double) (m%P.nmch)));
	else
		{
		x=P.mcwidth*fabs((double) (l/P.nmch-m/P.nmch));
		y=P.mcheight*fabs((double) (l%P.nmch-m%P.nmch));
		if(P.DoPeriodicBoundaries)
			{
			if(x>P.width*0.5) x=P.width-x;
			if(y>P.height*0.5) y=P.height-y;
			}
		return x*x+y*y;
		}
	}


//added function to find minimum distance
double dist2_mm_min(microcell *a, microcell *b)
{
	double x,y;
	int l,m,i,j;
	
	l=(int) (a-Mcells);
	m=(int) (b-Mcells);
	i=l;j=m;
	if(P.DoUTM_coords)
		{
		if(P.mcwidth*((double) abs(m/P.nmch-l/P.nmch))>PI)
			{
			if(m/P.nmch>l/P.nmch)
				j+=P.nmch;
			else if(m/P.nmch<l/P.nmch)
				i+=P.nmch;
			}
		else
			{
			if(m/P.nmch>l/P.nmch)
				i+=P.nmch;
			else if(m/P.nmch<l/P.nmch)
				j+=P.nmch;
			}
		if(m%P.nmch>l%P.nmch)
			i++;
		else if(m%P.nmch<l%P.nmch)
			j++;
		return dist2UTM(P.mcwidth*fabs((double) (i/P.nmch)),P.mcheight*fabs((double) (i%P.nmch)),
			P.mcwidth*fabs((double) (j/P.nmch)),P.mcheight*fabs((double) (j%P.nmch)));
		}
	else
		{
		if((P.DoPeriodicBoundaries)&&(P.mcwidth*((double) abs(m/P.nmch-l/P.nmch))>P.width*0.5))
			{
			if(m/P.nmch>l/P.nmch)
				j+=P.nmch;
			else if(m/P.nmch<l/P.nmch)
				i+=P.nmch;
			}
		else
			{
			if(m/P.nmch>l/P.nmch)
				i+=P.nmch;
			else if(m/P.nmch<l/P.nmch)
				j+=P.nmch;
			}
		if((P.DoPeriodicBoundaries)&&(P.height*((double) abs(m%P.nmch-l%P.nmch))>P.height*0.5))
			{
			if(m%P.nmch>l%P.nmch)
				j++;
			else if(m%P.nmch<l%P.nmch)
				i++;
			}
		else
			{
			if(m%P.nmch>l%P.nmch)
				i++;
			else if(m%P.nmch<l%P.nmch)
				j++;
			}
		x=P.mcwidth*fabs((double) (i/P.nmch-j/P.nmch));
		y=P.mcheight*fabs((double) (i%P.nmch-j%P.nmch));
		if(P.DoPeriodicBoundaries)
			{
			if(x>P.width*0.5) x=P.width-x;
			if(y>P.height*0.5) y=P.height-y;
			}
		return x*x+y*y;
		}
}

double dist2_raw(double ax,double ay,double bx,double by)
{
	double x,y;

	if(P.DoUTM_coords)
		return dist2UTM(ax,ay,bx,by);
	else
		{
		x=fabs(ax-bx);
		y=fabs(ay-by);
		if(P.DoPeriodicBoundaries)
			{
			if(x>P.width*0.5) x=P.width-x;
			if(y>P.height*0.5) y=P.height-y;
			}
		return x*x+y*y;
		}
}

int GetInputParameter(FILE *dat,FILE *dat2,char *SItemName,char *ItemType,void *ItemPtr, int NumItem,int NumItem2,int Offset)
{
	int FindFlag;
	
	FindFlag=GetInputParameter2(dat,dat2,SItemName,ItemType,ItemPtr,NumItem,NumItem2,Offset);
	if(!FindFlag)
		{
		fprintf(stderr,"\nUnable to find parameter `%s' in input file. Aborting program...\n",SItemName);
		exit(1);
		}
	return FindFlag;
}

int GetInputParameter2(FILE *dat,FILE *dat2,char *SItemName,char *ItemType,void *ItemPtr, int NumItem,int NumItem2,int Offset)
{
	int FindFlag=0;

	if(dat2) FindFlag=GetInputParameter3(dat2,SItemName,ItemType,ItemPtr,NumItem,NumItem2,Offset);
	if(!FindFlag)
		FindFlag=GetInputParameter3(dat,SItemName,ItemType,ItemPtr,NumItem,NumItem2,Offset);
	return FindFlag;
}

/*
	Reads a string (as per fscanf %s).
	Returns true if it succeeds, false on EOF, and does not return on error.
*/
bool readString(const char* SItemName, FILE* dat, char* buf) {
	int r = fscanf(dat, "%s", buf);
	if (r == 1) {
		return true;
	}
	else if (r == EOF) {
		if (ferror(dat)) {
			ERR_CRITICAL("fscanf failed for %s: %s.\n", SItemName, strerror(errno));
		}
		else {
			// EOF
			return false;
		}
	}
	else {
		ERR_CRITICAL("Unexpected fscanf result %d for %s.\n", r, SItemName);
	}
}

int GetInputParameter3(FILE* dat, const char* SItemName, const char* ItemType, void* ItemPtr, int NumItem, int NumItem2, int Offset)
{
	char match[10000] = "", ReadItemName[10000] = "", ItemName[10000];
	int FindFlag = 0, EndString, CurPos, i, j, n;

	n = 0;
	fseek(dat, 0, 0);
	sprintf(ItemName, "[%s]", SItemName);
	while (!FindFlag)
	{
		if (!readString(SItemName, dat, match)) return 0;
		FindFlag = (!strncmp(match, ItemName, strlen(match)));
		if (FindFlag)
		{
			CurPos = ftell(dat);
			strcpy(ReadItemName, match);
			EndString = (match[strlen(match) - 1] == ']');
			while ((!EndString) && (FindFlag))
			{
				if (!readString(SItemName, dat, match)) return 0;
				strcat(ReadItemName, " ");
				strcat(ReadItemName, match);
				FindFlag = (!strncmp(ReadItemName, ItemName, strlen(ReadItemName)));
				EndString = (ReadItemName[strlen(ReadItemName) - 1] == ']');
			}
			if (!EndString)
			{
				fseek(dat, CurPos, 0);
				FindFlag = 0;
			}
		}
	}
	if (FindFlag)
	{
		FindFlag = 0;
		if (!strcmp(ItemType, "%lf"))	n = 1;
		else if (!strcmp(ItemType, "%i"))	n = 2;
		else if (!strcmp(ItemType, "%s"))	n = 3;
		if (NumItem2 < 2)
		{
			if (NumItem == 1)
			{
				if (fscanf(dat, "%s", match) != 1) { ERR_CRITICAL("fscanf failed for %s\n", SItemName); }
				
				else if ((match[0] == '#') && (match[1] == '1') && (match[2] == '0'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP10;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP10;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '1') && (match[2] == '1'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP11;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP11;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '1'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP1;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP1;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '2'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP2;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP2;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '3'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP3;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP3;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '4'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP4;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP4;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '5'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP5;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP5;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '6'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP6;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP6;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '7'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP7;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP7;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '8'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP8;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP8;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '9'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP9;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP9;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] != '[') && (!feof(dat)))
				{
					FindFlag++;
					if (n == 1)
						sscanf(match, "%lf", (double*)ItemPtr);
					else if (n == 2)
						sscanf(match, "%i", (int*)ItemPtr);
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
			}
			else
			{
				for (CurPos = 0; CurPos < NumItem; CurPos++)
				{
					if (fscanf(dat, "%s", match) != 1) { ERR_CRITICAL("fscanf failed for %s\n", SItemName); }
					if ((match[0] != '[') && (!feof(dat)))
					{
						FindFlag++;
						if (n == 1)
							sscanf(match, "%lf", ((double*)ItemPtr) + CurPos + Offset);
						else if (n == 2)
							sscanf(match, "%i", ((int*)ItemPtr) + CurPos + Offset);
						else if (n == 3)
							sscanf(match, "%s", *(((char**)ItemPtr) + CurPos + Offset));
					}
					else
						CurPos = NumItem;
				}
			}
		}
		else
		{
			for (j = 0; j < NumItem; j++)
			{ //added these braces
				for (i = 0; i < NumItem2; i++)
				{
					if (fscanf(dat, "%s", match) != 1) { ERR_CRITICAL("fscanf failed for %s\n", SItemName); }
					if ((match[0] != '[') && (!feof(dat)))
					{
						FindFlag++;
						if (n == 1)
							sscanf(match, "%lf", ((double**)ItemPtr)[j + Offset] + i + Offset); //changed from [j+Offset]+i+Offset to +j+Offset+i, as ItemPtr isn't an array - 01/10: changed it back
						else
							sscanf(match, "%i", ((int**)ItemPtr)[j + Offset] + i + Offset);
					}
					else
					{
						i = NumItem2;
						j = NumItem;
					}
				}
				//Offset=Offset+(NumItem2-1); //added this line to get the correct offset in address position when incrementing j
			} //added these braces
		}
	}
	//	fprintf(stderr,"%s\n",SItemName);
	return FindFlag;

//	char match[10000]="",ReadItemName[10000]="",ItemName[10000];
//	int FindFlag=0,EndString,CurPos,i,j,n;
//
//	fseek(dat,0,0);
//	sprintf(ItemName,"[%s]",SItemName);
//	while((!FindFlag)&&(!feof(dat)))
//		{
//		fscanf(dat,"%s",match);
//		FindFlag=(!strncmp(match,ItemName,strlen(match)));
//		if(FindFlag)
//			{
//			CurPos=ftell(dat);
//			strcpy(ReadItemName,match);
//			EndString=(match[strlen(match)-1]==']');
//			while((!feof(dat))&&(!EndString)&&(FindFlag))
//				{
//				fscanf(dat,"%s",match);
//				strcat(ReadItemName," ");
//				strcat(ReadItemName,match);
//				FindFlag=(!strncmp(ReadItemName,ItemName,strlen(ReadItemName)));
//				EndString=(ReadItemName[strlen(ReadItemName)-1]==']');
//				}
//			if(!EndString)
//				{
//				fseek(dat,CurPos,0);
//				FindFlag=0;
//				}
//			}
//		}
//	if(FindFlag)
//		{
//		FindFlag=0;
//		if(!strcmp(ItemType,"%lg"))
//			n=1;
//		else if(!strcmp(ItemType,"%i"))
//			n=2;
//		else if(!strcmp(ItemType,"%s"))
//			n=3;
//		if(NumItem2<2)
//			{
//			if(NumItem==1)
//				{
//				fscanf(dat,"%s",match);
//				if((match[0]!='[')&&(!feof(dat)))
//					{
//					FindFlag++;
//					if(n==1)
//						sscanf(match,"%lg",(double *)ItemPtr);
//					else if(n==2)
//						sscanf(match,"%i",(int *)ItemPtr);
//					else if(n==3)
//						sscanf(match,"%s",(char *)ItemPtr);
//					}
//				}
//			else
//				{	
//				for(CurPos=0;CurPos<NumItem;CurPos++)
//					{
//					fscanf(dat,"%s",match);
//					if((match[0]!='[')&&(!feof(dat)))
//						{
//						FindFlag++;
//						if(n==1)
//							sscanf(match,"%lg",((double *)ItemPtr)+CurPos+Offset);
//						else if(n==2)
//							sscanf(match,"%i",((int *)ItemPtr)+CurPos+Offset);
//						else if(n==3)
//							sscanf(match,"%s",*(((char **)ItemPtr)+CurPos+Offset));
//						}
//					else
//						CurPos=NumItem;
//					}
//				}
//			}
//		else
//			{
//			for(j=0;j<NumItem;j++)
//				{ //added these braces
//					for(i=0;i<NumItem2;i++)
//						{
//						fscanf(dat,"%s",match);
//						if((match[0]!='[')&&(!feof(dat)))
//							{
//							FindFlag++;
//							if(n==1)
//								sscanf(match,"%lg",((double **)ItemPtr)[j+Offset]+i+Offset); //changed from [j+Offset]+i+Offset to +j+Offset+i, as ItemPtr isn't an array - 01/10: changed it back
// 							else
//								sscanf(match,"%i",((int **)ItemPtr)[j+Offset]+i+Offset);
//							}
//						else
//							{
//							i=NumItem2;
//							j=NumItem;
//							}
//						}
//					//Offset=Offset+(NumItem2-1); //added this line to get the correct offset in address position when incrementing j
//				} //added these braces
//			}
//		}
////	fprintf(stderr,"%s\n",SItemName);
//	return FindFlag;
}

#define PREVSCALE 50
#define IMM_PREVSCALE 5
#define BWCOLS 58
#define PREVCOLS 87

void CaptureBitmap(int ns,int tp)
{
	int i,j,x,y,f,mi;
	static double mx;
	static int fst=1;
	double prev,age;

	mi=(int) (P.bwidth*P.bheight);
	if(fst)
		{
		fst=0;
		mx=0;
		for(i=0;i<mi;i++) bmi[i]=0;
		for(i=0;i<P.N;i++)
			{
			x=((int) (Households[Hosts[i].hh].loc_x*P.scalex))-P.bminx;
			y=((int) (Households[Hosts[i].hh].loc_y*P.scaley))-P.bminy;
			if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
				{
				j=y*bmh->width+x;
				if((j<bmh->imagesize)&&(j>=0))
					{
					bmi[j]++;
					if(bmi[j]>mx) mx=(double) bmi[j];
					}
				}
			}
		mx=log(1.001*mx);
		for(i=0;i<P.NMC;i++)
			{
			f=0;
			if((i/P.nmch==(i+1)/P.nmch)&&((Mcells[i].country!=Mcells[i+1].country)||((P.DoAdunitBoundaryOutput)&&((AdUnits[Mcells[i].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor !=(AdUnits[Mcells[i+1].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor))))
			{
				f=1;
				//fprintf(stderr,"mcells %i, %i\n",i,i+1);
				//fprintf(stderr,"mcells countries %i, %i\n",Mcells[i].country,Mcells[i+1].country);
				//fprintf(stderr,"mcells adunits %i, %i\n",Mcells[i].adunit,Mcells[i+1].adunit);
				//fprintf(stderr,"mcells pop %i, %i\n",Mcells[i].n,Mcells[i+1].n);
			}
			else if((i>0)&&(i/P.nmch==(i-1)/P.nmch)&&(Mcells[i].country!=Mcells[i-1].country)) f=1;
			else if((i<P.NMC-P.nmch)&&((Mcells[i].country!=Mcells[i+P.nmch].country)||((P.DoAdunitBoundaryOutput)&&((AdUnits[Mcells[i].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor!=(AdUnits[Mcells[i+P.nmch].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor)))) f=1;
			else if((i>P.nmch)&&(Mcells[i].country!=Mcells[i-P.nmch].country)) f=1;
			if(f)
				{
				x=(int) (P.mcwidth*(((double) (i/P.nmch))+0.5)*P.scalex)-P.bminx;
				y=(int) (P.mcheight*(((double) (i%P.nmch))+0.5)*P.scaley)-P.bminy;
				if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
					{
					j=y*bmh->width+x;
					if((j<bmh->imagesize)&&(j>=0)) bmi[j]=-1;
					}
				}
			}
		for(i=0;i<P.bwidth/2;i++)
			{
			prev=floor(1.99999*((double) i)*PREVCOLS/((double) P.bwidth)*2);
			f=BWCOLS+((int) prev);
			for(j=0;j<10;j++)
				{
				bm[(j+P.bheight+5)*bmh->width+P.bwidth/4+i]=f;
				}
			}
		}
#pragma omp parallel for private(i,j,prev) schedule(static,5000)
	for(i=0;i<mi;i++)
		{
		if(bmi[i]==-1)
			bm[i]=BWCOLS-1; /* black for country boundary */
		else if((bmi4[i]==0)&&(bmi2[i]==0)&&(bmi3[i]>0))
			bm[i]=(unsigned char) (3*BWCOLS+BWCOLS*log(bmi3[i])/mx);  /* green for recovered */
		else if(bmi2[i]>0)
			bm[i]=(unsigned char) (BWCOLS+BWCOLS*log(bmi2[i])/mx); /* red for infected */
		else if(bmi4[i]>0)
			bm[i]=(unsigned char) (2*BWCOLS+BWCOLS*log(bmi[i])/mx); /* blue for treated */
		else if(bmi[i]>0)
			bm[i]=(unsigned char) (BWCOLS*log(bmi[i])/mx); /* grey for just people */
		else
			bm[i]=0;
		}
}


//void CaptureBitmap(int ns,int tp)
//{
//	int i,j,x,y,f,mi;
//	static double mx;
//	static int fst=1;
//	double prev,age;
//
//	mi=(int) (P.bwidth*P.bheight);
//	if(fst)
//		{
//		fst=0;
//		mx=0;
//		for(i=0;i<mi;i++) bmi[i]=0;
//		for(i=0;i<P.N;i++)
//			{
//			x=((int) (Households[Hosts[i].hh].loc_x*P.scalex))-P.bminx;
//			y=((int) (Households[Hosts[i].hh].loc_y*P.scaley))-P.bminy;
//			if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
//				{
//				j=y*bmh->width+x;
//				if((j<bmh->imagesize)&&(j>=0))
//					{
//					bmi[j]++;
//					if(bmi[j]>mx) mx=(double) bmi[j];
//					}
//				}
//			}
//		mx=log(1.001*mx);
//		for(i=0;i<P.NMC;i++)
//			{
//			f=0;
//			if((i/P.nmch==(i+1)/P.nmch)&&((Mcells[i].country!=Mcells[i+1].country)||((P.DoAdunitBoundaryOutput)&&((AdUnits[Mcells[i].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor !=(AdUnits[Mcells[i+1].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor)))) f=1;
//			else if((i>0)&&(i/P.nmch==(i-1)/P.nmch)&&(Mcells[i].country!=Mcells[i-1].country)) f=1;
//			else if((i<P.NMC-P.nmch)&&((Mcells[i].country!=Mcells[i+P.nmch].country)||((P.DoAdunitBoundaryOutput)&&((AdUnits[Mcells[i].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor!=(AdUnits[Mcells[i+P.nmch].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor)))) f=1;
//			else if((i>P.nmch)&&(Mcells[i].country!=Mcells[i-P.nmch].country)) f=1;
//			if(f)
//				{
//				x=(int) (P.mcwidth*(((double) (i/P.nmch))+0.5)*P.scalex)-P.bminx;
//				y=(int) (P.mcheight*(((double) (i%P.nmch))+0.5)*P.scaley)-P.bminy;
//				if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
//					{
//					j=y*bmh->width+x;
//					if((j<bmh->imagesize)&&(j>=0)) bmi[j]=-1;
//					}
//				}
//			}
//		for(i=0;i<P.bwidth/2;i++)
//			{
//			prev=floor(1.99999*((double) i)*PREVCOLS/((double) P.bwidth)*2);
//			f=BWCOLS+((int) prev);
//			for(j=0;j<10;j++)
//				{
//				bm[(j+P.bheight+5)*bmh->width+P.bwidth/4+i]=f;
//				}
//			}
//		}
//#pragma omp parallel for private(i,j,prev) schedule(static,5000)
//	for(i=0;i<mi;i++)
//		{
//		if(bmi[i]==-1)
//			bm[i]=BWCOLS-1; /* black for country boundary */
//		else if((bmi4[i]==0)&&(bmi2[i]==0)&&(bmi3[i]>0))
//			bm[i]=(unsigned char) (3*BWCOLS+BWCOLS*log(bmi3[i])/mx);  /* green for recovered */
//		else if(bmi2[i]>0)
//			bm[i]=(unsigned char) (BWCOLS+BWCOLS*log(bmi2[i])/mx); /* red for infected */
//		else if(bmi4[i]>0)
//			bm[i]=(unsigned char) (2*BWCOLS+BWCOLS*log(bmi[i])/mx); /* blue for treated */
//		else if(bmi[i]>0)
//			bm[i]=(unsigned char) (BWCOLS*log(bmi[i])/mx); /* grey for just people */
//		}
//}

void CaptureMeanBitmap(int ns)
	{
	int i,j,mi;
	static double mx;
	static int fst=1;

	mi=(int) bmh->imagesize;
	if(fst)
		{
		fst=0;
		mx=0;
		for(j=0;j<bmh->imagesize;j++)
			if(bmi[j]>mx) mx=bmi[j];
		mx=log(1.001*mx);
		}
#pragma omp parallel for private(i) schedule(static,500) //added i to private
	for(i=0;i<mi;i++)
		{
		bmi2[i]=ceil(TSMean[ns].bmi2[i]/((float) P.NRactual));
		bmi3[i]=ceil(TSMean[ns].bmi3[i]/((float) P.NRactual));
		bmi4[i]=ceil(TSMean[ns].bmi4[i]/((float) P.NRactual));
		}
#pragma omp parallel for private(i) schedule(static,500) //added i to private
	for(i=0;i<mi;i++)
		{
		if(bmi[i]==-1)
			bm[i]=BWCOLS-1; /* black for country boundary */
		else if((bmi4[i]==0)&&(bmi2[i]==0)&&(bmi3[i]>0))
			bm[i]=(unsigned char) (3*BWCOLS+BWCOLS*log(bmi3[i])/mx);  /* green for recovered */
		else if(bmi2[i]>0)
			bm[i]=(unsigned char) (BWCOLS+BWCOLS*log(bmi2[i])/mx); /* red for infected */
		else if(bmi4[i]>0)
			bm[i]=(unsigned char) (2*BWCOLS+BWCOLS*log(bmi4[i])/mx); /* blue for treated */
		else if(bmi[i]>0)
			bm[i]=(unsigned char) (BWCOLS*log(bmi[i])/mx); /* grey */
		}
}


void OutputBitmap(double t2,int tp)
{
	FILE *dat;
	char buf[1024],OutF[1024];
	int i,j;
	static int cn1=0,cn2=0,cn3=0,cn4=0;
	size_t a;

	if(tp==0)
		{
		j=cn1;
		cn1++;
		sprintf(OutF,"%s",OutFile);
		}
	else if(tp==1)
		{
		j=cn2;
		cn2++;
		sprintf(OutF,"Mean.%s",OutFile);
		}
	else if(tp==2)
		{
		j=cn3;
		cn3++;
		sprintf(OutF,"Min.%s",OutFile);
		}
	else if(tp==3)
		{
		j=cn4;
		cn4++;
		sprintf(OutF,"Max.%s",OutFile);
		}

#ifdef IMAGE_MAGICK
	using namespace Magick;
	fprintf(stderr,"\noutputing ImageMagick stuff");
	sprintf(buf,"%s.bmp",OutF);
	if(!(dat=fopen(buf,"wb"))) ERR_CRITICAL("Unable to open bitmap file\n");
	fprintf(dat,"BM");
	//fwrite_big((void *) &bmf,sizeof(unsigned char),(sizeof(bitmap_header)/sizeof(unsigned char))+bmh->imagesize,dat);
	fwrite_big((void *) bmf,sizeof(bitmap_header),1,dat);
	for(i=0;i<bmh->imagesize;i++) fputc(bm[i],dat);
	fclose(dat);
	Image bmap(buf);
	sprintf(buf,"%s.%d.png",OutF,j);
	ColorRGB white(1.0,1.0,1.0);
	bmap.transparent(white);
	bmap.write(buf);
#elif defined(WIN32_BM)	
//Windows specific bitmap manipulation code - could be recoded using LIBGD or another unix graphics library
	using namespace Gdiplus;  

	HBITMAP hbm;
	wchar_t wbuf[1024];
	static UINT palsize;
	static ColorPalette* palette;

//Add new bitmap to AVI
	if((P.OutputBitmap==1)&&(tp==0)) AddAviFrame(avi,bmpdib,(unsigned char *) (&bmh->palette[0][0]));

//This transfers HBITMAP to GDI+ Bitmap object
	Bitmap* gdip_bmp = Bitmap::FromHBITMAP(bmpdib,NULL);
//Now change White in palette (first entry) to be transparent
	if((cn1==1)&&(tp==0))
		{
		palsize = gdip_bmp->GetPaletteSize();
		palette = (ColorPalette*)malloc(palsize);
		}
	i=gdip_bmp->GetPalette(palette,palsize);
	palette->Flags=PaletteFlagsHasAlpha;
	palette->Entries[0] = 0x00ffffff; // Transparent white 
	gdip_bmp->SetPalette(palette);
//Now save as png
	sprintf(buf,"%s.%05i.png",OutF,j+1); //sprintf(buf,"%s.ge\\%s.%05i.png",OutFileBase,OutF,j+1);
	mbstowcs_s(&a, wbuf, strlen(buf)+1, buf, _TRUNCATE);
	gdip_bmp->Save(wbuf, &encoderClsid,NULL);
	delete gdip_bmp;
#else
	sprintf(buf,"%s.bmp",OutF);
	if(!(dat=fopen(buf,"wb"))) ERR_CRITICAL("Unable to open bitmap file\n");
	fprintf(dat,"BM");
	fwrite_big((void *) bmf,sizeof(unsigned char),sizeof(bitmap_header)/sizeof(unsigned char)+bmh->imagesize,dat);
	fclose(dat);
#endif
}

void InitBMHead()
{
	int i,j,k,k2,R,B,G;
	double x,f;

	fprintf(stderr,"Initialising bitmap\n");
	k=4*P.bwidth*P.bheight2;
	k2=sizeof(bitmap_header)/sizeof(unsigned char);

	if(!(bmf=(unsigned char *) malloc((k+k2)*sizeof(unsigned char))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	bm=&(bmf[k2]);
	bmp=&(bmf[12]);
	bmh=(bitmap_header *) bmf;
	bmh->spare=0;
	bmh->boffset=2+sizeof(bitmap_header);
	bmh->headersize=40;
	bmh->width=P.bwidth;
	bmh->height=P.bheight2;
	bmh->PlanesAndBitspp=8*65536+1;
	bmh->compr=0;
	bmh->imagesize=bmh->width*bmh->height;
	bmh->filesize=2+bmh->imagesize+((unsigned int) sizeof(bitmap_header));
	bmh->hres=bmh->vres=(int) (bmh->width*10);
	bmh->colours=BWCOLS*4;
	bmh->impcol=0;
	for(i=0;i<256;i++)
		bmh->palette[i][3]=0;
	for(j=0;j<BWCOLS;j++)
		{
		B=255-255*j/(BWCOLS-1);
		bmh->palette[j][0]=bmh->palette[j][1]=bmh->palette[j][2]=(unsigned char) B;
		bmh->palette[BWCOLS+j][0]=0;
		bmh->palette[BWCOLS+j][1]=0;
		bmh->palette[BWCOLS+j][2]=(unsigned char) B;
		bmh->palette[2*BWCOLS+j][0]=(unsigned char) B;
		bmh->palette[2*BWCOLS+j][1]=0;
		bmh->palette[2*BWCOLS+j][2]=0;
		bmh->palette[3*BWCOLS+j][0]=0;
		bmh->palette[3*BWCOLS+j][1]=(unsigned char) B;
		bmh->palette[3*BWCOLS+j][2]=0;
		}
	if(!(bmi=(float *) malloc(bmh->imagesize*sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	if(!(bmi2=(float *) malloc(bmh->imagesize*sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	if(!(bmi3=(float *) malloc(bmh->imagesize*sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	if(!(bmi4=(float *) malloc(bmh->imagesize*sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");

#ifdef WIN32_BM
	bmpdib=CreateDIBSection(GetDC(NULL),(BITMAPINFO *)bmp,DIB_RGB_COLORS,(void **) &bm,NULL,NULL);
	Gdiplus::GdiplusStartupInput gdiplusStartupInput;
	Gdiplus::GdiplusStartup(&m_gdiplusToken, &gdiplusStartupInput, NULL);
	
	using namespace Gdiplus;
	UINT  num = 0;          // number of image encoders
	UINT  size = 0;         // size of the image encoder array in bytes

	ImageCodecInfo* pImageCodecInfo = NULL;
	GetImageEncodersSize(&num, &size);
	pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
	GetImageEncoders(num, size, pImageCodecInfo);
	for(UINT j = 0; j < num; ++j)
		{
		if(wcscmp(pImageCodecInfo[j].MimeType, L"image/png") == 0 )
			{
			encoderClsid = pImageCodecInfo[j].Clsid;
			j=num;
			}    
		}
	free(pImageCodecInfo);
	char buf[1024];
	sprintf(buf,"%s.ge",OutFileBase);
	if(!(CreateDirectory(buf,NULL))) fprintf(stderr,"Unable to create directory %s\n",buf);

#endif


}


void HSB2RGB(double h1, double s1, double v1, int *r, int *g, int *b)
{
int i;
double  h,s,v,p, q, t;
double f;

  h=h1;s=s1;v=v1;
  if(s==0)
	{
    *r = *g = *b = v*255;  /* If sat=0 the color is grey as determined by value */
    return;
	}
  else
	{
    h *= 6;
    if(h>=6)  h=0;
    i = (int) h;
    f = h-floor(h);
    p = 255*v*(1-s);
    q = 255*v*(1-s*f);
    t = 255*v*(1-(s*(1-f)));
    switch(i)
		{
		case 0 : *r = v*255+0.5;     /* Add 0.5 to all values in order to round */
               *g = t+0.5;         /* instead of truncating */
               *b = p+0.5;
               return;
		case 1 : *r = q+0.5;
               *g = v*255+0.5;
               *b = p+0.5;
               return;
		case 2 : *r = p+0.5;
               *g = v*255+0.5;
               *b = t+0.5;
               return;
		case 3 : *r = p+0.5;
               *g = q+0.5;
               *b = v*255+0.5;
               return;
		case 4 : *r = t+0.5;
               *g = p+0.5;
               *b = v*255+0.5;
               return;
		case 5 : *r = v*255+0.5;
               *g = p+0.5;
               *b = q+0.5;
               return;
		}
    }
    
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
