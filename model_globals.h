param P;
person* Hosts;
household* Households;
popvar State, StateT[MAX_NUM_THREADS];
cell* Cells, ** CellLookup;
int* RevCellLookup;
microcell* Mcells, ** McellLookup;
place** Places;
adminunit AdUnits[MAX_ADUNITS];
results* TimeSeries, * TSMean, * TSVar, * TSMeanNE, * TSVarNE, * TSMeanE, * TSVarE;
results* TimeSeries_G, * TSMean_G, * TSVar_G, * TSMeanNE_G, * TSVarNE_G, * TSMeanE_G, * TSVarE_G;
results* TimeSeries_L, * TSMean_L, * TSVar_L, * TSMeanNE_L, * TSVarNE_L, * TSMeanE_L, * TSVarE_L;

//added declaration of pointer to events log: ggilani - 10/10/2014
events* InfEventLog;
int* nEvents;

//added counter for number of types ranf and ranf_mt are called
int count_ranf, count_ranfmt;

airport* Airports;
double (*Kernel)(double);
double* nKernel, * nKernelHR, ** PopDensity, * mcell_dens;
int* mcell_adunits, * mcell_num, * mcell_country;
double inftype[INFECT_TYPE_MASK], inftype_av[INFECT_TYPE_MASK], infcountry[MAX_COUNTRIES], infcountry_av[MAX_COUNTRIES], infcountry_num[MAX_COUNTRIES];
double indivR0[MAX_SEC_REC][MAX_GEN_REC], indivR0_av[MAX_SEC_REC][MAX_GEN_REC];
double inf_household[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], denom_household[MAX_HOUSEHOLD_SIZE + 1];
double inf_household_av[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], AgeDist[NUM_AGE_GROUPS], AgeDist2[NUM_AGE_GROUPS];
double vaccdose_dist[NUM_VACCDOSE_GROUPS], vaccdosering_dist[NUM_VACCDOSE_GROUPS], vaccdosecell_dist[NUM_VACCDOSECELL_GROUPS], vaccdistance_dist[NUM_VACCDIST_GROUPS], vaccpop_dist[NUM_POP_GROUPS];
double case_household[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], case_household_av[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1];
double PropPlaces[NUM_AGE_GROUPS * AGE_GROUP_WIDTH][NUM_PLACE_TYPES];
double PropPlacesC[NUM_AGE_GROUPS * AGE_GROUP_WIDTH][NUM_PLACE_TYPES], AirTravelDist[MAX_DIST];
double sinx[361], cosx[361], asin2sqx[1001];
double PeakHeightSum, PeakHeightSS, PeakTimeSum, PeakTimeSS;
bitmap_header* bmh, * bmh2, * bmh3;
void* BinFileBuf;
bin_file* BF;
unsigned char* bmf, * bm, * bmp, * bmf2, * bm2, * bmp2, * bmf3, * bm3, * bmp3;
float* bmi, * bmi2, * bmi3, * bmi4, * bmi5, * bmi6, * bmi2min, * bmi2max, * bmi2mean, * bmi3min, * bmi3max, * bmi3mean, * bmi7, * bmi7min, * bmi7max, * bmi7mean;

#ifdef WIN32_BM
HAVI avi, avi2, avi3;
HBITMAP bmpdib, bmpdib2, bmpdib3;
ULONG_PTR m_gdiplusToken;
CLSID  encoderClsid;
#endif
char OutFile[1024], OutFileBase[1024], OutDensFile[1024], SnapshotLoadFile[1024], SnapshotSaveFile[1024], RadiationFile[1024], RoadNetworkFile[1024]; //added radiation file for radiation model: ggilani 09/02/15

int ns, DoInitUpdateProbs, InterruptRun = 0;
int netbuf[NUM_PLACE_TYPES_NOAIR * 1000000], ** SamplingQueue;
unsigned int cntr = 0;
int PlaceDistDistrib[NUM_PLACE_TYPES][MAX_DIST], PlaceSizeDistrib[NUM_PLACE_TYPES][MAX_PLACE_SIZE];


FILE* KMLFile, * KMLFile2;

/* RANDLIB static variables */
long* Xcg1, * Xcg2;

