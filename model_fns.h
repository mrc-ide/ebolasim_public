void ReadParams(char*, char*);
void ReadInterventions(char*);
int GetXMLNode(FILE*, char*, char*, char*, int);
void ReadAirTravel(char*);
void SetupModel(char*, char*, char*, char*);
void SetupPopulation(char*, char*, char*);
void SetupAirports(void);
void SetupRoads(void); //added new function to take care of roads: ggilani - 12/02/15
void InitKernel(int, double);
void AssignHouseholdAges(int, int, int);
void AssignPeopleToPlaces(void);
void StratifyPlaces(void);
void LoadPeopleToPlaces(char*);
void SavePeopleToPlaces(char*);
void InitModel(int); //adding run number as a parameter for event log: ggilani - 15/10/2014
void SeedInfection(double, int*, int, int); //adding run number as a parameter for event log: ggilani - 15/10/2014
int CalcSeedResist(void);
void RunModel(int); //adding run number as a parameter for event log: ggilani - 15/10/2014
void TravelReturnSweep(double);
void TravelDepartSweep(double);
double CalcHouseInf(int, unsigned short int);
double CalcPlaceInf(int, int, unsigned short int);
double CalcSpatialInf(int, unsigned short int);
double CalcPersonInf(int, unsigned short int);
double CalcHouseSusc(int, unsigned short int, int, int);
double CalcPlaceSusc(int, int, unsigned short int, int, int);
double CalcSpatialSusc(int, unsigned short int, int, int);
double CalcPersonSusc(int, unsigned short int, int, int);
void InfectSweep(double, int); //added int as argument to InfectSweep to record run number: ggilani - 15/10/14
void IncubRecoverySweep(double, int); //added int as argument to record run number: ggilani - 15/10/14
int TreatSweep(double);
//void HospitalSweep(double); //added hospital sweep function: ggilani - 10/11/14
void HospitalSweepAdunits(double); //added hospitalisation by adunit sweep function: ggilani - 24/11/14
void ContactTracingSweep(double); // added function to update contact tracing number
void VaccSweep(double); //added function to process ring vaccination queue: ggilani - 21/08/19
void UpdateHospitals(double); //added function to update hospital parameters on each time step: ggilani - 11/03/2017
void UpdateContactTracing(double); //added function to update contact tracing capacity
void UpdateSDB(double); //added function to update safe burial capacity
void UpdateVaccination(double, int); //added function to update vaccination parameters at each time step: ggilani - 29/05/2019
void UpdateCaseDetection(double); //added function to update vaccination parameters at each time step: ggilani - 29/05/2019
void SaveAgeDistrib(void);
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
void RecordSample(double, int);
//adding function to record an event: ggilani - 10/10/2014
void RecordEvent(double, int, int, int, int); //added int as argument to InfectSweep to record run number: ggilani - 15/10/14
void DoInfect(int, double, int, int); //added int as argument to InfectSweep to record run number: ggilani - 15/10/14
void DoImmune(int);
void DoIncub(int, unsigned short int, int, int); //added int as argument to record run number: ggilani - 23/10/14
void DoDetectedCase(int, double, unsigned short int, int);
void DoCase(int, double, unsigned short int, int);
void DoFalseCase(int, double, unsigned short int, int);
void DoRecover(int, int, int); //added int as argument to record run number: ggilani - 23/10/14
void DoDeath(int, int, int); //added int as argument to record run number: ggilani - 23/10/14
void DoPlaceClose(int, int, unsigned short int, int, int);
void DoTreatCase(int, unsigned short int, int);
void DoProph(int, unsigned short int, int);
void DoPrivateTreatCase(int, unsigned short int, int);
void DoPrivateProph(int, unsigned short int, int);
void DoProphNoDelay(int, unsigned short int, int, int);
int DoVacc(int, int, int);
void DoVaccNoDelay(int, int);
double CalcPrevalenceDepTransmission(int);
void CalcOriginDestMatrix_adunit(void); //added function to calculate origin destination matrix: ggilani 28/01/15
int BedsAvailablePerAdUnit(double, int); //added function to check current number of beds available in an admin unit
void DetermineCellsWithCapitalCities(void); //added function to quickly check capital cities
double ExpKernel(double);
double PowerKernel(double);
double PowerKernelB(double);
double PowerKernelUS(double);
double PowerExpKernel(double);
double GaussianKernel(double);
double StepKernel(double);
double numKernel(double);
double dist2UTM(double, double, double, double);
double dist2(person*, person*);
double dist2_cc(cell*, cell*);
double dist2_cc_min(cell*, cell*);
double dist2_mm(microcell*, microcell*);
double dist2_mm_min(microcell*, microcell*); //added min distance from microcell to microcell for origin-destination matrix: ggilani 28/01/15
double dist2_raw(double, double, double, double);
int GetInputParameter(FILE*, FILE*, char*, char*, void*, int, int, int);
int GetInputParameter2(FILE*, FILE*, char*, char*, void*, int, int, int);
int GetInputParameter3(FILE*, const char*, const char*, void*, int, int, int);
//int GetInputParameter3(FILE *,char *,char *,void *, int,int,int);
void CaptureBitmap(int, int);
void CaptureMeanBitmap(int);
void OutputBitmap(double, int);
void InitBMHead();
void HSB2RGB(double, double, double, int*, int*, int*);
void HandleBreak(int);
/* RANDLIB functions */
long ignbin(long, double);
long ignpoi(double);
long ignbin_mt(long, double, int);
long ignpoi_mt(double, int);
double ranf(void);
double ranf_mt(int);
void setall(long, long);
double sexpo_mt(int);
double sexpo(void);
long mltmod(long, long, long);
double snorm(void);
double snorm_mt(int);
double fsign(double, double);
double gengam(double, double);
double gengam_mt(double, double, int);
//added some new beta, gamma generating functions: ggilani 27/11/14
double gen_norm(double, double);
double gen_norm_mt(double, double, int);
double gen_gamma(double, double);
double gen_gamma_mt(double, double, int);
double gen_beta(double, double);
double gen_beta_mt(double, double, int);
//added some new lognormal sampling functions: ggilani 09/02/17
double gen_lognormal(double, double);
double gen_lognormal_mt(double, double, int);
double sgamma(double);
double sgamma_mt(double, int);
void SampleWithoutReplacement(int, int, int);
double gammln(double);
double WeibullCDF(int, double, double);


#pragma once
