/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "model_constants.h"
#include "model_structs.h"
#include "model_fns.h"
#include "model_globals_ext.h"
#include "binio.h"



void ReadParams(char* ParamFile, char* PreParamFile)
{
	FILE* dat, * dat2, * dat3;
	double s, t, CumAgeDist[NUM_AGE_GROUPS + 1], AgeSuscScale;
	int i, j, k, nc, na;
	char buf[1024], * CountryNames[MAX_COUNTRIES], CountryNameBuf[128 * MAX_COUNTRIES];
	char** AdunitNames, * AdunitNamesBuf;

	if (!(dat = fopen(ParamFile, "r"))) ERR_CRITICAL("Unable to open parameter file\n");
	dat2 = fopen(PreParamFile, "r");

	AgeSuscScale = 1.0;
	GetInputParameter(dat, dat2, "Update timestep", "%lf", (void*)&(P.TimeStep), 1, 1, 0);
	GetInputParameter(dat, dat2, "Sampling timestep", "%lf", (void*)&(P.SampleStep), 1, 1, 0);
	if (P.TimeStep > P.SampleStep) ERR_CRITICAL("Update step must be smaller than sampling step\n");
	t = ceil(P.SampleStep / P.TimeStep - 1e-6);
	P.UpdatesPerSample = (int)t;
	P.TimeStep = P.SampleStep / t;
	P.TimeStepsPerDay = ceil(1.0 / P.TimeStep - 1e-6);
	P.TimeStepsPerDayInt = (int)P.TimeStepsPerDay;
	P.TimeStepsPerYear = P.TimeStepsPerDay * DAYS_PER_YEAR;
	fprintf(stderr, "Update step = %lf\nSampling step = %lf\nUpdates per sample=%i\nTimeStepsPerDay=%lf\n", P.TimeStep, P.SampleStep, P.UpdatesPerSample, P.TimeStepsPerDay);
	GetInputParameter(dat, dat2, "Sampling time", "%lf", (void*)&(P.SampleTime), 1, 1, 0);
	P.NumSamples = 1 + (int)ceil(P.SampleTime / P.SampleStep);
	GetInputParameter(dat, dat2, "Population size", "%i", (void*)&(P.N), 1, 1, 0);
	GetInputParameter(dat, dat2, "Number of realisations", "%i", (void*)&(P.NR), 1, 1, 0);
	if (!GetInputParameter2(dat, dat2, "Number of non-extinct realisations", "%i", (void*)&(P.NRN), 1, 1, 0)) P.NRN = P.NR;
	if (!GetInputParameter2(dat, dat2, "Maximum number of cases defining small outbreak", "%i", (void*)&(P.SmallEpidemicCases), 1, 1, 0)) P.SmallEpidemicCases = -1;
	GetInputParameter(dat, dat2, "Number of spatial cells", "%i", (void*)&(P.NC), 1, 1, 0);
	GetInputParameter(dat, dat2, "Number of micro-cells per spatial cell width", "%i", (void*)&(P.NMCL), 1, 1, 0);
	//added parameter to reset seeds after every run
	if (!GetInputParameter2(dat, dat2, "Reset seeds for every run", "%i", (void*)&(P.ResetSeeds), 1, 1, 0)) P.ResetSeeds = 0;
	if (P.ResetSeeds)
	{
		if (!GetInputParameter2(dat, dat2, "Keep same seeds for every run", "%i", (void*)&(P.KeepSameSeeds), 1, 1, 0)) P.KeepSameSeeds = 0; //added this to control which seeds are used: ggilani 27/11/19
		if (!GetInputParameter2(dat, dat2, "Use fixed input seeds", "%i", (void*)&(P.DoFixedSeeds), 1, 1, 0)) P.DoFixedSeeds = 0; //added this to allow a fixed list of seeds for Janetta: ggilani 08/03/2023
	}
	if (!GetInputParameter2(dat, dat2, "Reset seeds after intervention", "%i", (void*)&(P.ResetSeedsPostIntervention), 1, 1, 0)) P.ResetSeedsPostIntervention = 0;
	if (P.ResetSeedsPostIntervention)
	{
		if (!GetInputParameter2(dat, dat2, "Time to reset seeds after intervention", "%i", (void*)&(P.TimeToResetSeeds), 1, 1, 0)) P.TimeToResetSeeds = 1e6;
	}
	if (!GetInputParameter2(dat, dat2, "Include households", "%i", (void*)&(P.DoHouseholds), 1, 1, 0)) P.DoHouseholds = 0;
	if (P.DoHouseholds)
	{
		GetInputParameter(dat, dat2, "Household size distribution", "%lf", (void*)P.HouseholdSizeDistrib[0], MAX_HOUSEHOLD_SIZE, 1, 0);
		GetInputParameter(dat, dat2, "Household attack rate", "%lf", (void*)&(P.HouseholdTrans), 1, 1, 0);
		GetInputParameter(dat, dat2, "Household transmission denominator power", "%lf", (void*)&(P.HouseholdTransPow), 1, 1, 0);
		if (!GetInputParameter2(dat, dat2, "Median household income", "%lf", (void*)&(P.MedianIncome), 1, 1, 0)) P.MedianIncome = 10000;
		if (!GetInputParameter2(dat, dat2, "Household income weibull power", "%lf", (void*)&(P.IncomeWeibullPower), 1, 1, 0)) P.IncomeWeibullPower = 1.5;
		if (!GetInputParameter2(dat, dat2, "Output household file", "%i", (void*)&(P.DoHouseholdOutput), 1, 1, 0)) P.DoHouseholdOutput = 0;
	}
	else
	{
		P.HouseholdTrans = 0.0;
		P.HouseholdTransPow = 1.0;
		P.HouseholdSizeDistrib[0][0] = 1.0;
		for (i = 1; i < MAX_HOUSEHOLD_SIZE; i++)
			P.HouseholdSizeDistrib[0][i] = 0;
		P.MedianIncome = 10000;
		P.IncomeWeibullPower = 1.5;
	}
	for (i = 1; i < MAX_HOUSEHOLD_SIZE; i++)
		P.HouseholdSizeDistrib[0][i] = P.HouseholdSizeDistrib[0][i] + P.HouseholdSizeDistrib[0][i - 1];
	for (i = 0; i < MAX_HOUSEHOLD_SIZE; i++)
		P.HouseholdDenomLookup[i] = 1 / pow(((double)(i + 1)), P.HouseholdTransPow);
	if (!GetInputParameter2(dat, dat2, "Include administrative units within countries", "%i", (void*)&(P.DoAdUnits), 1, 1, 0)) P.DoAdUnits = 0;
	if (!GetInputParameter2(dat, dat2, "Divisor for countries", "%i", (void*)&(P.CountryDivisor), 1, 1, 0)) P.CountryDivisor = 1;
	if (!GetInputParameter2(dat, dat2, "Output country file", "%i", (void*)&(P.DoCountryOutput), 1, 1, 0)) P.DoCountryOutput = 0;
	if (P.DoAdUnits)
	{
		if (!(AdunitNames = (char**)malloc(3 * ADUNIT_LOOKUP_SIZE * sizeof(char*)))) ERR_CRITICAL("Unable to allocate temp storage\n");
		if (!(AdunitNamesBuf = (char*)malloc(3 * ADUNIT_LOOKUP_SIZE * 360 * sizeof(char)))) ERR_CRITICAL("Unable to allocate temp storage\n");

		for (i = 0; i < ADUNIT_LOOKUP_SIZE; i++)
		{
			P.AdunitLevel1Lookup[i] = -1;
			AdunitNames[3 * i] = AdunitNamesBuf + 3 * i * 360;
			AdunitNames[3 * i + 1] = AdunitNamesBuf + 3 * i * 360 + 60;
			AdunitNames[3 * i + 2] = AdunitNamesBuf + 3 * i * 360 + 160;
		}
		for (i = 0; i < MAX_COUNTRIES; i++) { CountryNames[i] = CountryNameBuf + 128 * i; CountryNames[i][0] = 0; }
		na = (GetInputParameter2(dat, dat2, "Codes and country/province names for admin units", "%s", (void*)AdunitNames, 3 * ADUNIT_LOOKUP_SIZE, 1, 0)) / 3;
		if ((na > 0) && (GetInputParameter2(dat, dat2, "Number of countries to include", "%i", (void*)&nc, 1, 1, 0)))
		{
			P.DoAdunitBoundaries = (nc > 0);
			nc = abs(nc);
			P.AdunitLevel1Divisor = 1;
			P.AdunitLevel1Mask = 1000000000;
			GetInputParameter(dat, dat2, "List of names of countries to include", "%s", (nc > 1) ? ((void*)CountryNames) : ((void*)CountryNames[0]), nc, 1, 0);
			P.NumAdunits = 0;
			for (i = 0; i < na; i++)
				for (j = 0; j < nc; j++)
				{
					if ((AdunitNames[3 * i + 1][0]) && (!strcmp(AdunitNames[3 * i + 1], CountryNames[j])) && (atoi(AdunitNames[3 * i]) != 0))
					{
						AdUnits[P.NumAdunits].id = atoi(AdunitNames[3 * i]);
						P.AdunitLevel1Lookup[AdUnits[P.NumAdunits].id] = P.NumAdunits;
						if (strlen(AdunitNames[3 * i + 1]) < 100) strcpy(AdUnits[P.NumAdunits].cnt_name, AdunitNames[3 * i + 1]);
						if (strlen(AdunitNames[3 * i + 2]) < 200) strcpy(AdUnits[P.NumAdunits].ad_name, AdunitNames[3 * i + 2]);
						//						fprintf(stderr,"%i %s %s ## ",AdUnits[P.NumAdunits].id,AdUnits[P.NumAdunits].cnt_name,AdUnits[P.NumAdunits].ad_name);
						P.NumAdunits++;
					}
				}
		}
		else
		{
			if (!GetInputParameter2(dat, dat2, "Number of level 1 administrative units to include", "%i", (void*)&(P.NumAdunits), 1, 1, 0)) P.NumAdunits = 0;
			if (!GetInputParameter2(dat, dat2, "Divisor for level 1 administrative units", "%i", (void*)&(P.AdunitLevel1Divisor), 1, 1, 0)) P.AdunitLevel1Divisor = 1;
			if (!GetInputParameter2(dat, dat2, "Mask for level 1 administrative units", "%i", (void*)&(P.AdunitLevel1Mask), 1, 1, 0)) P.AdunitLevel1Mask = 1000000000;
			if (P.NumAdunits > 0)
			{
				P.DoAdunitBoundaries = 1;
				int AdunitList[MAX_ADUNITS];
				if (P.DoAdunitBoundaries > MAX_ADUNITS) ERR_CRITICAL("MAX_ADUNITS too small.\n");
				GetInputParameter(dat, dat2, "List of level 1 administrative units to include", "%i", (void*)AdunitList, P.NumAdunits, 1, 0);
				for (i = 0; i < P.NumAdunits; i++)
				{
					AdUnits[i].id = AdunitList[i];
					P.AdunitLevel1Lookup[(AdUnits[i].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor] = i;
					for (j = 0; j < na; j++)
						if (atoi(AdunitNames[3 * j]) == AdunitList[i])
						{
							if (strlen(AdunitNames[3 * j + 1]) < 100) strcpy(AdUnits[i].cnt_name, AdunitNames[3 * j + 1]);
							if (strlen(AdunitNames[3 * j + 2]) < 200) strcpy(AdUnits[i].ad_name, AdunitNames[3 * j + 2]);
							j = na;
						}
				}
				if (!GetInputParameter2(dat, dat2, "Median income of level 1 administrative units", "%lf", (void*)&(P.MedianIncomeByAdunit[0]), P.DoAdunitBoundaries, 1, 0))
				{
					for (i = 0; i < P.NumAdunits; i++)
						P.MedianIncomeByAdunit[i] = P.MedianIncome;
				}
			}
			else
				P.DoAdunitBoundaries = 0;
			free(AdunitNames);
			free(AdunitNamesBuf);
		}
		if (!GetInputParameter2(dat, dat2, "Output incidence by administrative unit", "%i", (void*)&(P.DoAdunitOutput), 1, 1, 0)) P.DoAdunitOutput = 0;
		if (!GetInputParameter2(dat, dat2, "Draw administrative unit boundaries on maps", "%i", (void*)&(P.DoAdunitBoundaryOutput), 1, 1, 0)) P.DoAdunitBoundaryOutput = 0;
		if (!GetInputParameter2(dat, dat2, "Correct administrative unit populations", "%i", (void*)&(P.DoCorrectAdunitPop), 1, 1, 0)) P.DoCorrectAdunitPop = 0;
		if (!GetInputParameter2(dat, dat2, "Fix population size at specified value", "%i", (void*)&(P.DoSpecifyPop), 1, 1, 0)) P.DoSpecifyPop = 0;
		fprintf(stderr, "Using %i administrative units\n", P.NumAdunits);
		if (!GetInputParameter2(dat, dat2, "Divisor for administrative unit codes for boundary plotting on bitmaps", "%i", (void*)&(P.AdunitBitmapDivisor), 1, 1, 0)) P.AdunitBitmapDivisor = 1;
		if (!GetInputParameter2(dat, dat2, "Only output household to place distance distribution for one administrative unit", "%i", (void*)&(P.DoOutputPlaceDistForOneAdunit), 1, 1, 0)) P.DoOutputPlaceDistForOneAdunit = 0;
		if (P.DoOutputPlaceDistForOneAdunit)
		{
			if (!GetInputParameter2(dat, dat2, "Administrative unit for which household to place distance distribution to be output", "%i", (void*)&(P.OutputPlaceDistAdunit), 1, 1, 0)) P.DoOutputPlaceDistForOneAdunit = 0;
		}
	}
	else
	{
		P.DoAdunitBoundaries = P.DoAdunitBoundaryOutput = P.DoAdunitOutput = P.DoCorrectAdunitPop = P.DoSpecifyPop = 0; P.AdunitLevel1Divisor = 1; P.AdunitLevel1Mask = 1000000000; P.AdunitBitmapDivisor = P.AdunitLevel1Divisor;
	}
	//added flag for outputting summary results
	if (!GetInputParameter2(dat, dat2, "Output summary results", "%i", (void*)&(P.DoSummaryOutput), 1, 1, 0)) P.DoSummaryOutput = 0;

	if (!GetInputParameter2(dat, dat2, "Include age", "%i", (void*)&(P.DoAge), 1, 1, 0)) P.DoAge = 1;
	if (!P.DoAge)
	{
		P.DoAge = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.PropAgeGroup[0][i] = 1.0 / NUM_AGE_GROUPS;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
		{
			P.InitialImmunity[i] = 0;
			P.AgeInfectiousness[i] = P.AgeSusceptibility[i] = 1;
			P.RelativeSpatialContact[i] = P.RelativeTravelRate[i] = 1.0;
		}
	}
	else
	{

		if (!GetInputParameter2(dat, dat2, "Output age file", "%i", (void*)&(P.DoAgeOutput), 1, 1, 0)) P.DoAgeOutput = 0;
		if (P.DoHouseholds)
		{
			if (!GetInputParameter2(dat, dat2, "Initial immunity applied to all household members", "%i", (void*)&(P.DoWholeHouseholdImmunity), 1, 1, 0)) P.DoWholeHouseholdImmunity = 0;
		}
		else
			P.DoWholeHouseholdImmunity = 0;
		if (!GetInputParameter2(dat, dat2, "Initial immunity profile by age", "%lf", (void*)P.InitialImmunity, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.InitialImmunity[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Relative spatial contact rates by age", "%lf", (void*)P.RelativeSpatialContact, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.RelativeSpatialContact[i] = 1;
		if (!GetInputParameter2(dat, dat2, "Age-dependent infectiousness", "%lf", (void*)P.AgeInfectiousness, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.AgeInfectiousness[i] = 1.0;
		if (!GetInputParameter2(dat, dat2, "Age-dependent susceptibility", "%lf", (void*)P.AgeSusceptibility, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.AgeSusceptibility[i] = 1.0;
		GetInputParameter(dat, dat2, "Age distribution of population", "%lf", (void*)P.PropAgeGroup[0], NUM_AGE_GROUPS, 1, 0);
		t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			t += P.PropAgeGroup[0][i];
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.PropAgeGroup[0][i] /= t;
		t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			if (P.AgeSusceptibility[i] > t) t = P.AgeSusceptibility[i];  //peak susc has to be 1
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.AgeSusceptibility[i] /= t;
		AgeSuscScale = t;
		if (P.DoHouseholds) P.HouseholdTrans *= AgeSuscScale;
		if (!GetInputParameter2(dat, dat2, "Relative travel rates by age", "%lf", (void*)P.RelativeTravelRate, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.RelativeTravelRate[i] = 1;
		if (!GetInputParameter2(dat, dat2, "WAIFW matrix", "%lf", (void*)P.WAIFW_Matrix, NUM_AGE_GROUPS, NUM_AGE_GROUPS, 0))
		{
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				for (j = 0; j < NUM_AGE_GROUPS; j++)
					P.WAIFW_Matrix[i][j] = 1.0;
		}
		else
		{
			/* WAIFW matrix needs to be scaled to have max value of 1.
			1st index of matrix specifies host being infected, second the infector.
			Overall age variation in infectiousness/contact rates/susceptibility should be factored
			out of WAIFW_matrix and put in Age dep infectiousness/susceptibility for efficiency. */
			t = 0;
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				for (j = 0; j < NUM_AGE_GROUPS; j++)
					if (P.WAIFW_Matrix[i][j] > t) t = P.WAIFW_Matrix[i][j];
			if (t > 0)
			{
				for (i = 0; i < NUM_AGE_GROUPS; i++)
					for (j = 0; j < NUM_AGE_GROUPS; j++)
						P.WAIFW_Matrix[i][j] /= t;
			}
			else
			{
				for (i = 0; i < NUM_AGE_GROUPS; i++)
					for (j = 0; j < NUM_AGE_GROUPS; j++)
						P.WAIFW_Matrix[i][j] = 1.0;
			}
		}
		P.DoDeath = 0;
		t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			t += P.AgeInfectiousness[i] * P.PropAgeGroup[0][i];
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.AgeInfectiousness[i] /= t;
	}
	if (!GetInputParameter2(dat, dat2, "Include spatial transmission", "%i", (void*)&(P.DoSpatial), 1, 1, 0)) P.DoSpatial = 1;
	GetInputParameter(dat, dat2, "Kernel type", "%i", (void*)&(P.MoveKernelType), 1, 1, 0);
	GetInputParameter(dat, dat2, "Kernel scale", "%lf", (void*)&(P.MoveKernelScale), 1, 1, 0);
	if (P.KernelOffsetScale != 1)
	{
		P.MoveKernelScale *= P.KernelOffsetScale;
	}
	if (!GetInputParameter2(dat, dat2, "Kernel 3rd param", "%lf", (void*)&(P.MoveKernelP3), 1, 1, 0)) P.MoveKernelP3 = 0;
	if (!GetInputParameter2(dat, dat2, "Kernel 4th param", "%lf", (void*)&(P.MoveKernelP4), 1, 1, 0)) P.MoveKernelP4 = 0;
	if (!GetInputParameter2(dat, dat2, "Kernel Shape", "%lf", (void*)&(P.MoveKernelShape), 1, 1, 0)) P.MoveKernelShape = 1.0;
	if (P.KernelPowerScale != 1)
	{
		P.MoveKernelShape *= P.KernelPowerScale;
	}
	if (!GetInputParameter2(dat, dat2, "Airport Kernel Type", "%i", (void*)&(P.AirportKernelType), 1, 1, 0)) P.AirportKernelType = P.MoveKernelType;
	if (!GetInputParameter2(dat, dat2, "Airport Kernel Scale", "%lf", (void*)&(P.AirportKernelScale), 1, 1, 0)) P.AirportKernelScale = P.MoveKernelScale;
	if (!GetInputParameter2(dat, dat2, "Airport Kernel Shape", "%lf", (void*)&(P.AirportKernelShape), 1, 1, 0)) P.AirportKernelShape = P.MoveKernelShape;
	if (!GetInputParameter2(dat, dat2, "Airport Kernel 3rd param", "%lf", (void*)&(P.AirportKernelP3), 1, 1, 0)) P.AirportKernelP3 = P.MoveKernelP3;
	if (!GetInputParameter2(dat, dat2, "Airport Kernel 4th param", "%lf", (void*)&(P.AirportKernelP4), 1, 1, 0)) P.AirportKernelP4 = P.MoveKernelP4;
	if (!GetInputParameter2(dat, dat2, "Include places", "%i", (void*)&(P.DoPlaces), 1, 1, 0)) P.DoPlaces = P.PlaceTypeNum = P.DoAirports = 0;
	if (P.DoPlaces)
	{
		GetInputParameter(dat, dat2, "Number of types of places", "%i", (void*)&(P.PlaceTypeNum), 1, 1, 0);
		if (P.PlaceTypeNum == 0) P.DoPlaces = P.DoAirports = 0;
	}
	if (P.DoPlaces)
	{
		if (P.PlaceTypeNum > NUM_PLACE_TYPES) ERR_CRITICAL("Too many place types\n");
		GetInputParameter(dat, dat2, "Minimum age for age group 1 in place types", "%i", (void*)P.PlaceTypeAgeMin, P.PlaceTypeNum, 1, 0);
		GetInputParameter(dat, dat2, "Maximum age for age group 1 in place types", "%i", (void*)P.PlaceTypeAgeMax, P.PlaceTypeNum, 1, 0);
		GetInputParameter(dat, dat2, "Proportion of age group 1 in place types", "%lf", (void*)&(P.PlaceTypePropAgeGroup), P.PlaceTypeNum, 1, 0);
		if (!GetInputParameter2(dat, dat2, "Proportion of age group 2 in place types", "%lf", (void*)&(P.PlaceTypePropAgeGroup2), P.PlaceTypeNum, 1, 0))
		{
			for (i = 0; i < NUM_PLACE_TYPES; i++)
			{
				P.PlaceTypePropAgeGroup2[i] = 0;
				P.PlaceTypeAgeMin2[i] = 0;
				P.PlaceTypeAgeMax2[i] = 1000;
			}
		}
		else
		{
			GetInputParameter(dat, dat2, "Minimum age for age group 2 in place types", "%i", (void*)P.PlaceTypeAgeMin2, P.PlaceTypeNum, 1, 0);
			GetInputParameter(dat, dat2, "Maximum age for age group 2 in place types", "%i", (void*)P.PlaceTypeAgeMax2, P.PlaceTypeNum, 1, 0);
		}
		if (!GetInputParameter2(dat, dat2, "Proportion of age group 3 in place types", "%lf", (void*)&(P.PlaceTypePropAgeGroup3), P.PlaceTypeNum, 1, 0))
		{
			for (i = 0; i < NUM_PLACE_TYPES; i++)
			{
				P.PlaceTypePropAgeGroup3[i] = 0;
				P.PlaceTypeAgeMin3[i] = 0;
				P.PlaceTypeAgeMax3[i] = 1000;
			}
		}
		else
		{
			GetInputParameter(dat, dat2, "Minimum age for age group 3 in place types", "%i", (void*)P.PlaceTypeAgeMin3, P.PlaceTypeNum, 1, 0);
			GetInputParameter(dat, dat2, "Maximum age for age group 3 in place types", "%i", (void*)P.PlaceTypeAgeMax3, P.PlaceTypeNum, 1, 0);
		}
		if (!GetInputParameter2(dat, dat2, "Kernel shape params for place types", "%lf", (void*)&(P.PlaceTypeKernelShape), P.PlaceTypeNum, 1, 0))
		{
			for (i = 0; i < NUM_PLACE_TYPES; i++)
			{
				P.PlaceTypeKernelShape[i] = P.MoveKernelShape;
				P.PlaceTypeKernelScale[i] = P.MoveKernelScale;
			}
		}
		else
			GetInputParameter(dat, dat2, "Kernel scale params for place types", "%lf", (void*)&(P.PlaceTypeKernelScale), P.PlaceTypeNum, 1, 0);
		if (!GetInputParameter2(dat, dat2, "Kernel 3rd param for place types", "%lf", (void*)&(P.PlaceTypeKernelP3), P.PlaceTypeNum, 1, 0))
		{
			for (i = 0; i < NUM_PLACE_TYPES; i++)
			{
				P.PlaceTypeKernelP3[i] = P.MoveKernelP3;
				P.PlaceTypeKernelP4[i] = P.MoveKernelP4;
			}
		}
		else
			GetInputParameter(dat, dat2, "Kernel 4th param for place types", "%lf", (void*)&(P.PlaceTypeKernelP4), P.PlaceTypeNum, 1, 0);
		if (!GetInputParameter2(dat, dat2, "Number of closest places people pick from (0=all) for place types", "%i", (void*)&(P.PlaceTypeNearestNeighb), P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeNearestNeighb[i] = 0;
		if (P.DoAdUnits)
		{
			if (!GetInputParameter2(dat, dat2, "Degree to which crossing administrative unit boundaries to go to places is inhibited", "%lf", (void*)&(P.InhibitInterAdunitPlaceAssignment), P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++)
					P.InhibitInterAdunitPlaceAssignment[i] = 0;
		}
		if (P.DoHouseholds)
		{
			if (!GetInputParameter2(dat, dat2, "Do place type members have household contacts", "%i", (void*)&(P.PlaceTypeDoHousehold), P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++)
					P.PlaceTypeDoHousehold[i] = 1;
		}
		//add a condition about household and place membership overlapping - ggilani 13/02/17
		if (P.DoHouseholds)
		{
			if (!GetInputParameter2(dat, dat2, "Do place and household membership overlap", "%i", (void*)&(P.PlaceHouseholdOverlap), 1, 1, 0))
			{
				P.PlaceHouseholdOverlap = 0;
			}
			if (P.PlaceHouseholdOverlap)
			{
				if (!GetInputParameter2(dat, dat2, "Divisor for households into places", "%lf", (void*)&(P.PlaceHouseholdDivisor), 1, 1, 0)) P.PlaceHouseholdDivisor = 2.0;
			}
		}
		if (!GetInputParameter2(dat, dat2, "Include air travel", "%i", (void*)&(P.DoAirports), 1, 1, 0)) P.DoAirports = 0;
#ifdef FAST_US
		P.DoAirports = 0;
#endif
		if (P.DoAirports)
		{
			if (!GetInputParameter2(dat, dat2, "Scaling factor for input file to convert to daily traffic", "%lf", (void*)&(P.AirportTrafficScale), 1, 1, 0)) P.AirportTrafficScale = 1.0;
			if (!GetInputParameter2(dat, dat2, "Proportion of hotel attendees who are local", "%lf", (void*)&(P.HotelPropLocal), 1, 1, 0)) P.HotelPropLocal = 0;
			if (!GetInputParameter2(dat, dat2, "Distribution of duration of air journeys", "%lf", (void*)&(P.JourneyDurationDistrib), MAX_TRAVEL_TIME, 1, 0))
			{
				P.JourneyDurationDistrib[0] = 1;
				for (i = 0; i < MAX_TRAVEL_TIME; i++)
					P.JourneyDurationDistrib[i] = 0;
			}
			if (!GetInputParameter2(dat, dat2, "Distribution of duration of local journeys", "%lf", (void*)&(P.LocalJourneyDurationDistrib), MAX_TRAVEL_TIME, 1, 0))
			{
				P.LocalJourneyDurationDistrib[0] = 1;
				for (i = 0; i < MAX_TRAVEL_TIME; i++)
					P.LocalJourneyDurationDistrib[i] = 0;
			}
			P.MeanJourneyTime = P.MeanLocalJourneyTime = 0;
			for (i = 0; i < MAX_TRAVEL_TIME; i++)
			{
				P.MeanJourneyTime += ((double)(i)) * P.JourneyDurationDistrib[i];
				P.MeanLocalJourneyTime += ((double)(i)) * P.LocalJourneyDurationDistrib[i];
			}
			fprintf(stderr, "Mean duration of local journeys = %lf days\n", P.MeanLocalJourneyTime);
			for (i = 1; i < MAX_TRAVEL_TIME; i++)
			{
				P.JourneyDurationDistrib[i] += P.JourneyDurationDistrib[i - 1];
				P.LocalJourneyDurationDistrib[i] += P.LocalJourneyDurationDistrib[i - 1];
			}
			for (i = j = 0; i <= 1024; i++)
			{
				s = ((double)i) / 1024;
				while (P.JourneyDurationDistrib[j] < s)j++;
				P.InvJourneyDurationDistrib[i] = j;
			}
			for (i = j = 0; i <= 1024; i++)
			{
				s = ((double)i) / 1024;
				while (P.LocalJourneyDurationDistrib[j] < s)j++;
				P.InvLocalJourneyDurationDistrib[i] = j;
			}
		}
		GetInputParameter(dat, dat2, "Mean size of place types", "%lf", (void*)P.PlaceTypeMeanSize, P.PlaceTypeNum, 1, 0);
		GetInputParameter(dat, dat2, "Param 1 of place group size distribution", "%lf", (void*)P.PlaceTypeGroupSizeParam1, P.PlaceTypeNum, 1, 0);
		if (!GetInputParameter2(dat, dat2, "Power of place size distribution", "%lf", (void*)P.PlaceTypeSizePower, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeSizePower[i] = 0;
		//added to enable lognormal distribution - ggilani 09/02/17
		if (!GetInputParameter2(dat, dat2, "Standard deviation of place size distribution", "%lf", (void*)P.PlaceTypeSizeSD, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeSizeSD[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Offset of place size distribution", "%lf", (void*)P.PlaceTypeSizeOffset, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeSizeOffset[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Maximum of place size distribution", "%lf", (void*)P.PlaceTypeSizeMax, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeSizeMax[i] = 1e20;
		GetInputParameter(dat, dat2, "Proportion of between group place links", "%lf", (void*)P.PlaceTypePropBetweenGroupLinks, P.PlaceTypeNum, 1, 0);
		if (!GetInputParameter2(dat, dat2, "Kernel type for place types", "%i", (void*)P.PlaceTypeKernelType, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeKernelType[i] = P.MoveKernelType;
		GetInputParameter(dat, dat2, "Place overlap matrix", "%lf", (void*)P.PlaceExclusivityMatrix, P.PlaceTypeNum * P.PlaceTypeNum, 1, 0); //changed from P.PlaceTypeNum,P.PlaceTypeNum,0);
/* Note P.PlaceExclusivityMatrix not used at present - places assumed exclusive (each person belongs to 0 or 1 place) */
		GetInputParameter(dat, dat2, "Relative transmission rates for place types", "%lf", (void*)P.PlaceTypeTrans, P.PlaceTypeNum, 1, 0);
		for (i = 0; i < P.PlaceTypeNum; i++) P.PlaceTypeTrans[i] *= AgeSuscScale;
	}
	if (!GetInputParameter2(dat, dat2, "Daily seasonality coefficients", "%lf", (void*)P.Seasonality, DAYS_PER_YEAR, 1, 0))
	{
		P.DoSeasonality = 0;
		for (i = 0; i < DAYS_PER_YEAR; i++)
			P.Seasonality[i] = 1;
	}
	else
	{
		P.DoSeasonality = 1;
		s = 0;
		for (i = 0; i < DAYS_PER_YEAR; i++)
			s += P.Seasonality[i];
		s += 1e-20;
		s /= DAYS_PER_YEAR;
		for (i = 0; i < DAYS_PER_YEAR; i++)
			P.Seasonality[i] /= s;
	}
	if (!GetInputParameter2(dat, dat2, "Mutation probability per infection of step increase in infectiousness", "%lf", (void*)&(P.EvolInfectMutationRate), 1, 1, 0)) P.EvolInfectMutationRate = 0;
	if (!GetInputParameter2(dat, dat2, "Step size of increase in infectiousness due to mutation", "%lf", (void*)&(P.EvolInfectStep), 1, 1, 0)) P.EvolInfectStep = 0;
	if (!GetInputParameter2(dat, dat2, "Maximum increase in infectiousness due to mutation", "%lf", (void*)&(P.EvolInfectMax), 1, 1, 0)) P.EvolInfectMax = 0;

	if (!GetInputParameter2(dat, dat2, "Number of seed locations", "%i", (void*)&(P.NumSeedLocations), 1, 1, 0)) P.NumSeedLocations = 1;
	if (P.NumSeedLocations > MAX_NUM_SEED_LOCATIONS)
	{
		fprintf(stderr, "Too many seed locations\n");
		P.NumSeedLocations = MAX_NUM_SEED_LOCATIONS;
	}
	GetInputParameter(dat, dat2, "Initial number of infecteds", "%i", (void*)P.NumInitialInfections, P.NumSeedLocations, 1, 0);
	GetInputParameter(dat, dat2, "Location of initial infecteds", "%lf", (void*)&(P.LocationInitialInfection[0][0]), P.NumSeedLocations * 2, 1, 0);
	if (!GetInputParameter2(dat, dat2, "Minimum population in microcell of initial infection", "%i", (void*)&(P.MinPopDensForInitialInfection), 1, 1, 0)) P.MinPopDensForInitialInfection = 0;
	GetInputParameter(dat, dat2, "Maximum population in microcell of initial infection", "%i", (void*)&(P.MaxPopDensForInitialInfection), 1, 1, 0);
	GetInputParameter(dat, dat2, "Randomise initial infection location", "%i", (void*)&(P.DoRandomInitialInfectionLoc), 1, 1, 0);
	GetInputParameter(dat, dat2, "All initial infections located in same microcell", "%i", (void*)&(P.DoAllInitialInfectioninSameLoc), 1, 1, 0);
	if (P.DoAdUnits)
	{
		if (!GetInputParameter2(dat, dat2, "Administrative unit to seed initial infection into", "%i", (void*)&(P.InitialInfectionsAdminUnit[0]), P.NumSeedLocations, 1, 0))
			for (i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnit[i] = 0;
	}
	else
	{
		for (i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnit[i] = 0;
	}
	GetInputParameter(dat, dat2, "Reproduction number", "%lf", (void*)&(P.R0), 1, 1, 0);
	GetInputParameter(dat, dat2, "Infectious period", "%lf", (void*)&(P.InfectiousPeriod), 1, 1, 0);
	if (!GetInputParameter2(dat, dat2, "Assume SIS model", "%i", (void*)&(P.DoSIS), 1, 1, 0)) P.DoSIS = 0;
	if (!GetInputParameter2(dat, dat2, "SD of individual variation in infectiousness", "%lf", (void*)&(P.InfectiousnessSD), 1, 1, 0)) P.InfectiousnessSD = 0;
	if (P.InfectiousnessSD > 0)
	{
		P.InfectiousnessGamA = P.InfectiousnessGamR = 1 / (P.InfectiousnessSD * P.InfectiousnessSD);
		//GetInputParameter(dat,dat2,"Infectiousness Beta A","%lf",(void *) &(P.InfectiousnessBetaA),1,1,0);
		//GetInputParameter(dat,dat2,"Infectiousness Beta B","%lf",(void *) &(P.InfectiousnessBetaB),1,1,0);
	}
	if (P.DoSIS)
	{
		if (!GetInputParameter2(dat, dat2, "Factor by which susceptibility is multiplied per infection", "%lf", (void*)&(P.SuscReductionFactorPerInfection), 1, 1, 0)) P.SuscReductionFactorPerInfection = 1;
	}
	else
		P.SuscReductionFactorPerInfection = 0;

	if (!GetInputParameter2(dat, dat2, "Output R0 file", "%i", (void*)&(P.DoROutput), 1, 1, 0)) P.DoROutput = 0;
	if (!GetInputParameter2(dat, dat2, "Output infection type file", "%i", (void*)&(P.DoInftypeOutput), 1, 1, 0)) P.DoInftypeOutput = 0;

	if (!GetInputParameter2(dat, dat2, "Model time varying infectiousness", "%i", (void*)&(P.DoInfectiousnessProfile), 1, 1, 0)) P.DoInfectiousnessProfile = 0;
	if (!GetInputParameter2(dat, dat2, "Power of scaling of spatial R0 with density", "%lf", (void*)&(P.R0DensityScalePower), 1, 1, 0)) P.R0DensityScalePower = 0;
	if (P.DoInfectiousnessProfile)
	{
		if (!GetInputParameter2(dat, dat2, "Infectiousness profile", "%lf", (void*)P.infectious_prof, INFPROF_RES, 1, 0))
		{
			for (i = 0; i < INFPROF_RES; i++)
				P.infectious_prof[i] = 1;
		}
		k = (int)ceil(P.InfectiousPeriod / P.TimeStep);
		if (k >= MAX_INFECTIOUS_STEPS) ERR_CRITICAL("MAX_INFECTIOUS_STEPS not big enough\n");
		s = 0;
		P.infectious_prof[INFPROF_RES] = 0;
		for (i = 0; i < k; i++)
		{
			t = (((double)i) * P.TimeStep / P.InfectiousPeriod * INFPROF_RES);
			j = (int)t;
			t -= (double)j;
			if (j < INFPROF_RES)
				s += (P.infectiousness[i] = P.infectious_prof[j] * (1 - t) + P.infectious_prof[j + 1] * t);
			else
				s += (P.infectiousness[i] = P.infectious_prof[INFPROF_RES]);
		}
		s /= ((double)k);
		for (i = 0; i <= k; i++) P.infectiousness[i] /= s;
		for (i = 0; i <= CDF_RES; i++) P.infectious_icdf[i] = exp(-1.0);
	}
	else
	{
		if (!GetInputParameter2(dat, dat2, "Infectious period inverse CDF", "%lf", (void*)P.infectious_icdf, CDF_RES + 1, 1, 0))
		{
			P.infectious_icdf[CDF_RES] = 100;
			for (i = 0; i < CDF_RES; i++)
				P.infectious_icdf[i] = -log(1 - ((double)i) / CDF_RES);
		}
		k = (int)ceil(P.InfectiousPeriod * P.infectious_icdf[CDF_RES] / P.TimeStep);
		if (k >= MAX_INFECTIOUS_STEPS) ERR_CRITICAL("MAX_INFECTIOUS_STEPS not big enough\n");
		for (i = 0; i < k; i++) P.infectiousness[i] = 1.0;
		P.infectiousness[k] = 0;
		for (i = 0; i <= CDF_RES; i++) P.infectious_icdf[i] = exp(-P.infectious_icdf[i]);
	}
	//if(!GetInputParameter2(dat,dat2,"Proportion of cases hospitalised","%lf",(void *) &(P.ProportionHospitalised),1,1,0)) P.ProportionHospitalised=0;
	if (!GetInputParameter2(dat, dat2, "Include mortality", "%i", (void*)&(P.DoMortality), 1, 1, 0)) P.DoMortality = 0;
	if (P.DoMortality) //changes to mortality to allow for probability of death/recovery to be decided when leaving the infectious class: ggilani - 23/10/14
	{
		if (!GetInputParameter2(dat, dat2, "Include event-dependent mortality", "%i", (void*)&(P.DoEventMortality), 1, 1, 0)) P.DoEventMortality = 0;
		if (P.DoEventMortality)
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
	if (!GetInputParameter2(dat, dat2, "Include funeral transmission", "%i", (void*)&(P.DoFuneralTransmission), 1, 1, 0)) P.DoFuneralTransmission = 0;
	if (P.DoFuneralTransmission) //Funeral transmission parameters: ggilani - 26/10/14
	{
		if (!GetInputParameter2(dat, dat2, "Duration of funeral transmission", "%lf", (void*)&(P.FuneralTransmissionDuration), 1, 1, 0)) P.FuneralTransmissionDuration = 0;
		if (!GetInputParameter2(dat, dat2, "Relative infectiousness of a funeral transmission", "%lf", (void*)&(P.RelativeInfectiousnessFuneral), 1, 1, 0)) P.RelativeInfectiousnessFuneral = 1;
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
	if (!GetInputParameter2(dat, dat2, "Include hospitalisation", "%i", (void*)&(P.DoHospitalisation), 1, 1, 0)) P.DoHospitalisation = 0; //Hospitalisation parameters: ggilani 28/10/14
	if (P.DoHospitalisation)
	{

		//leave the terminology the same at the moment but consider making these more general. These relate to both hospitals and ETUs
		if (!GetInputParameter2(dat, dat2, "Mean time to hospitalisation", "%lf", (void*)&(P.HospitalisationTime), 1, 1, 0)) P.HospitalisationTime = 0;
		if (!GetInputParameter2(dat, dat2, "Hospital waiting time", "%lf", (void*)&(P.HospWaitingTime), 1, 1, 0)) P.HospWaitingTime = 0.25; //To ensure acceptance 
		if (!GetInputParameter2(dat, dat2, "Time to hospitalisation inverse CDF", "%lf", (void*)P.hospital_icdf, CDF_RES + 1, 1, 0))
		{
			P.hospital_icdf[CDF_RES] = 1e10;
			for (i = 0; i < CDF_RES; i++)
				P.hospital_icdf[i] = -log(1 - ((double)i) / CDF_RES);
		}
		for (i = 0; i <= CDF_RES; i++)
			P.hospital_icdf[i] = exp(-P.hospital_icdf[i]);
		if (!GetInputParameter2(dat, dat2, "Number of times to hospitalisation", "%i", (void*)&(P.NMeanTimeToHosp), 1, 1, 0)) P.NMeanTimeToHosp = 0;
		if (P.NMeanTimeToHosp > 0)
		{
			if (P.NMeanTimeToHosp >= MAX_CHANGE_POINTS) ERR_CRITICAL("MAX_CHANGE_POINTS too small\n");
			GetInputParameter(dat, dat2, "Change points times to hospitalisation", "%lf", (void*)P.ChangePointMeanTimeToHosp, P.NMeanTimeToHosp, 1, 0);
			GetInputParameter(dat, dat2, "Times to hospitalisation", "%lf", (void*)P.MeanTimeToHosp, P.NMeanTimeToHosp, 1, 0);
		}
		P.CurrIndMeanTimeToHosp = 0;

		//
		if (!GetInputParameter2(dat, dat2, "Proportion of cases seeking care before outbreak declared", "%lf", (void*)&P.PropHospSeekPreOutbreak, 1, 1, 0)) P.PropHospSeekPreOutbreak = 0.5;
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
		if (P.DoAdUnits)
		{
			int AdunitETUCapacity[MAX_ADUNITS];
			double AdunitETUTime[MAX_ADUNITS]; //time when beds become available in each admin unit
			if (!GetInputParameter2(dat, dat2, "Do ETUs by admin unit", "%i", (void*)&(P.DoETUByAdUnit), 1, 1, 0)) P.DoETUByAdUnit = 0;

			if (P.DoETUByAdUnit)
			{
				if (!GetInputParameter2(dat, dat2, "Initial bed capacity per ETU admin unit", "%i", (void*)AdunitETUCapacity, P.NumAdunits, 1, 0))
				{
					for (i = 0; i < P.NumAdunits; i++) AdunitETUCapacity[i] = 0;
				}
				if (!GetInputParameter2(dat, dat2, "Time ETU beds available per admin unit", "%lf", (void*)AdunitETUTime, P.NumAdunits, 1, 0))
				{
					for (i = 0; i < P.NumAdunits; i++) AdunitETUTime[i] = 1e10;
				}

				//save these to admin unit variables
				for (i = 0; i < P.NumAdunits; i++)
				{
					//AdUnits[i].totalBeds=AdunitHospCapacity[i];
					AdUnits[i].totalETUBeds = 0;
					AdUnits[i].initialETUBeds = AdunitETUCapacity[i];
					AdUnits[i].timeETUBedsAvailable = AdunitETUTime[i];
					AdUnits[i].currentETUBeds = 0;
					AdUnits[i].nextETUBeds = 0;
					AdUnits[i].nextTimeToETUBeds = 1e10;
					AdUnits[i].flagETUBeds = 0;
				}

				if (!GetInputParameter2(dat, dat2, "Do reactive ETU bed allocation", "%i", (void*)&(P.DoReactETUBeds), 1, 1, 0)) P.DoReactETUBeds = 0;

				if (P.DoReactETUBeds)
				{
					// for reactive provisioning of beds
					if (!GetInputParameter2(dat, dat2, "Reactive ETU bed allocation trigger incidence per cell", "%lf", (void*)&(P.ETUCellIncThresh), 1, 1, 0)) P.ETUCellIncThresh = 1000000000;
					if (!GetInputParameter2(dat, dat2, "Reactive ETU bed allocation start time", "%lf", (void*)&(P.ETUTimeStartBase), 1, 1, 0)) P.ETUTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
					//if (!GetInputParameter2(dat, dat2, "Time to start reactive bed policy", "%i", (void *) &(P.StartTimeReactiveBeds), 1, 1, 0)) P.StartTimeReactiveBeds= 1e6;
					if (!GetInputParameter2(dat, dat2, "Initial number of detected cases to trigger ETU beds in an admin unit", "%i", (void*)&(P.InitCasesToETUBeds), 1, 1, 0)) P.InitCasesToETUBeds = 1e6;
					if (!GetInputParameter2(dat, dat2, "Initial time to install ETU beds", "%lf", (void*)&(P.InitDelayToETUBeds), 1, 1, 0)) P.InitDelayToETUBeds = 0;
					if (!GetInputParameter2(dat, dat2, "Initial number of ETU beds", "%i", (void*)&(P.InitNumETUBeds), 1, 1, 0)) P.InitNumETUBeds = 0;
					if (!GetInputParameter2(dat, dat2, "Hospital capacity reached before increasing ETU bed numbers", "%lf", (void*)&(P.CapacityToMoreETUBeds), 1, 1, 0)) P.CapacityToMoreETUBeds = 0;
					if (!GetInputParameter2(dat, dat2, "Subsequent time to install ETU beds", "%lf", (void*)&(P.SubDelayToETUBeds), 1, 1, 0)) P.SubDelayToETUBeds = 0;
					if (!GetInputParameter2(dat, dat2, "Subsequent number of additional ETU beds", "%i", (void*)&(P.SubNumETUBeds), 1, 1, 0)) P.SubNumETUBeds = 0;
				}
				if (!GetInputParameter2(dat, dat2, "Output ETU capacity", "%i", (void*)&(P.DoOutputETUCapacity), 1, 1, 0)) P.DoOutputETUCapacity = 0;
			}
		}
		else P.DoETUByAdUnit = 0;

		//This section is going to be specific to hospitals rather than ETUs
		if (!GetInputParameter2(dat, dat2, "Include hospitals as places", "%i", (void*)&(P.IncludeHospitalPlaceType), 1, 1, 0)) P.IncludeHospitalPlaceType = 0;
		if (P.IncludeHospitalPlaceType)
		{
			if (!GetInputParameter2(dat, dat2, "Output keyworker file", "%i", (void*)&(P.DoKeyworkerOutput), 1, 1, 0)) P.DoKeyworkerOutput = 0;
			GetInputParameter(dat, dat2, "Hospital place type number", "%i", (void*)&(P.HospPlaceTypeNum), 1, 1, 0);
			GetInputParameter(dat, dat2, "Max number of beds for Ebola cases in hospital", "%i", (void*)&(P.HospCaseCapacity), 1, 1, 0);
			GetInputParameter(dat, dat2, "Number of HCWs per 1000 population", "%lf", (void*)&(P.HCWPerThousand), 1, 1, 0);
			GetInputParameter(dat, dat2, "Proportion of HCWs and FLWs vaccinated before outbreak", "%lf", (void*)&(P.PropHCWFLWVacc), 1, 1, 0);
			if (!GetInputParameter2(dat, dat2, "Relative susceptibility of HCWs and FLWs vaccinated before outbreak", "%lf", (void*)&(P.HCWVaccSuscDrop), 1, 1, 0)) P.HCWVaccSuscDrop = -1;
			if (!GetInputParameter2(dat, dat2, "Revaccinate HCW and FLWs after outbreak alert", "%i", (void*)&(P.RevaccHCWs), 1, 1, 0)) P.RevaccHCWs = 0;
			GetInputParameter(dat, dat2, "Number of days between HCW and FLW vaccination and outbreak", "%i", (void*)&(P.DayHCWFLWVacc), 1, 1, 0);
			GetInputParameter(dat, dat2, "Min age for HCWs and FLWs", "%i", (void*)&(P.MinAgeHCWFLW), 1, 1, 0);
			GetInputParameter(dat, dat2, "Max age for HCWs and FLWs", "%i", (void*)&(P.MaxAgeHCWFLW), 1, 1, 0);
			if (!GetInputParameter2(dat, dat2, "Include FLWs", "%i", (void*)&(P.IncludeFLWs), 1, 1, 0)) P.IncludeFLWs = 0;
			if (P.IncludeFLWs)
			{
				GetInputParameter(dat, dat2, "Number of FLWs per 1000 population", "%lf", (void*)&(P.FLWPerThousand), 1, 1, 0);
				GetInputParameter(dat, dat2, "Relative suspectibility of FLWs", "%lf", (void*)&(P.RelSuscFLW), 1, 1, 0);
			}
			if (!GetInputParameter2(dat, dat2, "Relative susceptibility of HCWs and FLWs due to PPE post outbreak detection", "%lf", (void*&)P.RelSuscPPE, 1, 1, 0)) P.RelSuscPPE = 1;
		}

		if (!GetInputParameter2(dat, dat2, "Relative infectiousness of an ETU case", "%lf", (void*)&(P.RelativeInfectiousnessETU), 1, 1, 0)) P.RelativeInfectiousnessETU = 1;
	}

	if (!GetInputParameter2(dat, dat2, "Do interrupt interventions", "%i", (void*)&(P.DoInterruptIntervention), 1, 1, 0)) P.DoInterruptIntervention = 0; //Interruption of interventions parameters: ggilani 28/10/14
	if (P.DoInterruptIntervention)
	{
		if (!GetInputParameter2(dat, dat2, "Number of days of interruptions", "%i", (void*)&(P.NDaysInterrupt), 1, 1, 0)) P.NDaysInterrupt = 0;
		if (P.NDaysInterrupt > 0)
		{
			if (P.NDaysInterrupt >= MAX_CHANGE_POINTS) ERR_CRITICAL("MAX_CHANGE_POINTS too small\n");
			GetInputParameter(dat, dat2, "Days of interruptions", "%i", (void*)&(P.DaysInterruptIntervention), P.NDaysInterrupt, 1, 0);
		}
	}

	if (!GetInputParameter2(dat, dat2, "Do prevalence dependent transmission", "%i", (void*)&(P.DoPrevDependTrans), 1, 1, 0)) P.DoPrevDependTrans = 0;
	//added prevalent dependent transmission - ggilani, 22/12/14
	if (P.DoPrevDependTrans)
	{
		// determine type of prevalence dependent transmission - either by total number of cases so far (which will never go down i.e. once we're over the threshold, people's behaviour is permanently changed) or by current number of infections (i.e. it can go down again)
		if (!GetInputParameter2(dat, dat2, "Do prevalence dependent transmission by total cases", "%i", (void*)&(P.DoPrevDepTransTotalCases), 1, 1, 0)) P.DoPrevDepTransTotalCases = 0;
		if (!GetInputParameter2(dat, dat2, "Do prevalence dependent transmission by current infections", "%i", (void*)&(P.DoPrevDepTransCurrInf), 1, 1, 0)) P.DoPrevDepTransCurrInf = 0;
		// we can only do one version though
		if ((P.DoPrevDepTransTotalCases) && (P.DoPrevDepTransCurrInf))
		{
			ERR_CRITICAL("Only one type of prevalent dependent transmission can be specified!\n");
		}
		if (!GetInputParameter2(dat, dat2, "Threshold for prevalence dependent transmission", "%lf", (void*)&(P.PrevDependTransThresh), 1, 1, 0)) P.PrevDependTransThresh = 1.0;
		if (!GetInputParameter2(dat, dat2, "Threshold for prevalence dependent transmission", "%lf", (void*)&(P.PrevDependRelativeSusc), 1, 1, 0)) P.PrevDependRelativeSusc = 1.0;

	}
	if (!GetInputParameter2(dat, dat2, "Include latent period", "%i", (void*)&(P.DoLatent), 1, 1, 0)) P.DoLatent = 0;
	if (P.DoLatent)
	{
		GetInputParameter(dat, dat2, "Latent period", "%lf", (void*)&(P.LatentPeriod), 1, 1, 0);
		if (!GetInputParameter2(dat, dat2, "Latent period inverse CDF", "%lf", (void*)P.latent_icdf, CDF_RES + 1, 1, 0))
		{
			P.latent_icdf[CDF_RES] = 1e10;
			for (i = 0; i < CDF_RES; i++)
				P.latent_icdf[i] = -log(1 - ((double)i) / CDF_RES);
		}
		for (i = 0; i <= CDF_RES; i++)
			P.latent_icdf[i] = exp(-P.latent_icdf[i]);
	}

	if (!GetInputParameter2(dat, dat2, "Include symptoms", "%i", (void*)&(P.DoSymptoms), 1, 1, 0)) P.DoSymptoms = 0;
	if (!P.DoSymptoms)
	{
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.ProportionSymptomatic[i] = 0;
		P.FalsePositiveRate = 0;
		P.SymptInfectiousness = 1.0;
		P.LatentToSymptDelay = 0;
	}
	else
	{
		if (P.DoAge)
			GetInputParameter(dat, dat2, "Proportion symptomatic by age group", "%lf", (void*)P.ProportionSymptomatic, NUM_AGE_GROUPS, 1, 0);
		else
		{
			GetInputParameter(dat, dat2, "Proportion symptomatic", "%lf", (void*)P.ProportionSymptomatic, 1, 1, 0);
			for (i = 1; i < NUM_AGE_GROUPS; i++)
				P.ProportionSymptomatic[i] = P.ProportionSymptomatic[0];
		}
		GetInputParameter(dat, dat2, "Delay from end of latent period to start of symptoms", "%lf", (void*)&(P.LatentToSymptDelay), 1, 1, 0);
		GetInputParameter(dat, dat2, "Relative rate of random contacts if symptomatic", "%lf", (void*)&(P.SymptSpatialContactRate), 1, 1, 0);
		GetInputParameter(dat, dat2, "Symptomatic infectiousness relative to asymptomatic", "%lf", (void*)&(P.SymptInfectiousness), 1, 1, 0);
		if (!GetInputParameter2(dat, dat2, "Model symptomatic withdrawal to home as true absenteeism", "%i", (void*)&P.DoRealSymptWithdrawal, 1, 1, 0)) P.DoRealSymptWithdrawal = 0;
		if (P.DoPlaces)
		{
			GetInputParameter(dat, dat2, "Relative level of place attendance if symptomatic", "%lf", (void*)P.SymptPlaceTypeContactRate, P.PlaceTypeNum, 1, 0);
			if (P.DoRealSymptWithdrawal)
			{
				for (j = 0; j < NUM_PLACE_TYPES; j++)
				{
					P.SymptPlaceTypeWithdrawalProp[j] = 1.0 - P.SymptPlaceTypeContactRate[j];
					P.SymptPlaceTypeContactRate[j] = 1.0;
				}
			}
			else
				for (j = 0; j < NUM_PLACE_TYPES; j++) P.SymptPlaceTypeWithdrawalProp[j] = 0.0;
		}
		if (!GetInputParameter2(dat, dat2, "Maximum age of child at home for whom one adult also stays at home", "%i", (void*)&P.CaseAbsentChildAgeCutoff, 1, 1, 0)) P.CaseAbsentChildAgeCutoff = 0;
		if (!GetInputParameter2(dat, dat2, "Proportion of children at home for whom one adult also stays at home", "%lf", (void*)&P.CaseAbsentChildPropAdultCarers, 1, 1, 0)) P.CaseAbsentChildPropAdultCarers = 0;
#ifdef ABSENTEEISM_PLACE_CLOSURE
		P.CaseAbsenteeismDelay = 0;  // Set to zero for tracking absenteeism
#else
		if (!GetInputParameter2(dat, dat2, "Delay in starting place absenteeism for cases who withdraw", "%lf", (void*)&P.CaseAbsenteeismDelay, 1, 1, 0)) P.CaseAbsenteeismDelay = 0;
#endif
		if (!GetInputParameter2(dat, dat2, "Duration of place absenteeism for cases who withdraw", "%lf", (void*)&P.CaseAbsenteeismDuration, 1, 1, 0)) P.CaseAbsenteeismDuration = 7;

		if (!GetInputParameter2(dat, dat2, "False positive rate", "%lf", (void*)&(P.FalsePositiveRate), 1, 1, 0)) P.FalsePositiveRate = 0.0;
		if (!GetInputParameter2(dat, dat2, "False positive per capita incidence", "%lf", (void*)&(P.FalsePositivePerCapitaIncidence), 1, 1, 0)) P.FalsePositivePerCapitaIncidence = 0.0;
		if (!GetInputParameter2(dat, dat2, "False positive relative incidence by age", "%lf", (void*)P.FalsePositiveAgeRate, NUM_AGE_GROUPS, 1, 0))
			for (j = 0; j < NUM_AGE_GROUPS; j++) P.FalsePositiveAgeRate[j] = 1.0;
	}
	if (!GetInputParameter2(dat, dat2, "Bounding box for bitmap", "%lf", (void*)&(P.BoundingBox[0]), 4, 1, 0))
	{
		P.BoundingBox[0] = P.BoundingBox[1] = 0.0;
		P.BoundingBox[2] = P.BoundingBox[3] = 1.0;
	}
	if (!GetInputParameter2(dat, dat2, "Spatial domain for simulation", "%lf", (void*)&(P.SpatialBoundingBox[0]), 4, 1, 0))
	{
		P.SpatialBoundingBox[0] = P.SpatialBoundingBox[1] = 0.0;
		P.SpatialBoundingBox[2] = P.SpatialBoundingBox[3] = 1.0;
	}
	if (!GetInputParameter2(dat, dat2, "Grid size", "%lf", (void*)&(P.cwidth), 1, 1, 0)) P.cwidth = 1.0 / 120.0;
	if (!GetInputParameter2(dat, dat2, "Use long/lat coord system", "%i", (void*)&(P.DoUTM_coords), 1, 1, 0)) P.DoUTM_coords = 0;
	if (!GetInputParameter2(dat, dat2, "Bitmap scale", "%lf", (void*)&(P.BitmapScale), 1, 1, 0)) P.BitmapScale = 1.0;
	if (!GetInputParameter2(dat, dat2, "Bitmap y:x aspect scaling", "%lf", (void*)&(P.BitmapAspectScale), 1, 1, 0)) P.BitmapAspectScale = 1.0;
	if (!GetInputParameter2(dat, dat2, "Bitmap movie frame interval", "%i", (void*)&(P.BitmapMovieFrame), 1, 1, 0)) P.BitmapMovieFrame = 250;
	if (!GetInputParameter2(dat, dat2, "Output bitmap", "%i", (void*)&(P.OutputBitmap), 1, 1, 0)) P.OutputBitmap = 0;
	if (!GetInputParameter2(dat, dat2, "Output bitmap detected", "%i", (void*)&(P.OutputBitmapDetected), 1, 1, 0)) P.OutputBitmapDetected = 0;
	if (!GetInputParameter2(dat, dat2, "Output immunity on bitmap", "%i", (void*)&(P.DoImmuneBitmap), 1, 1, 0)) P.DoImmuneBitmap = 0;
	if (!GetInputParameter2(dat, dat2, "Output infection tree", "%i", (void*)&(P.DoInfectionTree), 1, 1, 0)) P.DoInfectionTree = 0;
	if (!GetInputParameter2(dat, dat2, "Do one generation", "%i", (void*)&(P.DoOneGen), 1, 1, 0)) P.DoOneGen = 0;
	if (!GetInputParameter2(dat, dat2, "Output every realisation", "%i", (void*)&(P.OutputAll), 1, 1, 0)) P.OutputAll = 1;
	if (!GetInputParameter2(dat, dat2, "Maximum number to sample for correlations", "%i", (void*)&(P.MaxCorrSample), 1, 1, 0)) P.MaxCorrSample = 1000000000;
	if (!GetInputParameter2(dat, dat2, "Assume SI model", "%i", (void*)&(P.DoSI), 1, 1, 0)) P.DoSI = 0;
	if (!GetInputParameter2(dat, dat2, "Assume periodic boundary conditions", "%i", (void*)&(P.DoPeriodicBoundaries), 1, 1, 0)) P.DoPeriodicBoundaries = 0;
	if (!GetInputParameter2(dat, dat2, "Only output non-extinct realisations", "%i", (void*)&(P.OutputNonExtinct), 1, 1, 0)) P.OutputNonExtinct = 0;

	if (!GetInputParameter2(dat, dat2, "Output control file", "%i", (void*)&(P.DoControlOutput), 1, 1, 0)) P.DoControlOutput = 0;
	if (!GetInputParameter2(dat, dat2, "Use cases per thousand threshold for area controls", "%i", (void*)&(P.DoPerCapitaTriggers), 1, 1, 0)) P.DoPerCapitaTriggers = 0;
	if (!GetInputParameter2(dat, dat2, "Use global triggers for interventions", "%i", (void*)&(P.DoGlobalTriggers), 1, 1, 0)) P.DoGlobalTriggers = 0;
	if (!GetInputParameter2(dat, dat2, "Divisor for per-capita area threshold (default 1000)", "%i", (void*)&(P.IncThreshPop), 1, 1, 0)) P.IncThreshPop = 1000;
	if (!GetInputParameter2(dat, dat2, "Divisor for per-capita global threshold (default 1000)", "%i", (void*)&(P.GlobalIncThreshPop), 1, 1, 0)) P.GlobalIncThreshPop = 1000;


	if (!GetInputParameter2(dat, dat2, "Number of sampling intervals over which cumulative incidence measured for global trigger", "%i", (void*)&(P.TriggersSamplingInterval), 1, 1, 0)) P.TriggersSamplingInterval = 10000000;
	//if (!GetInputParameter2(dat, dat2, "Number of undetected infections before first case detected", "%i", (void*)&(P.NumUndetectedInfPreOutbreakAlert), 1, 1, 0)) P.NumUndetectedInfPreOutbreakAlert = 0;
	if (!GetInputParameter2(dat, dat2, "Proportion of cases detected after surveillance alert", "%lf", (void*)&(P.PostAlertControlPropCasesId), 1, 1, 0)) P.PostAlertControlPropCasesId = 1;
	//if(!GetInputParameter2(dat,dat2,"Proportion of cases detected before surveillance alert","%lf",(void *) &(P.PreAlertControlPropCasesId),1,1,0)) P.PreAlertControlPropCasesId=1;
	//if(P.PreControlClusterIdCaseThreshold==0)
	//	{
	if (!GetInputParameter2(dat, dat2, "Number of detected cases before outbreak alert triggered", "%i", (void*)&(P.PreControlClusterIdCaseThreshold), 1, 1, 0)) P.PreControlClusterIdCaseThreshold = 0;
	//	}
	if (!GetInputParameter2(dat, dat2, "Only use confirmed cases to trigger alert", "%i", (void*)&(P.DoEarlyCaseDiagnosis), 1, 1, 0)) P.DoEarlyCaseDiagnosis = 0;
	if (!GetInputParameter2(dat, dat2, "Target country", "%i", (void*)&(P.TargetCountry), 1, 1, 0)) P.TargetCountry = 1;
	fprintf(stderr, "Target country=%i\n", P.TargetCountry);
	//added extra target countries to help keep track of country wise parameters: ggilani
	if (!GetInputParameter2(dat, dat2, "Target country 2", "%i", (void*)&(P.TargetCountry2), 1, 1, 0)) P.TargetCountry2 = 1;
	if (!GetInputParameter2(dat, dat2, "Target country 3", "%i", (void*)&(P.TargetCountry3), 1, 1, 0)) P.TargetCountry3 = 1;
	if (!GetInputParameter2(dat, dat2, "Restrict treatment to target country", "%i", (void*)&(P.RestrictTreatToTarget), 1, 1, 0)) P.RestrictTreatToTarget = 0;
	if (!GetInputParameter2(dat, dat2, "Only treat mixing groups within places", "%i", (void*)&(P.DoPlaceGroupTreat), 1, 1, 0)) P.DoPlaceGroupTreat = 0;

	if (!GetInputParameter2(dat, dat2, "Treatment trigger incidence per cell", "%lf", (void*)&(P.TreatCellIncThresh), 1, 1, 0)) P.TreatCellIncThresh = 1000000000;
	if (!GetInputParameter2(dat, dat2, "Relative susceptibility of treated individual", "%lf", (void*)&(P.TreatSuscDrop), 1, 1, 0)) P.TreatSuscDrop = 1;
	if (!GetInputParameter2(dat, dat2, "Relative infectiousness of treated individual", "%lf", (void*)&(P.TreatInfDrop), 1, 1, 0)) P.TreatInfDrop = 1;
	if (!GetInputParameter2(dat, dat2, "Proportion of symptomatic cases resulting in death prevented by treatment", "%lf", (void*)&(P.TreatDeathDrop), 1, 1, 0)) P.TreatDeathDrop = 0;
	if (!GetInputParameter2(dat, dat2, "Proportion of symptomatic cases prevented by treatment", "%lf", (void*)&(P.TreatSympDrop), 1, 1, 0)) P.TreatSympDrop = 0;
	if (!GetInputParameter2(dat, dat2, "Delay to treat cell", "%lf", (void*)&(P.TreatDelayMean), 1, 1, 0)) P.TreatDelayMean = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of course of treatment", "%lf", (void*)&(P.TreatCaseCourseLength), 1, 1, 0)) P.TreatCaseCourseLength = 5;
	if (!GetInputParameter2(dat, dat2, "Duration of course of prophylaxis", "%lf", (void*)&(P.TreatProphCourseLength), 1, 1, 0)) P.TreatProphCourseLength = 10;
	if (!GetInputParameter2(dat, dat2, "Proportion of detected cases treated", "%lf", (void*)&(P.TreatPropCases), 1, 1, 0)) P.TreatPropCases = 1;
	if (!GetInputParameter2(dat, dat2, "Proportion of detected cases with private stockpile treated", "%lf", (void*)&(P.PrivateTreatPropCases), 1, 1, 0)) P.PrivateTreatPropCases = 0;
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(dat, dat2, "Delay to treat households with private stockpiles", "%lf", (void*)&(P.PrivateTreatDelayMean), 1, 1, 0)) P.PrivateTreatDelayMean = 0;
		if (!GetInputParameter2(dat, dat2, "Proportion of households of cases treated", "%lf", (void*)&(P.TreatPropCaseHouseholds), 1, 1, 0)) P.TreatPropCaseHouseholds = 0;
		if (!GetInputParameter2(dat, dat2, "Proportion of households with private stockpiles treated", "%lf", (void*)&(P.PrivateTreatPropCaseHouseholds), 1, 1, 0)) P.PrivateTreatPropCaseHouseholds = 0;
		if (!GetInputParameter2(dat, dat2, "Duration of household prophylaxis policy", "%lf", (void*)&(P.TreatHouseholdsDuration), 1, 1, 0)) P.TreatHouseholdsDuration = USHRT_MAX / P.TimeStepsPerDay;
	}
	if (!GetInputParameter2(dat, dat2, "Proportion treated", "%lf", (void*)&(P.TreatPropRadial), 1, 1, 0)) P.TreatPropRadial = 1.0;
	if (!GetInputParameter2(dat, dat2, "Proportion treated in radial prophylaxis", "%lf", (void*)&(P.TreatPropRadial), 1, 1, 0)) P.TreatPropRadial = 1.0;
	if (!GetInputParameter2(dat, dat2, "Treatment radius", "%lf", (void*)&(P.TreatRadius), 1, 1, 0)) P.TreatRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of place/geographic prophylaxis policy", "%lf", (void*)&(P.TreatPlaceGeogDuration), 1, 1, 0)) P.TreatPlaceGeogDuration = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Treatment start time", "%lf", (void*)&(P.TreatTimeStartBase), 1, 1, 0)) P.TreatTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(dat, dat2, "Proportion of places treated after case detected", "%lf", (void*)P.TreatPlaceProbCaseId, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.TreatPlaceProbCaseId[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Proportion of people treated in targeted places", "%lf", (void*)P.TreatPlaceTotalProp, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.TreatPlaceTotalProp[i] = 0;
	}
	if (!GetInputParameter2(dat, dat2, "Maximum number of doses available", "%lf", (void*)&(P.TreatMaxCoursesBase), 1, 1, 0)) P.TreatMaxCoursesBase = 1e20;
	if (!GetInputParameter2(dat, dat2, "Start time of additional treatment production", "%lf", (void*)&(P.TreatNewCoursesStartTime), 1, 1, 0)) P.TreatNewCoursesStartTime = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Rate of additional treatment production (courses per day)", "%lf", (void*)&(P.TreatNewCoursesRate), 1, 1, 0)) P.TreatNewCoursesRate = 0;
	if (!GetInputParameter2(dat, dat2, "Maximum number of people targetted with radial prophylaxis per case", "%i", (void*)&(P.TreatMaxCoursesPerCase), 1, 1, 0)) P.TreatMaxCoursesPerCase = 1000000000;
	if (!GetInputParameter2(dat, dat2, "Proportion of households buying private stockpile", "%lf", (void*)&(P.PropPrivateStockpile), 1, 1, 0)) P.PropPrivateStockpile = 0;
	if (!GetInputParameter2(dat, dat2, "Income ordering of private stockpile (-1=low to high, 0=random, 1=high to low)", "%i", (void*)&(P.PrivateStockpileOrderByIncome), 1, 1, 0)) P.PrivateStockpileOrderByIncome = 0;

	if (!GetInputParameter2(dat, dat2, "Number of resistance levels", "%i", (void*)&(P.EvolResistNumTypes), 1, 1, 0)) P.EvolResistNumTypes = 0;
	P.EvolResistNumTypes++;
	if (P.EvolResistNumTypes > MAX_NUM_RESIST_TYPES) ERR_CRITICAL("Too many resistance levels\n");
	if (!GetInputParameter2(dat, dat2, "Mutation probability per treated infection of step increase in resistance", "%lf", (void*)&(P.EvolResistTreatMutationRate), 1, 1, 0)) P.EvolResistTreatMutationRate = 0;
	if (!GetInputParameter2(dat, dat2, "Mutation probability per prophylaxed infection of step increase in resistance", "%lf", (void*)&(P.EvolResistProphMutationRate), 1, 1, 0)) P.EvolResistProphMutationRate = 0;
	P.EvolResistRelInf[0] = 1.0;
	P.EvolResistRelTreatInfDrop[0] = 0.0;
	P.EvolResistRelProphSusc[0] = 0.0;
	P.EvolResistRelTreatSympDrop[0] = 0.0;
	P.EvolResistRelTreatDeathDrop[0] = 0.0;
	P.EvolResistSeedProp[0] = 1.0;
	if (!GetInputParameter2(dat, dat2, "Proportion of seed infections with different resistance levels", "%lf", (void*)(P.EvolResistSeedProp + 1), P.EvolResistNumTypes - 1, 1, 0))
	{
		for (i = 1; i < P.EvolResistNumTypes; i++) P.EvolResistSeedProp[i] = 0.0;
	}
	else
	{
		for (i = 1; i < P.EvolResistNumTypes; i++) P.EvolResistSeedProp[0] -= P.EvolResistSeedProp[i];
	}
	if (!GetInputParameter2(dat, dat2, "Relative infectiousness of resistance levels", "%lf", (void*)(P.EvolResistRelInf + 1), P.EvolResistNumTypes - 1, 1, 0))
	{
		for (i = 1; i < P.EvolResistNumTypes; i++) P.EvolResistRelInf[i] = 1.0;
	}
	if (!GetInputParameter2(dat, dat2, "Reduction in drop in infectiousness from treatment caused by resistance levels", "%lf", (void*)(P.EvolResistRelTreatInfDrop + 1), P.EvolResistNumTypes - 1, 1, 0))
	{
		for (i = 1; i < P.EvolResistNumTypes; i++) P.EvolResistRelTreatInfDrop[i] = 0.0;
	}
	if (!GetInputParameter2(dat, dat2, "Reduction in drop in susceptiblity from prophylaxis caused by resistance levels", "%lf", (void*)(P.EvolResistRelProphSusc + 1), P.EvolResistNumTypes - 1, 1, 0))
	{
		for (i = 1; i < P.EvolResistNumTypes; i++) P.EvolResistRelProphSusc[i] = 0.0;
	}
	if (!GetInputParameter2(dat, dat2, "Reduction in drop in symptoms from treatment caused by resistance levels", "%lf", (void*)(P.EvolResistRelTreatSympDrop + 1), P.EvolResistNumTypes - 1, 1, 0))
	{
		for (i = 1; i < P.EvolResistNumTypes; i++) P.EvolResistRelTreatSympDrop[i] = 0.0;
	}
	if (!GetInputParameter2(dat, dat2, "Reduction in drop in deaths from treatment caused by resistance levels", "%lf", (void*)(P.EvolResistRelTreatDeathDrop + 1), P.EvolResistNumTypes - 1, 1, 0))
	{
		for (i = 1; i < P.EvolResistNumTypes; i++) P.EvolResistRelTreatDeathDrop[i] = 0.0;
	}
	for (i = 0; i < P.EvolResistNumTypes; i++)
	{
		P.EvolResistRelTreatInfDrop[i] = 1 - (1 - P.TreatInfDrop) * (1 - P.EvolResistRelTreatInfDrop[i]);
		P.EvolResistRelProphSusc[i] = 1 - (1 - P.TreatSuscDrop) * (1 - P.EvolResistRelProphSusc[i]);
		P.EvolResistRelTreatSympDrop[i] = 1 - P.TreatSympDrop * (1 - P.EvolResistRelTreatSympDrop[i]);
		P.EvolResistRelTreatDeathDrop[i] = 1 - P.TreatDeathDrop * (1 - P.EvolResistRelTreatDeathDrop[i]);
	}
	if (P.DoAdUnits)
	{
		if (!GetInputParameter2(dat, dat2, "Treat administrative units rather than rings", "%i", (void*)&(P.TreatByAdminUnit), 1, 1, 0)) P.TreatByAdminUnit = 0;
		if (!GetInputParameter2(dat, dat2, "Administrative unit divisor for treatment", "%i", (void*)&(P.TreatAdminUnitDivisor), 1, 1, 0)) P.TreatAdminUnitDivisor = 1;
		if ((P.TreatAdminUnitDivisor == 0) || (P.TreatByAdminUnit == 0)) { P.TreatByAdminUnit = 0; P.TreatAdminUnitDivisor = 1; }
	}
	else
	{
		P.TreatAdminUnitDivisor = 1; P.TreatByAdminUnit = 0;
	}

	if (!GetInputParameter2(dat, dat2, "Vaccination trigger incidence per cell", "%lf", (void*)&(P.VaccCellIncThresh), 1, 1, 0)) P.VaccCellIncThresh = 1000000000;
	if (P.VaccCaseScale > 1 && P.VaccCellIncThresh == 1)
	{
		P.VaccCellIncThresh = P.VaccCellIncThresh * P.VaccCaseScale;
	}
	if (!GetInputParameter2(dat, dat2, "Relative susceptibility of vaccinated individual", "%lf", (void*)&(P.VaccSuscDrop), 1, 1, 0)) P.VaccSuscDrop = 1;
	if (P.HCWVaccSuscDrop < 0)
	{
		P.HCWVaccSuscDrop = P.VaccSuscDrop;
	}
	if (!GetInputParameter2(dat, dat2, "Relative susceptibility of individual vaccinated after switch time", "%lf", (void*)&(P.VaccSuscDrop2), 1, 1, 0)) P.VaccSuscDrop2 = 1;
	if (!GetInputParameter2(dat, dat2, "Switch time at which vaccine efficacy increases", "%lf", (void*)&(P.VaccTimeEfficacySwitch), 1, 1, 0)) P.VaccTimeEfficacySwitch = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Decay rate of vaccine efficacy (per year)", "%lf", (void*)&(P.VaccEfficacyDecay), 1, 1, 0)) P.VaccEfficacyDecay = 0;
	P.VaccEfficacyDecay /= DAYS_PER_YEAR;
	if (!GetInputParameter2(dat, dat2, "Relative infectiousness of vaccinated individual", "%lf", (void*)&(P.VaccInfDrop), 1, 1, 0)) P.VaccInfDrop = 1;
	if (!GetInputParameter2(dat, dat2, "Proportion of symptomatic cases resulting in death prevented by vaccination", "%lf", (void*)&(P.VaccMortDrop), 1, 1, 0)) P.VaccMortDrop = 0;
	if (!GetInputParameter2(dat, dat2, "Proportion of symptomatic cases prevented by vaccination", "%lf", (void*)&(P.VaccSympDrop), 1, 1, 0)) P.VaccSympDrop = 0;
	if (!GetInputParameter2(dat, dat2, "Delay to vaccinate", "%lf", (void*)&(P.VaccDelayMean), 1, 1, 0)) P.VaccDelayMean = 0;
	if ((P.VaccDelayMean == 1) & (P.VaccDelayScale != 1))
	{
		P.VaccDelayMean = P.VaccDelayMean * P.VaccDelayScale; //added this to also us to scale delay to vaccination time from command line - ggilani 28/02/17
	}
	if (!GetInputParameter2(dat, dat2, "Delay from vaccination to full protection", "%lf", (void*)&(P.VaccTimeToEfficacy), 1, 1, 0)) P.VaccTimeToEfficacy = 0;
	if ((P.VaccTimeToEfficacy == 1) & (P.VaccEffTimeScale != 1))
	{
		P.VaccTimeToEfficacy = P.VaccTimeToEfficacy * P.VaccEffTimeScale; //added this to also us to scale time to vaccination protection from command line - ggilani 28/02/17
	}
	if (!GetInputParameter2(dat, dat2, "Years between rounds of vaccination", "%lf", (void*)&(P.VaccCampaignInterval), 1, 1, 0)) P.VaccCampaignInterval = 1e10;
	if (!GetInputParameter2(dat, dat2, "Base vaccine doses per day", "%i", (void*)&(P.BaseVaccDosePerDay), 1, 1, 0)) P.BaseVaccDosePerDay = -1;
	if (!GetInputParameter2(dat, dat2, "Max vaccine doses per day", "%i", (void*)&(P.MaxVaccDosePerDay), 1, 1, 0)) P.MaxVaccDosePerDay = -1;
	if (!GetInputParameter2(dat, dat2, "Base geo vaccine doses per day", "%i", (void*)&(P.BaseVaccGeoDosePerDay), 1, 1, 0)) P.BaseVaccGeoDosePerDay = -1;
	if (!GetInputParameter2(dat, dat2, "Max geo vaccine doses per day", "%i", (void*)&(P.MaxVaccGeoDosePerDay), 1, 1, 0)) P.MaxVaccGeoDosePerDay = -1;
	if (!GetInputParameter2(dat, dat2, "Reset vaccination queue each day", "%i", (void*)&(P.ResetVaccQueue), 1, 1, 0)) P.ResetVaccQueue = 0;
	P.VaccCampaignInterval *= DAYS_PER_YEAR;
	if (!GetInputParameter2(dat, dat2, "Maximum number of rounds of vaccination", "%i", (void*)&(P.VaccMaxRounds), 1, 1, 0)) P.VaccMaxRounds = 1;
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(dat, dat2, "Proportion of households of cases vaccinated", "%lf", (void*)&(P.VaccPropCaseHouseholds), 1, 1, 0)) P.VaccPropCaseHouseholds = 0;
		if (!GetInputParameter2(dat, dat2, "Duration of household vaccination policy", "%lf", (void*)&(P.VaccHouseholdsDuration), 1, 1, 0)) P.VaccHouseholdsDuration = USHRT_MAX / P.TimeStepsPerDay;
	}
	if (P.DoHouseholds && P.DoPlaces)
	{
		if (!GetInputParameter2(dat, dat2, "Do ring vaccination", "%i", (void*)&(P.DoRingVaccination), 1, 1, 0)) P.DoRingVaccination = 0;
		if (P.DoRingVaccination)
		{
			//if(!GetInputParameter2(dat,dat2,"Number of rings to vaccinate","%i",(void *) &(P.NVaccRings),1,1,0)) P.NVaccRings=1;
			//if(P.VaccRingScale>1&&P.NVaccRings==1)
			//{
			//	P.NVaccRings=P.NVaccRings*P.VaccRingScale;
			//}

			if (!GetInputParameter2(dat, dat2, "Ring vaccination trigger incidence per cell", "%lf", (void*)&(P.RingVaccCellIncThresh), 1, 1, 0)) P.RingVaccCellIncThresh = 1000000000;
			if (!GetInputParameter2(dat, dat2, "Ring vaccination start time", "%lf", (void*)&(P.RingVaccTimeStartBase), 1, 1, 0)) P.RingVaccTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;

			if (!GetInputParameter2(dat, dat2, "Number of proportion of ring to vaccinate", "%i", (void*)&(P.NPropRingVacc), 1, 1, 0)) P.NPropRingVacc = 0;
			if (P.NPropRingVacc > 0)
			{
				if (P.NPropRingVacc >= MAX_CHANGE_POINTS) ERR_CRITICAL("MAX_CHANGE_POINTS too small\n");
				GetInputParameter(dat, dat2, "Change points proportion of ring to vaccinate", "%lf", (void*)P.ChangePointPropRingVacc, P.NPropRingVacc, 1, 0);
				GetInputParameter(dat, dat2, "List of proportions of ring to vaccinate", "%lf", (void*)P.ListPropRingVacc, P.NPropRingVacc, 1, 0);
			}
			P.CurrIndPropRingVacc = 0;

			//add time to increase to third ring
			if (!GetInputParameter2(dat, dat2, "Time to increase vaccination rings", "%lf", (void*)&(P.TimeToIncVaccRing), 1, 1, 0)) P.TimeToIncVaccRing = 1e10;
			//add time to efficiency for the third ring
			if (!GetInputParameter2(dat, dat2, "Delay from vaccination to full protection third ring", "%lf", (void*)&(P.VaccTimeToEfficacyThirdVaccRing), 1, 1, 0)) P.VaccTimeToEfficacyThirdVaccRing = 0;

			if (!GetInputParameter2(dat, dat2, "Proportion of ring to vaccinate", "%lf", (void*)&(P.PropRingVacc), 1, 1, 0)) P.PropRingVacc = 1;
			if (P.PropRingVacc == 1 && P.VaccPropScale != 1)
			{
				P.PropRingVacc = P.PropRingVacc * P.VaccPropScale;
			}
			if (!GetInputParameter2(dat, dat2, "Minimum vaccination age", "%i", (void*)&(P.MinVaccAge), 1, 1, 0)) P.MinVaccAge = 0;
		}
	}
	if (!GetInputParameter2(dat, dat2, "Days with no cases after which capacity is removed", "%lf", (void*)&(P.DaysToRemoveCapacity), 1, 1, 0)) P.DaysToRemoveCapacity = USHRT_MAX / P.TimeStepsPerDay;

	if (!GetInputParameter2(dat, dat2, "Vaccination start time", "%lf", (void*)&(P.VaccTimeStartBase), 1, 1, 0)) P.VaccTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Proportion of population vaccinated", "%lf", (void*)&(P.VaccProp), 1, 1, 0)) P.VaccProp = 0;
	if (!GetInputParameter2(dat, dat2, "Time taken to reach max vaccination coverage (in years)", "%lf", (void*)&(P.VaccCoverageIncreasePeriod), 1, 1, 0)) P.VaccCoverageIncreasePeriod = 0;
	P.VaccCoverageIncreasePeriod *= DAYS_PER_YEAR;
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
	if (!GetInputParameter2(dat, dat2, "Minimum radius from case to vaccinate", "%lf", (void*)&(P.VaccMinRadius), 1, 1, 0)) P.VaccMinRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Maximum number of vaccine courses available", "%lf", (void*)&(P.VaccMaxCoursesBase), 1, 1, 0)) P.VaccMaxCoursesBase = 1e20;
	if (!GetInputParameter2(dat, dat2, "Start time of additional vaccine production", "%lf", (void*)&(P.VaccNewCoursesStartTimeBase), 1, 1, 0)) P.VaccNewCoursesStartTimeBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Start time of boosted vaccine production", "%lf", (void*)&(P.VaccNewCoursesBoostStartTimeBase), 1, 1, 0)) P.VaccNewCoursesBoostStartTimeBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "End time of additional vaccine production", "%lf", (void*)&(P.VaccNewCoursesEndTime), 1, 1, 0)) P.VaccNewCoursesEndTime = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Do daily vaccine replenishment", "%i", (void*)&(P.DoVaccDailyReplenishment), 1, 1, 0)) P.DoVaccDailyReplenishment = 0;
	if (P.DoVaccDailyReplenishment)
	{
		if (!GetInputParameter2(dat, dat2, "Rate of additional vaccine production (courses per day)", "%lf", (void*)&(P.VaccNewCoursesRate), 1, 1, 0)) P.VaccNewCoursesRate = 0;

	}
	if (!GetInputParameter2(dat, dat2, "Do bulk vaccine replenishment", "%i", (void*)&(P.DoVaccBulkReplenishment), 1, 1, 0)) P.DoVaccBulkReplenishment = 0;
	if (P.DoVaccBulkReplenishment)
	{
		if (!GetInputParameter2(dat, dat2, "Init fraction of vaccine doses per delivery", "%lf", (void*)&(P.VaccNewCoursesInitFracBulk), 1, 1, 0)) P.VaccNewCoursesInitFracBulk = 0;
		P.VaccNewCoursesInitBulk = P.VaccNewCoursesInitFracBulk * P.VaccMaxCoursesBase;
		if (!GetInputParameter2(dat, dat2, "Init time between vaccine doses deliveries (in days)", "%lf", (void*)&(P.VaccNewCoursesInitDelay), 1, 1, 0)) P.VaccNewCoursesInitDelay = 0;
		if (!GetInputParameter2(dat, dat2, "Boosted number of vaccine doses per delivery", "%lf", (void*)&(P.VaccNewCoursesBulk), 1, 1, 0)) P.VaccNewCoursesBulk = 0;
		if (!GetInputParameter2(dat, dat2, "Boosted time between vaccine doses deliveries (in days)", "%lf", (void*)&(P.VaccNewCoursesDelay), 1, 1, 0)) P.VaccNewCoursesDelay = 0;
	}
	if (!GetInputParameter2(dat, dat2, "Apply mass rather than reactive vaccination", "%i", (void*)&(P.DoMassVacc), 1, 1, 0)) P.DoMassVacc = 0;
	if (!GetInputParameter2(dat, dat2, "Priority age range for mass vaccination", "%i", (void*)P.VaccPriorityGroupAge, 2, 1, 0)) { P.VaccPriorityGroupAge[0] = 1; P.VaccPriorityGroupAge[1] = 0; }
	if (P.DoAdUnits)
	{
		if (!GetInputParameter2(dat, dat2, "Vaccinate administrative units rather than rings", "%i", (void*)&(P.VaccByAdminUnit), 1, 1, 0)) P.VaccByAdminUnit = 0;
		if (!GetInputParameter2(dat, dat2, "Administrative unit divisor for vaccination", "%i", (void*)&(P.VaccAdminUnitDivisor), 1, 1, 0)) P.VaccAdminUnitDivisor = 1;
		if ((P.VaccAdminUnitDivisor == 0) || (P.VaccByAdminUnit == 0)) P.VaccAdminUnitDivisor = 1;
	}
	else
	{
		P.VaccAdminUnitDivisor = 1; P.VaccByAdminUnit = 0;
	}


	if (!GetInputParameter2(dat, dat2, "Movement restrictions trigger incidence per cell", "%i", (void*)&(P.MoveRestrCellIncThresh), 1, 1, 0)) P.MoveRestrCellIncThresh = 1000000000;
	if (!GetInputParameter2(dat, dat2, "Delay to start movement restrictions", "%lf", (void*)&(P.MoveDelayMean), 1, 1, 0)) P.MoveDelayMean = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of movement restrictions", "%lf", (void*)&(P.MoveRestrDuration), 1, 1, 0)) P.MoveRestrDuration = 7;
	if (!GetInputParameter2(dat, dat2, "Residual movements after restrictions", "%lf", (void*)&(P.MoveRestrEffect), 1, 1, 0)) P.MoveRestrEffect = 0;
	if (!GetInputParameter2(dat, dat2, "Minimum radius of movement restrictions", "%lf", (void*)&(P.MoveRestrRadius), 1, 1, 0)) P.MoveRestrRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Movement restrictions start time", "%lf", (void*)&(P.MoveRestrTimeStartBase), 1, 1, 0)) P.MoveRestrTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Impose blanket movement restrictions", "%i", (void*)&(P.DoBlanketMoveRestr), 1, 1, 0)) P.DoBlanketMoveRestr = 0;
	if (P.DoAdUnits)
	{
		if (!GetInputParameter2(dat, dat2, "Movement restrictions in administrative units rather than rings", "%i", (void*)&(P.MoveRestrByAdminUnit), 1, 1, 0)) P.MoveRestrByAdminUnit = 0;
		if (!GetInputParameter2(dat, dat2, "Administrative unit divisor for movement restrictions", "%i", (void*)&(P.MoveRestrAdminUnitDivisor), 1, 1, 0)) P.MoveRestrAdminUnitDivisor = 1;
		if ((P.MoveRestrAdminUnitDivisor == 0) || (P.MoveRestrByAdminUnit == 0)) P.MoveRestrAdminUnitDivisor = 1;
	}
	else
	{
		P.MoveRestrAdminUnitDivisor = 1; P.MoveRestrByAdminUnit = 0;
	}

	if (!GetInputParameter2(dat, dat2, "Include delay to case detection", "%i", (void*)&(P.DoDetectDelay), 1, 1, 0)) P.DoDetectDelay = 0;
	if (!GetInputParameter2(dat, dat2, "Do clustered case detection by household", "%i", (void*)&(P.DoClusterCaseDetection), 1, 1, 0)) P.DoClusterCaseDetection = 0;
	if (P.DoDetectDelay)
	{
		if (!GetInputParameter2(dat, dat2, "Mean detection delay", "%lf", (void*)&(P.DetectTime), 1, 1, 0)) P.DetectTime = 0;
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
	if (!GetInputParameter2(dat, dat2, "Include contact tracing", "%i", (void*)&(P.DoContactTracing), 1, 1, 0)) P.DoContactTracing = 0;
	if (P.DoContactTracing)
	{
		if (!GetInputParameter2(dat, dat2, "Contact tracing trigger incidence per cell", "%lf", (void*)&(P.ContactTracingCellIncThresh), 1, 1, 0)) P.ContactTracingCellIncThresh = 1000000000;
		if (!GetInputParameter2(dat, dat2, "Contact tracing start time", "%lf", (void*)&(P.ContactTracingTimeStartBase), 1, 1, 0)) P.ContactTracingTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
		GetInputParameter(dat, dat2, "Relative infectiousness of contact traced case", "%lf", (void*)&(P.RelativeInfectiousnessContactTraced), 1, 1, 0);
		GetInputParameter(dat, dat2, "Duration of contact tracing", "%lf", (void*)&(P.contactTraceDuration), 1, 1, 0);
		//GetInputParameter(dat,dat2,"Cumulative cases before contact tracing begins","%i",(void *) &(P.contactTraceCaseThreshold),1,1,0);
		//if(!GetInputParameter2(dat,dat2,"Cumulative cases when contact tracing increases","%i",(void *) &(P.contactTraceCaseThresholdInc),1,1,0)) P.contactTraceCaseThresholdInc=1e9;
		//couple of flags to help with contact tracing
		P.GlobalContractTracingStarted = P.GlobalContractTracingIncreased = 0;
		if (P.DoAdUnits)
		{
			//int AdunitCTCapacity[MAX_ADUNITS];
			//int AdunitCTCapacityInc[MAX_ADUNITS];
			//int AdunitCTThreshold[MAX_ADUNITS];
			if (!GetInputParameter(dat, dat2, "Contact tracing capacity per admin unit", "%i", (void*)&(P.AdunitCTCapacity), 1, 1, 0)) P.AdunitCTCapacity = 0;
			if (!GetInputParameter(dat, dat2, "Contact tracing threshold per admin unit", "%i", (void*)&(P.AdunitCTThreshold), 1, 1, 0)) P.AdunitCTThreshold = 0;
			if (!GetInputParameter2(dat, dat2, "Contact tracing increased capacity per admin unit", "%i", (void*)&(P.AdunitCTCapacityInc), 1, 1, 0)) P.AdunitCTCapacityInc = 0;
			if (!GetInputParameter2(dat, dat2, "Contact tracing capacity reached before increasing teams", "%lf", (void*)&(P.CapacityToMoreCT), 1, 1, 0)) P.CapacityToMoreCT = 1;
			if (!GetInputParameter2(dat, dat2, "Subsequent time to increase teams", "%lf", (void*)&(P.DelayToCT), 1, 1, 0)) P.DelayToCT = 0;

			for (i = 0; i < P.NumAdunits; i++)
			{
				//Some terrible hard coding to assign different thresholds to Guinea and Liberia&Sierra Leone!! Replace as soon as possible!! ggilani: 19/11/14
				if ((int)(AdUnits[i].id / P.CountryDivisor) == P.TargetCountry) //i.e. if the admin unit is in Guinea, assign contact tracing threshold 1
				{
					AdUnits[i].contactTraceCaseThreshold = P.CT_scale1 * P.AdunitCTThreshold;
					AdUnits[i].contactTraceCapacity = P.CT_scale1 * P.AdunitCTCapacity; //scaling up contact tracing capacity if needed
					AdUnits[i].contactTraceCapacityInc = P.CTinc_scale1 * P.AdunitCTCapacityInc; //scaling up increased contact tracing capacity if needed
					AdUnits[i].nextTimeToCT = 0;
				}
				else //else assign contact tracing threshold 2
				{
					AdUnits[i].contactTraceCaseThreshold = P.CT_thresh2 * P.AdunitCTThreshold;
					AdUnits[i].contactTraceCapacity = P.CT_scale2 * P.AdunitCTCapacity; //scaling up contact tracing capacity if needed
					AdUnits[i].contactTraceCapacityInc = P.CTinc_scale2 * P.AdunitCTCapacityInc;
				}
				//also set the flag for whether contact tracing has begun to zero, and set day on which contact tracing begins to zero for each admin unit
				//AdUnits[i].contactTraceThresholdCrossed=0;
				//AdUnits[i].contactTraceStartDay=1e9;
				//AdUnits[i].contactTraceCurrent=0; //added this to keep track of total number of people being contact traced at any given time - ggilani 07/06/17
			}
		}
		else
		{
			if (!GetInputParameter2(dat, dat2, "Countrywide capacity for contact tracing", "%i", (void*)&(P.contactTraceCapacity), 1, 1, 0)) P.contactTraceCapacity = 0;
			if (!GetInputParameter(dat, dat2, "Cumulative cases per country before contact tracing begins", "%i", (void*)&(P.contactTraceCaseThreshold), 1, 1, 0)) P.contactTraceCaseThreshold = 1;
		}
		if (P.DoContactTracing) //  check this!
		{
			if (!GetInputParameter2(dat, dat2, "Proportion of contacts to trace", "%lf", (void*)&(P.propContactTraced), 1, 1, 0)) P.propContactTraced = 1; //so if we don't specify this, everyone will be contact traced
			if (!GetInputParameter2(dat, dat2, "Proportion of contacts lost to follow up", "%lf", (void*)&(P.propContactLost), 1, 1, 0)) P.propContactLost = 0;
			if (!GetInputParameter2(dat, dat2, "Time to hospitalisation for contact traced case", "%lf", (void*)&(P.HospitalisationTime_contactTrace), 1, 1, 0)) P.HospitalisationTime_contactTrace = 1;
			if (!GetInputParameter2(dat, dat2, "Number of times to hospitalisation contact traced", "%i", (void*)&(P.NMeanTimeToHospCT), 1, 1, 0)) P.NMeanTimeToHosp = 0;
			if (P.NMeanTimeToHospCT > 0)
			{
				if (P.NMeanTimeToHospCT >= MAX_CHANGE_POINTS) ERR_CRITICAL("MAX_CHANGE_POINTS too small\n");
				GetInputParameter(dat, dat2, "Change points times to hospitalisation contact traced", "%lf", (void*)P.ChangePointMeanTimeToHospCT, P.NMeanTimeToHospCT, 1, 0);
				GetInputParameter(dat, dat2, "Times to hospitalisation contact traced", "%lf", (void*)P.MeanTimeToHospCT, P.NMeanTimeToHospCT, 1, 0);
			}
			P.CurrIndMeanTimeToHospCT = 0;
		}
	}
	//Moved the number of rings to here as it should be included for both ring vaccination and new contact tracing: ggilani 06/06/17
	if (P.DoRingVaccination || P.DoContactTracing)
	{
		if (!GetInputParameter2(dat, dat2, "Number of rings to vaccinate", "%i", (void*)&(P.NVaccRings), 1, 1, 0)) P.NVaccRings = 1;
		if (P.VaccRingScale > 1 && P.NVaccRings == 1)
		{
			P.NVaccRings = P.NVaccRings * P.VaccRingScale;
		}
		if (!GetInputParameter2(dat, dat2, "Probability of establishing vaccination ring/contact tracing", "%lf", (void*)&(P.ProbEstablishRing), 1, 1, 0)) P.ProbEstablishRing = 1;
	}

	//Capital City effect
	if (!GetInputParameter2(dat, dat2, "Include capital city effect", "%i", (void*)&(P.DoCapitalCityEffect), 1, 1, 0)) P.DoCapitalCityEffect = 0;
	if ((P.DoCapitalCityEffect) && (P.DoAdUnits)) //requires admin units
	{
		//read in target capital cities
		if (!GetInputParameter2(dat, dat2, "Capital city admin unit", "%i", (void*)&(P.CapitalCityAdunit), 1, 1, 0)) P.CapitalCityAdunit = 0;
		if (!GetInputParameter2(dat, dat2, "Capital city admin unit 2", "%i", (void*)&(P.CapitalCityAdunit2), 1, 1, 0)) P.CapitalCityAdunit2 = 0;
		if (!GetInputParameter2(dat, dat2, "Capital city admin unit 3", "%i", (void*)&(P.CapitalCityAdunit3), 1, 1, 0)) P.CapitalCityAdunit3 = 0;
		if (!GetInputParameter2(dat, dat2, "Do capital city distance effect", "%i", (void*)&(P.DoCapitalCityDistanceEffect), 1, 1, 0)) P.DoCapitalCityDistanceEffect = 0;
		if (!GetInputParameter2(dat, dat2, "Capital city distance effect", "%lf", (void*)&(P.CapitalCityDistanceEffect), 1, 1, 0)) P.CapitalCityDistanceEffect = 1.0;
		if (!GetInputParameter2(dat, dat2, "Do capital city population effect", "%i", (void*)&(P.DoCapitalCityPopEffect), 1, 1, 0)) P.DoCapitalCityPopEffect = 0;
		if (!GetInputParameter2(dat, dat2, "Capital city population effect", "%lf", (void*)&(P.CapitalCityPopEffect), 1, 1, 0)) P.CapitalCityPopEffect = 1.0;
		if (!GetInputParameter2(dat, dat2, "Do capital city additive effect", "%i", (void*)&(P.DoCapitalCityAddEffect), 1, 1, 0)) P.DoCapitalCityAddEffect = 0;
		//if(!GetInputParameter2(dat,dat2,"Capital city additive effect","%lf",(void *) &(P.CapitalCityAddEffect),1,1,0)) P.CapitalCityAddEffect=0.0;

	}

	if (!GetInputParameter2(dat, dat2, "Trigger incidence per cell for place closure", "%i", (void*)&(P.PlaceCloseCellIncThresh), 1, 1, 0)) P.PlaceCloseCellIncThresh = 1000000000;
	if (!GetInputParameter2(dat, dat2, "Delay to start place closure", "%lf", (void*)&(P.PlaceCloseDelayMean), 1, 1, 0)) P.PlaceCloseDelayMean = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of place closure", "%lf", (void*)&(P.PlaceCloseDurationBase), 1, 1, 0)) P.PlaceCloseDurationBase = 7;
	if (!GetInputParameter2(dat, dat2, "Duration of second place closure", "%lf", (void*)&(P.PlaceCloseDuration2), 1, 1, 0)) P.PlaceCloseDuration2 = 7;
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(dat, dat2, "Proportion of places remaining open after closure by place type", "%lf", (void*)P.PlaceCloseEffect, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.PlaceCloseEffect[i] = 1;
	}
	if (P.DoHouseholds)
		if (!GetInputParameter2(dat, dat2, "Relative houshold contact rate after closure", "%lf", (void*)&P.PlaceCloseHouseholdRelContact, 1, 1, 0)) P.PlaceCloseHouseholdRelContact = 1;
	if (!GetInputParameter2(dat, dat2, "Relative spatial contact rate after closure", "%lf", (void*)&P.PlaceCloseSpatialRelContact, 1, 1, 0))
	{
		P.PlaceCloseSpatialRelContact = 1;
	}
	if (!GetInputParameter2(dat, dat2, "Minimum radius for place closure", "%lf", (void*)&(P.PlaceCloseRadius), 1, 1, 0)) P.PlaceCloseRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Place closure start time", "%lf", (void*)&(P.PlaceCloseTimeStartBase), 1, 1, 0)) P.PlaceCloseTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Place closure second start time", "%lf", (void*)&(P.PlaceCloseTimeStartBase2), 1, 1, 0)) P.PlaceCloseTimeStartBase2 = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Places close only once", "%i", (void*)&(P.DoPlaceCloseOnceOnly), 1, 1, 0)) P.DoPlaceCloseOnceOnly = 0;
	if (!GetInputParameter2(dat, dat2, "Place closure incidence threshold", "%i", (void*)&(P.PlaceCloseIncTrig), 1, 1, 0)) P.PlaceCloseIncTrig = 1;
	if (!GetInputParameter2(dat, dat2, "Place closure fractional incidence threshold", "%lf", (void*)&(P.PlaceCloseFracIncTrig), 1, 1, 0)) P.PlaceCloseFracIncTrig = 0;
	if ((P.DoAdUnits) && (P.DoPlaces))
	{
		if (!GetInputParameter2(dat, dat2, "Place closure in administrative units rather than rings", "%i", (void*)&(P.PlaceCloseByAdminUnit), 1, 1, 0)) P.PlaceCloseByAdminUnit = 0;
		if (!GetInputParameter2(dat, dat2, "Administrative unit divisor for place closure", "%i", (void*)&(P.PlaceCloseAdminUnitDivisor), 1, 1, 0)) P.PlaceCloseAdminUnitDivisor = 1;
		if (!GetInputParameter2(dat, dat2, "Place types to close for admin unit closure (0/1 array)", "%i", (void*)&(P.PlaceCloseAdunitPlaceTypes), P.PlaceTypeNum, 1, 0))
			for (i = 0; i < P.PlaceTypeNum; i++) P.PlaceCloseAdunitPlaceTypes[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Cumulative proportion of place members needing to become sick for admin unit closure", "%lf", (void*)&(P.PlaceCloseCasePropThresh), 1, 1, 0)) P.PlaceCloseCasePropThresh = 2;
		if (!GetInputParameter2(dat, dat2, "Proportion of places in admin unit needing to pass threshold for place closure", "%lf", (void*)&(P.PlaceCloseAdunitPropThresh), 1, 1, 0)) P.PlaceCloseAdunitPropThresh = 2;
		if ((P.PlaceCloseAdminUnitDivisor < 1) || (P.PlaceCloseByAdminUnit == 0)) P.PlaceCloseAdminUnitDivisor = 1;
	}
	else
	{
		P.PlaceCloseAdminUnitDivisor = 1; P.PlaceCloseByAdminUnit = 0;
	}

	if (!GetInputParameter2(dat, dat2, "Trigger incidence per cell for social distancing", "%i", (void*)&(P.SocDistCellIncThresh), 1, 1, 0)) P.SocDistCellIncThresh = 1000000000;
	if (!GetInputParameter2(dat, dat2, "Duration of social distancing", "%lf", (void*)&(P.SocDistDuration), 1, 1, 0)) P.SocDistDuration = 7;
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(dat, dat2, "Relative place contact rate given social distancing by place type", "%lf", (void*)P.SocDistPlaceEffect, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.SocDistPlaceEffect[i] = 1;
		if (!GetInputParameter2(dat, dat2, "Relative place contact rate given enhanced social distancing by place type", "%lf", (void*)P.ESocDistPlaceEffect, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.ESocDistPlaceEffect[i] = 1;
	}
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(dat, dat2, "Relative houshold contact rate given social distancing", "%lf", (void*)&P.SocDistHouseholdEffect, 1, 1, 0)) P.SocDistHouseholdEffect = 1;
		if (!GetInputParameter2(dat, dat2, "Relative houshold contact rate given enhanced social distancing", "%lf", (void*)&P.ESocDistHouseholdEffect, 1, 1, 0)) P.ESocDistHouseholdEffect = 1;
	}
	if (!GetInputParameter2(dat, dat2, "Relative spatial contact rate given social distancing", "%lf", (void*)&P.SocDistSpatialEffect, 1, 1, 0)) P.SocDistSpatialEffect = 1;
	if (!GetInputParameter2(dat, dat2, "Minimum radius for social distancing", "%lf", (void*)&(P.SocDistRadius), 1, 1, 0)) P.SocDistRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Social distancing start time", "%lf", (void*)&(P.SocDistTimeStartBase), 1, 1, 0)) P.SocDistTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Proportion compliant with enhanced social distancing", "%lf", (void*)&(P.ESocProportionCompliant), 1, 1, 0)) P.ESocProportionCompliant = 0;
	if (!GetInputParameter2(dat, dat2, "Relative spatial contact rate given enhanced social distancing", "%lf", (void*)&P.ESocDistSpatialEffect, 1, 1, 0)) P.ESocDistSpatialEffect = 1;


	if (!GetInputParameter2(dat, dat2, "Airport closure effectiveness", "%lf", (void*)&(P.AirportCloseEffectiveness), 1, 1, 0)) P.AirportCloseEffectiveness = 0;
	P.AirportCloseEffectiveness = 1.0 - P.AirportCloseEffectiveness;
	if (!GetInputParameter2(dat, dat2, "Airport closure start time", "%lf", (void*)&(P.AirportCloseTimeStartBase), 1, 1, 0)) P.AirportCloseTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Airport closure duration", "%lf", (void*)&(P.AirportCloseDuration), 1, 1, 0)) P.AirportCloseDuration = USHRT_MAX / P.TimeStepsPerDay;

	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(dat, dat2, "Household quarantine start time", "%lf", (void*)&(P.HQuarantineTimeStartBase), 1, 1, 0)) P.HQuarantineTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
		if (!GetInputParameter2(dat, dat2, "Delay to start household quarantine", "%lf", (void*)&(P.HQuarantineHouseDelay), 1, 1, 0)) P.HQuarantineHouseDelay = 0;
		if (!GetInputParameter2(dat, dat2, "Length of time households are quarantined", "%lf", (void*)&(P.HQuarantineHouseDuration), 1, 1, 0)) P.HQuarantineHouseDuration = 0;
		if (!GetInputParameter2(dat, dat2, "Duration of household quarantine policy", "%lf", (void*)&(P.HQuarantinePolicyDuration), 1, 1, 0)) P.HQuarantinePolicyDuration = USHRT_MAX / P.TimeStepsPerDay;
		if (!GetInputParameter2(dat, dat2, "Relative household contact rate after quarantine", "%lf", (void*)&(P.HQuarantineHouseEffect), 1, 1, 0)) P.HQuarantineHouseEffect = 1;
		if (P.DoPlaces)
		{
			if (!GetInputParameter2(dat, dat2, "Residual place contacts after household quarantine by place type", "%lf", (void*)P.HQuarantinePlaceEffect, P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++) P.HQuarantinePlaceEffect[i] = 1;
		}
		if (!GetInputParameter2(dat, dat2, "Residual spatial contacts after household quarantine", "%lf", (void*)&(P.HQuarantineSpatialEffect), 1, 1, 0)) P.HQuarantineSpatialEffect = 1;
		if (!GetInputParameter2(dat, dat2, "Household level compliance with quarantine", "%lf", (void*)&(P.HQuarantinePropHouseCompliant), 1, 1, 0)) P.HQuarantinePropHouseCompliant = 1;
		if (!GetInputParameter2(dat, dat2, "Individual level compliance with quarantine", "%lf", (void*)&(P.HQuarantinePropIndivCompliant), 1, 1, 0)) P.HQuarantinePropIndivCompliant = 1;
	}
	else
		P.HQuarantineTimeStartBase = 1e10;
	if (!GetInputParameter2(dat, dat2, "Case isolation start time", "%lf", (void*)&(P.CaseIsolationTimeStartBase), 1, 1, 0)) P.CaseIsolationTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Proportion of detected cases isolated", "%lf", (void*)&(P.CaseIsolationProp), 1, 1, 0)) P.CaseIsolationProp = 0;
	if (!GetInputParameter2(dat, dat2, "Delay to start case isolation", "%lf", (void*)&(P.CaseIsolationDelay), 1, 1, 0)) P.CaseIsolationDelay = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of case isolation", "%lf", (void*)&(P.CaseIsolationDuration), 1, 1, 0)) P.CaseIsolationDuration = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of case isolation policy", "%lf", (void*)&(P.CaseIsolationPolicyDuration), 1, 1, 0)) P.CaseIsolationPolicyDuration = 1e10;
	if (!GetInputParameter2(dat, dat2, "Residual contacts after case isolation", "%lf", (void*)&(P.CaseIsolationEffectiveness), 1, 1, 0)) P.CaseIsolationEffectiveness = 1;
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(dat, dat2, "Residual household contacts after case isolation", "%lf", (void*)&(P.CaseIsolationHouseEffectiveness), 1, 1, 0))
			P.CaseIsolationHouseEffectiveness = P.CaseIsolationEffectiveness;
	}
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(dat, dat2, "Number of key workers randomly distributed in the population", "%i", (void*)&(P.KeyWorkerPopNum), 1, 1, 0)) P.KeyWorkerPopNum = 0;
		if (!GetInputParameter2(dat, dat2, "Number of key workers in different places by place type", "%i", (void*)P.KeyWorkerPlaceNum, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.KeyWorkerPlaceNum[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Proportion of staff who are key workers per chosen place by place type", "%lf", (void*)P.KeyWorkerPropInKeyPlaces, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.KeyWorkerPropInKeyPlaces[i] = 1.0;
		if (!GetInputParameter2(dat, dat2, "Trigger incidence per cell for key worker prophylaxis", "%i", (void*)&(P.KeyWorkerProphCellIncThresh), 1, 1, 0)) P.KeyWorkerProphCellIncThresh = 1000000000;
		if (!GetInputParameter2(dat, dat2, "Key worker prophylaxis start time", "%lf", (void*)&(P.KeyWorkerProphTimeStartBase), 1, 1, 0)) P.KeyWorkerProphTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
		if (!GetInputParameter2(dat, dat2, "Duration of key worker prophylaxis", "%lf", (void*)&(P.KeyWorkerProphDuration), 1, 1, 0)) P.KeyWorkerProphDuration = 0;
		if (!GetInputParameter2(dat, dat2, "Time interval from start of key worker prophylaxis before policy restarted", "%lf", (void*)&(P.KeyWorkerProphRenewalDuration), 1, 1, 0)) P.KeyWorkerProphRenewalDuration = P.KeyWorkerProphDuration;
		if (P.DoHouseholds)
		{
			if (!GetInputParameter2(dat, dat2, "Proportion of key workers whose households are also treated as key workers", "%lf", (void*)&(P.KeyWorkerHouseProp), 1, 1, 0)) P.KeyWorkerHouseProp = 0;
		}
		if (!GetInputParameter2(dat, dat2, "Minimum radius for key worker prophylaxis", "%lf", (void*)&(P.KeyWorkerProphRadius), 1, 1, 0)) P.KeyWorkerProphRadius = 0;
	}
	else
	{
		P.KeyWorkerPopNum = 0;
		P.KeyWorkerProphTimeStartBase = 1e10;
	}
	if (!GetInputParameter2(dat, dat2, "Initial rate of importation of infections", "%lf", (void*)&(P.InfectionImportRate1), 1, 1, 0)) P.InfectionImportRate1 = 0;
	if (!GetInputParameter2(dat, dat2, "Changed rate of importation of infections", "%lf", (void*)&(P.InfectionImportRate2), 1, 1, 0)) P.InfectionImportRate2 = 0;
	if (!GetInputParameter2(dat, dat2, "Time when infection rate changes", "%lf", (void*)&(P.InfectionImportChangeTime), 1, 1, 0)) P.InfectionImportChangeTime = 1e10;
	if (!GetInputParameter2(dat, dat2, "Imports via air travel", "%i", (void*)&(P.DoImportsViaAirports), 1, 1, 0)) P.DoImportsViaAirports = 0;
	if (!GetInputParameter2(dat, dat2, "Length of importation time profile provided", "%i", (void*)&(P.DurImportTimeProfile), 1, 1, 0)) P.DurImportTimeProfile = 0;
	if (P.DurImportTimeProfile > 0)
	{
		if (P.DurImportTimeProfile >= MAX_DUR_IMPORT_PROFILE) ERR_CRITICAL("MAX_DUR_IMPORT_PROFILE too small\n");
		GetInputParameter(dat, dat2, "Daily importation time profile", "%lf", (void*)P.ImportInfectionTimeProfile, P.DurImportTimeProfile, 1, 0);
	}
	//added some more code to allow us to specify where to import to
	if (!GetInputParameter2(dat, dat2, "Import to specific location", "%i", (void*)&(P.DoImportToSpecLocation), 1, 1, 0)) P.DoImportToSpecLocation = 0;
	if (!GetInputParameter2(dat, dat2, "Location for imported cases", "%lf", (void*)&(P.ImportLocation[0]), 2, 1, 0))
	{
		P.ImportLocation[0] = -1000;
		P.ImportLocation[1] = -1000;
	}
	//Added this to parameter list so that recording infection events (and the number to record) can easily be turned off and on: ggilani - 10/10/2014
	if (!GetInputParameter2(dat, dat2, "Record infection events", "%i", (void*)&(P.DoRecordInfEvents), 1, 1, 0)) P.DoRecordInfEvents = 0;
	if (P.DoRecordInfEvents)
	{
		if (!GetInputParameter2(dat, dat2, "Max number of infection events to record", "%i", (void*)&(P.MaxInfEvents), 1, 1, 0)) P.MaxInfEvents = 1000;
		if (!GetInputParameter2(dat, dat2, "Record infection events per run", "%i", (void*)&(P.RecordInfEventsPerRun), 1, 1, 0)) P.RecordInfEventsPerRun = 0;
	}
	else
	{
		P.MaxInfEvents = 0;
	}
	//Include a limit to the number of infections to simulate, if this happens before time runs out
	if (!GetInputParameter2(dat, dat2, "Limit number of infections", "%i", (void*)&(P.LimitNumInfections), 1, 1, 0)) P.LimitNumInfections = 0;
	if (P.LimitNumInfections)
	{
		if (!GetInputParameter2(dat, dat2, "Max number of infections", "%i", (void*)&(P.MaxNumInfections), 1, 1, 0)) P.MaxNumInfections = 60000;
	}
	//Add cross border contact parameter
	if (!GetInputParameter2(dat, dat2, "Proportion of infections allowed across country borders", "%lf", (void*)&(P.PropCrossBorderInf), 1, 1, 0)) P.PropCrossBorderInf = 1;
	if (P.BC_scale != 1.0)
	{
		P.PropCrossBorderInf *= P.BC_scale; //scale cross border contact if necessary
	}
	//Add origin-destination matrix parameter
	if (!GetInputParameter2(dat, dat2, "Output origin destination matrix", "%i", (void*)&(P.DoOriginDestinationMatrix), 1, 1, 0)) P.DoOriginDestinationMatrix = 0;
	//Some parameters relating to road networks: ggilani 12/02/15
	if (P.DoRoadNetwork)
	{
		if (!GetInputParameter2(dat, dat2, "Do distance road effect", "%i", (void*)&(P.DoRoadDistanceEffect), 1, 1, 0)) P.DoRoadDistanceEffect = 0;
		if (!GetInputParameter2(dat, dat2, "Distance scaling for road accessibility", "%lf", (void*)&(P.RoadAccessDistance), 1, 1, 0)) P.RoadAccessDistance = 1;
		if (!GetInputParameter2(dat, dat2, "Do population road effect", "%i", (void*)&(P.DoRoadPopEffect), 1, 1, 0)) P.DoRoadPopEffect = 0;
		if (!GetInputParameter2(dat, dat2, "Population scaling for road accessibility", "%lf", (void*)&(P.RoadAccessPop), 1, 1, 0)) P.RoadAccessPop = 1;
		if (!GetInputParameter2(dat, dat2, "Maximum road type to consider", "%i", (void*)&(P.MaxRoadType), 1, 1, 0)) P.MaxRoadType = 2; //2 corresponds to primary roads in the road network file
		if (!GetInputParameter2(dat, dat2, "Maximum cells neighbouring road", "%i", (void*)&(P.MaxRoadNeighbour), 1, 1, 0)) P.MaxRoadNeighbour = 2;
	}

	fclose(dat);

	if (P.DoOneGen != 0) P.DoOneGen = 1;
	P.ColourPeriod = 2000;
	P.MoveRestrRadius2 = P.MoveRestrRadius * P.MoveRestrRadius;
	P.SocDistRadius2 = P.SocDistRadius * P.SocDistRadius;
	P.VaccRadius2 = P.VaccRadius * P.VaccRadius;
	if (!P.LimitGeoVaccDosesPerCase)
	{
		P.VaccRadiusHighDensity2 = P.VaccRadiusHighDensity * P.VaccRadiusHighDensity;
	}
	P.VaccMinRadius2 = P.VaccMinRadius * P.VaccMinRadius;
	P.TreatRadius2 = P.TreatRadius * P.TreatRadius;
	P.PlaceCloseRadius2 = P.PlaceCloseRadius * P.PlaceCloseRadius;
	P.KeyWorkerProphRadius2 = P.KeyWorkerProphRadius * P.KeyWorkerProphRadius;
	if (P.TreatRadius2 == 0) P.TreatRadius2 = -1;
	if (P.VaccRadius2 == 0) P.VaccRadius2 = -1;
	if (P.PlaceCloseRadius2 == 0) P.PlaceCloseRadius2 = -1;
	if (P.MoveRestrRadius2 == 0) P.MoveRestrRadius2 = -1;
	if (P.SocDistRadius2 == 0) P.SocDistRadius2 = -1;
	if (P.KeyWorkerProphRadius2 == 0) P.KeyWorkerProphRadius2 = -1;
	if (P.TreatCellIncThresh < 1) P.TreatCellIncThresh = 1;
	if (P.MoveRestrCellIncThresh < 1) P.MoveRestrCellIncThresh = 1;
	if (P.PlaceCloseCellIncThresh < 1) P.PlaceCloseCellIncThresh = 1;
	if (P.KeyWorkerProphCellIncThresh < 1) P.KeyWorkerProphCellIncThresh = 1;
	P.usHQuarantineHouseDuration = ((unsigned short int) (P.HQuarantineHouseDuration * P.TimeStepsPerDay));
	P.usVaccTimeToEfficacy = ((unsigned short int) (P.VaccTimeToEfficacy * P.TimeStepsPerDay));
	P.usVaccTimeToEfficacyThirdRing = ((unsigned short int) (P.VaccTimeToEfficacyThirdVaccRing * P.TimeStepsPerDay));
	P.usVaccTimeEfficacySwitch = ((unsigned short int) (P.VaccTimeEfficacySwitch * P.TimeStepsPerDay));
	P.usCaseIsolationDelay = ((unsigned short int) (P.CaseIsolationDelay * P.TimeStepsPerDay));
	P.usCaseIsolationDuration = ((unsigned short int) (P.CaseIsolationDuration * P.TimeStepsPerDay));
	P.usCaseAbsenteeismDuration = ((unsigned short int) (P.CaseAbsenteeismDuration * P.TimeStepsPerDay));
	P.usCaseAbsenteeismDelay = ((unsigned short int) (P.CaseAbsenteeismDelay * P.TimeStepsPerDay));
	if (P.DoUTM_coords)
	{
		for (i = 0; i <= 1000; i++)
		{
			asin2sqx[i] = asin(sqrt(((double)(i)) / 1000));
			asin2sqx[i] = asin2sqx[i] * asin2sqx[i];
		}
		for (t = 0; t <= 360; t++)
		{
			sinx[(int)t] = sin(PI * t / 180);
			cosx[(int)t] = cos(PI * t / 180);
		}
	}
	fprintf(stderr, "Parameters read\n");
}

void ReadInterventions(char* IntFile)
{
	FILE* dat, * dat2;
	double r, s, t, dt, startt, stopt;
	int n, i, j, k, au, ni, f, nsr;
	char buf[65536], txt[65536];
	intervention CurInterv;

	fprintf(stderr, "Reading intervention file.\n");
	if (!(dat = fopen(IntFile, "r"))) ERR_CRITICAL("Unable to open intervention file\n");
	fscanf(dat, "%*[^<]"); // needs to be separate line because start of file
	fscanf(dat, "<%[^>]", txt);
	if (strcmp(txt, "\?xml version=\"1.0\" encoding=\"ISO-8859-1\"\?") != 0) ERR_CRITICAL("Intervention file not XML.\n");
	fscanf(dat, "%*[^<]<%[^>]", txt);
	if (strcmp(txt, "InterventionSettings") != 0) ERR_CRITICAL("Intervention has no top level.\n");
	ni = 0;
	while (!feof(dat))
	{
		fscanf(dat, "%*[^<]<%[^>]", txt);
		if (strcmp(txt, "intervention") == 0)
		{
			ni++;
			fscanf(dat, "%*[^<]<%[^>]", txt);
			if (strcmp(txt, "parameters") != 0) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			if (!GetXMLNode(dat, "Type", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			if (strcmp(txt, "Treatment") == 0)
				CurInterv.InterventionType = 0;
			else if (strcmp(txt, "Vaccination") == 0)
				CurInterv.InterventionType = 1;
			else if (strcmp(txt, "ITN") == 0)
				CurInterv.InterventionType = 2;
			else if (strcmp(txt, "IRS") == 0)
				CurInterv.InterventionType = 3;
			else if (strcmp(txt, "GM") == 0)
				CurInterv.InterventionType = 4;
			else if (strcmp(txt, "MSAT") == 0)
				CurInterv.InterventionType = 5;
			else
				sscanf(txt, "%i", &CurInterv.InterventionType);
			if (!GetXMLNode(dat, "AUThresh", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%i", &CurInterv.DoAUThresh);
			if (!GetXMLNode(dat, "StartTime", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lg", &CurInterv.StartTime);
			startt = CurInterv.StartTime;
			if (!GetXMLNode(dat, "StopTime", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lg", &CurInterv.StopTime);
			stopt = CurInterv.StopTime;
			if (!GetXMLNode(dat, "MinDuration", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lg", &CurInterv.MinDuration);
			CurInterv.MinDuration *= DAYS_PER_YEAR;
			if (!GetXMLNode(dat, "RepeatInterval", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lg", &CurInterv.RepeatInterval);
			CurInterv.RepeatInterval *= DAYS_PER_YEAR;
			if (!GetXMLNode(dat, "MaxPrevAtStart", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lg", &CurInterv.StartThresholdHigh);
			if (!GetXMLNode(dat, "MinPrevAtStart", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lg", &CurInterv.StartThresholdLow);
			if (!GetXMLNode(dat, "MaxPrevAtStop", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lg", &CurInterv.StopThreshold);
			if (GetXMLNode(dat, "NoStartAfterMinDur", "parameters", txt, 1))
				sscanf(txt, "%i", &CurInterv.NoStartAfterMin);
			else
				CurInterv.NoStartAfterMin = 0;
			if (!GetXMLNode(dat, "Level", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lg", &CurInterv.Level);
			if (GetXMLNode(dat, "LevelCellVar", "parameters", txt, 1))
				sscanf(txt, "%lg", &CurInterv.LevelCellVar);
			else
				CurInterv.LevelCellVar = 0;
			if (GetXMLNode(dat, "LevelAUVar", "parameters", txt, 1))
				sscanf(txt, "%lg", &CurInterv.LevelAUVar);
			else
				CurInterv.LevelCellVar = 0;
			if (GetXMLNode(dat, "LevelCountryVar", "parameters", txt, 1))
				sscanf(txt, "%lg", &CurInterv.LevelCountryVar);
			else
				CurInterv.LevelCellVar = 0;
			if (GetXMLNode(dat, "LevelClustering", "parameters", txt, 1))
				sscanf(txt, "%lg", &CurInterv.LevelClustering);
			else
				CurInterv.LevelClustering = 0;
			if (GetXMLNode(dat, "ControlParam", "parameters", txt, 1))
				sscanf(txt, "%lg", &CurInterv.ControlParam);
			else
				CurInterv.ControlParam = 0;
			if (GetXMLNode(dat, "TimeOffset", "parameters", txt, 1))
				sscanf(txt, "%lg", &CurInterv.TimeOffset);
			else
				CurInterv.TimeOffset = 0;

			if (!GetXMLNode(dat, "MaxRounds", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%u", &CurInterv.MaxRounds);
			if (!GetXMLNode(dat, "MaxResource", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%u", &CurInterv.MaxResource);
			if (GetXMLNode(dat, "NumSequentialReplicas", "parameters", txt, 1))
				sscanf(txt, "%i", &nsr);
			else
				nsr = 0;
			do { fscanf(dat, "%*[^<]<%[^>]", txt); } while ((strcmp(txt, "/intervention") != 0) && (strcmp(txt, "/parameters") != 0) && (!feof(dat)));
			if (strcmp(txt, "/parameters") != 0) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			fscanf(dat, "%*[^<]<%[^>]", txt);
			if ((strcmp(txt, "adunits") != 0) && (strcmp(txt, "countries") != 0)) ERR_CRITICAL("Incomplete adunits/countries specification in intervention file\n");
			if (strcmp(txt, "adunits") == 0)
			{
				while (GetXMLNode(dat, "A", "adunits", buf, 0))
				{
					sscanf(buf, "%s", txt);
					j = atoi(txt);
					if (j == 0)
					{
						f = 1; au = -1;
						do
						{
							au++; f = strcmp(txt, AdUnits[au].ad_name);
						} while ((f) && (au < P.NumAdunits));
						if (!f)
						{
							r = fabs(CurInterv.Level) + (2.0 * ranf() - 1) * CurInterv.LevelAUVar;
							if ((CurInterv.Level < 1) && (r > 1))
								r = 1;
							else if (r < 0)
								r = 0;
							for (k = 0; k <= nsr; k++)
							{
								AdUnits[au].InterventionList[AdUnits[au].NI] = CurInterv;
								AdUnits[au].InterventionList[AdUnits[au].NI].Level = r;
								AdUnits[au].InterventionList[AdUnits[au].NI].StartTime = startt + ((double)k) * (stopt - startt);
								AdUnits[au].InterventionList[AdUnits[au].NI].StopTime = stopt + ((double)k) * (stopt - startt);
								AdUnits[au].NI++;
							}
						}
					}
					else
					{
						k = (j % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor;
						au = P.AdunitLevel1Lookup[k];
						if ((au >= 0) && (AdUnits[au].id / P.AdunitLevel1Divisor == j / P.AdunitLevel1Divisor))
						{
							r = CurInterv.Level + (2.0 * ranf() - 1) * CurInterv.LevelAUVar;
							if ((CurInterv.Level < 1) && (r > 1))
								r = 1;
							else if (r < 0)
								r = 0;
							for (k = 0; k <= nsr; k++)
							{
								AdUnits[au].InterventionList[AdUnits[au].NI] = CurInterv;
								AdUnits[au].InterventionList[AdUnits[au].NI].Level = r;
								AdUnits[au].InterventionList[AdUnits[au].NI].StartTime = startt + ((double)k) * (stopt - startt);
								AdUnits[au].InterventionList[AdUnits[au].NI].StopTime = stopt + ((double)k) * (stopt - startt);
								AdUnits[au].NI++;
							}
						}
					}
				}
			}
			else
			{
				while (GetXMLNode(dat, "C", "countries", buf, 0))
				{
					s = (2.0 * ranf() - 1) * CurInterv.LevelCountryVar;
					sscanf(buf, "%s", txt);
					j = atoi(txt);
					for (au = 0; au < P.NumAdunits; au++)
						if (((j == 0) && (strcmp(txt, AdUnits[au].cnt_name) == 0)) || ((j > 0) && (j == AdUnits[au].cnt_id)))
						{
							r = CurInterv.Level + (2.0 * ranf() - 1) * CurInterv.LevelAUVar + s;
							if ((CurInterv.Level < 1) && (r > 1))
								r = 1;
							else if (r < 0)
								r = 0;
							for (k = 0; k <= nsr; k++)
							{
								AdUnits[au].InterventionList[AdUnits[au].NI] = CurInterv;
								AdUnits[au].InterventionList[AdUnits[au].NI].Level = r;
								AdUnits[au].InterventionList[AdUnits[au].NI].StartTime = startt + ((double)k) * (stopt - startt);
								AdUnits[au].InterventionList[AdUnits[au].NI].StopTime = stopt + ((double)k) * (stopt - startt);
								AdUnits[au].NI++;
							}
						}
				}
			}
			fscanf(dat, "%*[^<]<%[^>]", txt);
			if (strcmp(txt, "/intervention") != 0) ERR_CRITICAL("Incorrect intervention specification in intervention file\n");
		}
	}
	if (strcmp(txt, "/InterventionSettings") != 0) ERR_CRITICAL("Intervention has no top level closure.\n");
	fprintf(stderr, "%i interventions read\n", ni);
	fclose(dat);
}

int GetXMLNode(FILE* dat, char* NodeName, char* ParentName, char* Value, int ResetFilePos)
{
	// ResetFilePos=1 leaves dat cursor in same position as when function was called. 0 leaves it at end of NodeName closure
	// GetXMLNode returns 1 if NodeName found, 0 otherwise. If NodeName not found, ParentName closure must be

	char buf[65536], CloseNode[2048], CloseParent[2048];
	int CurPos, ret;

	sprintf(CloseParent, "/%s", ParentName);
	CurPos = ftell(dat);
	do
	{
		fscanf(dat, "%*[^<]<%[^>]", buf);
	} while ((strcmp(buf, CloseParent) != 0) && (strcmp(buf, NodeName) != 0) && (!feof(dat)));
	if (strcmp(buf, CloseParent) == 0)
		ret = 0;
	else
	{
		if (strcmp(buf, NodeName) != 0) ERR_CRITICAL("Incomplete node specification in XML file\n");
		fscanf(dat, ">%[^<]", buf);
		if (strlen(buf) < 2048) strcpy(Value, buf);
		//		fprintf(stderr,"# %s=%s\n",NodeName,Value);
		fscanf(dat, "<%[^>]", buf);
		sprintf(CloseNode, "/%s", NodeName);
		if (strcmp(buf, CloseNode) != 0) ERR_CRITICAL("Incomplete node specification in XML file\n");
		ret = 1;
	}
	if (ResetFilePos) fseek(dat, CurPos, 0);
	return ret;
}

void ReadAirTravel(char* AirTravelFile)
{
	int i, j, k, l;
	float sc, t, t2;
	float* buf;
	double traf;
	char outname[1024];
	FILE* dat;

	fprintf(stderr, "Reading airport data...\nAirports with no connections = ");
	if (!(dat = fopen(AirTravelFile, "r"))) ERR_CRITICAL("Unable to open airport file\n");
	fscanf(dat, "%i %i", &P.Nairports, &P.Air_popscale);
	sc = ((double)P.N) / ((double)P.Air_popscale);
	if (P.Nairports > MAX_AIRPORTS) ERR_CRITICAL("Too many airports\n");
	if (P.Nairports < 2) ERR_CRITICAL("Too few airports\n");
	if (!(buf = (float*)calloc(P.Nairports + 1, sizeof(float)))) ERR_CRITICAL("Unable to allocate airport storage\n");
	if (!(Airports = (airport*)calloc(P.Nairports, sizeof(airport)))) ERR_CRITICAL("Unable to allocate airport storage\n");
	for (i = 0; i < P.Nairports; i++)
	{
		fscanf(dat, "%f %f %lg", &(Airports[i].loc_x), &(Airports[i].loc_y), &traf);
		traf *= (P.AirportTrafficScale * sc);
		if ((Airports[i].loc_x < P.SpatialBoundingBox[0]) || (Airports[i].loc_x > P.SpatialBoundingBox[2])
			|| (Airports[i].loc_y < P.SpatialBoundingBox[1]) || (Airports[i].loc_y > P.SpatialBoundingBox[3]))
		{
			Airports[i].loc_x = Airports[i].loc_y = -1;
			Airports[i].total_traffic = 0;
		}
		else
		{
			//fprintf(stderr,"(%f,%f) ",Airports[i].loc_x,Airports[i].loc_y);
			Airports[i].loc_x -= P.SpatialBoundingBox[0];
			Airports[i].loc_y -= P.SpatialBoundingBox[1];
			Airports[i].total_traffic = traf;
		}
		t = 0;
		for (j = k = 0; j < P.Nairports; j++)
		{
			fscanf(dat, "%f", buf + j);
			if (buf[j] > 0) { k++; t += buf[j]; }
		}
		Airports[i].num_connected = k;
		if (Airports[i].num_connected > 0)
		{
			if (!(Airports[i].prop_traffic = (float*)calloc(Airports[i].num_connected, sizeof(float)))) ERR_CRITICAL("Unable to allocate airport storage\n");
			if (!(Airports[i].conn_airports = (unsigned short int*) calloc(Airports[i].num_connected, sizeof(unsigned short int)))) ERR_CRITICAL("Unable to allocate airport storage\n");
			for (j = k = 0; j < P.Nairports; j++)
				if (buf[j] > 0)
				{
					Airports[i].conn_airports[k] = j;
					Airports[i].prop_traffic[k] = buf[j] / t;
					k++;
				}
		}
		else
		{
			if (Airports[i].total_traffic > 0)
				fprintf(stderr, "#%i# ", i);
			else
				fprintf(stderr, "%i ", i);
		}
	}
	fclose(dat);
	free(buf);
	fprintf(stderr, "\nAirport data read OK.\n");
	for (i = 0; i < P.Nairports; i++)
	{
		/*		fprintf(stderr,"(%f %i|",Airports[i].total_traffic,Airports[i].num_connected);
		*/		t = 0; k = 0;
		for (j = Airports[i].num_connected - 1; j >= 0; j--)
		{
			if ((Airports[i].prop_traffic[j] > 0) && (Airports[Airports[i].conn_airports[j]].total_traffic == 0))
			{
				t += Airports[i].prop_traffic[j];
				Airports[i].num_connected--;
				if (j < Airports[i].num_connected)
				{
					Airports[i].prop_traffic[j] = Airports[i].prop_traffic[Airports[i].num_connected];
					Airports[i].conn_airports[j] = Airports[i].conn_airports[Airports[i].num_connected];
				}
				Airports[i].prop_traffic[Airports[i].num_connected] = 0;
				Airports[i].conn_airports[Airports[i].num_connected] = 0;
			}
			else if (Airports[i].prop_traffic[j] > 0)
				k = 1;
		}
		/*		fprintf(stderr,"%f %i ",t,k);
		*/		t = 1.0 - t;
		if (k)
		{
			Airports[i].total_traffic *= t;
			t2 = 0;
			for (j = 0; j < Airports[i].num_connected; j++)
			{
				Airports[i].prop_traffic[j] = t2 + Airports[i].prop_traffic[j];
				t2 = Airports[i].prop_traffic[j];
			}
			for (j = 0; j < Airports[i].num_connected; j++)
				Airports[i].prop_traffic[j] /= t2;
			/*			if((Airports[i].num_connected>0)&&(Airports[i].prop_traffic[Airports[i].num_connected-1]!=1))
							fprintf(stderr,"<%f> ",Airports[i].prop_traffic[Airports[i].num_connected-1]);
			*/
		}
		else
		{
			Airports[i].total_traffic = 0; Airports[i].num_connected = 0;
		}
		if (Airports[i].num_connected > 0)
		{
			for (j = k = 0; k < 128; k++)
			{
				t = ((double)k) / 128;
				while (Airports[i].prop_traffic[j] < t) j++;
				Airports[i].Inv_prop_traffic[k] = j;
			}
			Airports[i].Inv_prop_traffic[128] = Airports[i].num_connected - 1;
		}
		/*		fprintf(stderr,"%f) ",Airports[i].total_traffic);
		*/
	}
	fprintf(stderr, "Airport data clipped OK.\n");
	for (i = 0; i < MAX_DIST; i++) AirTravelDist[i] = 0;
	for (i = 0; i < P.Nairports; i++)
		if (Airports[i].total_traffic > 0)
		{
			for (j = 0; j < Airports[i].num_connected; j++)
			{
				k = (int)Airports[i].conn_airports[j];
				traf = floor(sqrt(dist2_raw(Airports[i].loc_x, Airports[i].loc_y, Airports[k].loc_x, Airports[k].loc_y)) / OUTPUT_DIST_SCALE);
				l = (int)traf;
				//fprintf(stderr,"%(%i) ",l);
				if (l < MAX_DIST)
					AirTravelDist[l] += Airports[i].total_traffic * Airports[i].prop_traffic[j];
			}
		}
	sprintf(outname, "%s.airdist.xls", OutFile);
	if (!(dat = fopen(outname, "w"))) ERR_CRITICAL("Unable to open air travel output file\n");
	fprintf(dat, "dist\tfreq\n");
	for (i = 0; i < MAX_DIST; i++)
		fprintf(dat, "%i\t%lg\n", i, AirTravelDist[i]);
	fclose(dat);
}


int GetInputParameter(FILE* dat, FILE* dat2, char* SItemName, char* ItemType, void* ItemPtr, int NumItem, int NumItem2, int Offset)
{
	int FindFlag;

	FindFlag = GetInputParameter2(dat, dat2, SItemName, ItemType, ItemPtr, NumItem, NumItem2, Offset);
	if (!FindFlag)
	{
		fprintf(stderr, "\nUnable to find parameter `%s' in input file. Aborting program...\n", SItemName);
		exit(1);
	}
	return FindFlag;
}

int GetInputParameter2(FILE* dat, FILE* dat2, char* SItemName, char* ItemType, void* ItemPtr, int NumItem, int NumItem2, int Offset)
{
	int FindFlag = 0;

	if (dat2) FindFlag = GetInputParameter3(dat2, SItemName, ItemType, ItemPtr, NumItem, NumItem2, Offset);
	if (!FindFlag)
		FindFlag = GetInputParameter3(dat, SItemName, ItemType, ItemPtr, NumItem, NumItem2, Offset);
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
				else if ((match[0] == '#') && (match[1] == '1') && (match[2] == '2'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP12;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP12;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '1') && (match[2] == '3'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP13;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP13;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '1') && (match[2] == '4'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP14;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP14;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '1') && (match[2] == '5'))
				{
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP15;
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP15;
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

