#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <time.h>


// #define COUNTRY_WA
#define DO_OMP_PARALLEL
#define MAX_NUM_THREADS 16
#define CACHE_LINE_SIZE 256
#include <omp.h>

#ifndef UNIX
#define WIN32_BM
#endif
#define STRICT
#ifdef WIN32_BM
#define _WIN32_WINNT 0x0400
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <vfw.h>
#include <gdiplus.h>
#include "avi_utils.h"
#endif
#ifdef IMAGE_MAGICK
#include "Magick++.h"
#endif


#define ERR_CRITICAL(x) {fprintf(stderr,x);exit(1);}
/*
  #define HOST_TREATED(x) (0)
  #define HOST_TO_BE_TREATED(x) (0)
  #define PLACE_TREATED(x,y) (0)
*/
#define FAST_US

#define HOST_TREATED(x) ((Hosts[x].treat_stop_time>ts)&&(Hosts[x].treat_start_time<=ts))
#define HOST_TO_BE_TREATED(x) (Hosts[x].treat_stop_time>ts)
#define PLACE_TREATED(x,y) (Places[x][y].treat_end_time>ts)
#define PLACE_CLOSED(x,y) ((Places[x][y].close_start_time<=ts)&&(Places[x][y].close_end_time>ts))
#define PLACE_TO_BE_CLOSED(x,y) (Places[x][y].close_end_time>ts)
#define HOST_TO_BE_VACCED(x) (Hosts[x].vacc_start_time<USHRT_MAX-1)
#define HOST_VACCED(x) ((Hosts[x].vaccRing==3)?(Hosts[x].vacc_start_time+P.usVaccTimeToEfficacyThirdRing<=ts):(Hosts[x].vacc_start_time+P.usVaccTimeToEfficacy<=ts))
//#define HOST_VACCED(x) (Hosts[x].vacc_start_time+P.usVaccTimeToEfficacy<=ts)
#define HOST_VACCED_SWITCH(x) (Hosts[x].vacc_start_time>=P.usVaccTimeEfficacySwitch)
#define HOST_QUARANTINED(x) ((Hosts[x].quar_comply==1)&&(Hosts[x].quar_start_time+P.usHQuarantineHouseDuration>ts)&&(Hosts[x].quar_start_time<=ts))
#define HOST_TO_BE_QUARANTINED(x) ((Hosts[x].quar_start_time+P.usHQuarantineHouseDuration>ts)&&(Hosts[x].quar_comply<2))
#define HOST_ISOLATED(x) ((Hosts[x].isolation_start_time+P.usCaseIsolationDelay<=ts)&&(Hosts[x].isolation_start_time+P.usCaseIsolationDelay+P.usCaseIsolationDuration>ts))
#define HOST_ABSENT(x) ((Hosts[x].absent_start_time<=ts)&&(Hosts[x].absent_stop_time>ts))

#define PI 3.1415926535
#define EARTHRADIUS 6366707.0
#define BIT15 32768
#define BIT10 1024
#define BIT20 1048576

#define OUTPUT_DIST_SCALE 1000
#define MAX_PLACE_SIZE 20000
#define MAX_NUM_SEED_LOCATIONS 1000

#define MAX_BACKG_COLS 32

#define GEMMA_HELPFUL_LABEL 2
#define DANNY_HELPFUL_LABEL 1
#define PERSON_NOT_IN_PLACE_AT_THAT_TIME -1

#ifdef COUNTRY_UK
#define MAX_DIST 2000
#define NK_HR 400
#define NKR 2000000
#elif defined(COUNTRY_THAILAND)
#define MAX_DIST 2000
#define NK_HR 400
#define NKR 2000000

#else
#define MAX_DIST 1000
#define NK_HR 2500
#define NKR 4000000
#endif

#define CDF_RES 20
#define INFPROF_RES 28
#define RECOVERY_RES 100 //resolution of recovery probability function: ggilani - 23/10/14
/* Max household size of 16 */



#define HH_DISABLED 16384
#define CONTACT_MASK 16383

#define MAX_HOUSEHOLD_SIZE 15
#define MAX_NUM_RESIST_TYPES 3
#define MAX_INTERVENTION_TYPES 1
#define MAX_INTERVENTIONS_PER_ADUNIT 10

#define NUM_PLACE_TYPES 4
#define NUM_PLACE_TYPES_NOAIR 4
#define HOTEL_PLACE_TYPE 5
#define MAX_ADUNITS 500
#define ADUNIT_LOOKUP_SIZE 800000
#define MAX_COUNTRIES 100

#define PLACE_GROUPS_PER_INT 2
#define NUM_AGE_GROUPS 17
#define AGE_GROUP_WIDTH 5
#define DAYS_PER_YEAR 364
#define INFECT_TYPE_MASK 16
#define MAX_GEN_REC 20
#define MAX_SEC_REC 500
#define INF_QUEUE_SCALE 5
#define MAX_TRAVEL_TIME 14
#define MAX_TRAVEL_RATIO 50
#define MAX_INFECTIOUS_STEPS 255
#define MAX_DUR_IMPORT_PROFILE 10245
#define MAX_CHANGE_POINTS 100
#define MAX_RING_SIZE 1000000 //added this for vaccination code
#define NUM_VACCDOSE_GROUPS 11
#define VACCDIST_WIDTH 1000
#define POPDIST_WIDTH 1
#define VACCDOSE_WIDTH 40
#define NUM_VACCDIST_GROUPS 11
#define NUM_POP_GROUPS 11
#define VACCDOSECELL_WIDTH 400
#define NUM_VACCDOSECELL_GROUPS 11
#define VACCDOSE_GROUP(x) (((x/VACCDOSE_WIDTH>=NUM_VACCDOSE_GROUPS))?(NUM_VACCDOSE_GROUPS-1):(x/VACCDOSE_WIDTH))
#define VACCDOSECELL_GROUP(x) (((x/VACCDOSECELL_WIDTH>=NUM_VACCDOSECELL_GROUPS))?(NUM_VACCDOSECELL_GROUPS-1):(x/VACCDOSECELL_WIDTH))
#define VACCDIST_GROUP(x) (((x/VACCDIST_WIDTH>=NUM_VACCDIST_GROUPS))?(NUM_VACCDIST_GROUPS-1):(x/VACCDIST_WIDTH))
#define POPDIST_GROUP(x) (((x/POPDIST_WIDTH>=NUM_POP_GROUPS))?(NUM_POP_GROUPS-1):(x/POPDIST_WIDTH))


#define MAX_AIRPORTS 5000
#define NNA 10
#define MAX_DIST_AIRPORT_TO_HOTEL 200000.0
#define MIN_HOTELS_PER_AIRPORT 20
#define HOTELS_PER_1000PASSENGER 10
#define PLACE_CLOSE_ROUND_HOUSEHOLD 1

//added for radiation model
#define RAD_INVCDF 100
#define RAD_TRIPS_PP 5
#define RAD_TEMP_MEM_SCALE 1

//added for fixing seeds for janetta - 08/03/2023
#define MAX_FIXED_SEEDS 1000

/*
  #define NO_TREAT_PROPH_CASES
*/

#define HOST_AGE_YEAR(x) (Hosts[x].age)
#define HOST_AGE_GROUP(x) (Hosts[x].age/AGE_GROUP_WIDTH)

#define MAXINTFILE 10
