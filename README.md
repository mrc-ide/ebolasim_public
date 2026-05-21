# ebolasim_public

ebolasim is Imperial College's individual-based stochastic simulation of Ebola Virus Disease (EVD) transmission.

It has been used to simulate intervention scenarios during the West Africa Ebola outbreak, the 2018-2020 North Kivu outbreak, and more recently it has been used to inform the size of the global Ervebo vaccine stockpile size.

It includes a range of interventions, including contact tracing, hospitals and Ebola treatment centres, safe and dignified burials, and vaccination.

## Data inputs and outputs

### Input files

In the most minimal form, the simulation needs 7 input parameters:
- a parameter file (identified by the prefix 'P/:')
- an output path and file name (identified by 'O/:')
- a density file, providing information of the population of each cell in the simulation (identified by 'D/:')
- four seeds, which are used to seed the random number generator in the set-up and running phases

If fewer than 7 inputs are provided, the model will show an error message indicating the required syntax. Note that the four seeds should always be the final inputs

Other optional input files and parameters can also be provided and specified before the seeds. These are:
- Pre-parameter file ('/PP:). This can be useful to separate out variables that do not change between scenarios(e.g. population and simulation parameters) from intervention variables that do change, and can prevent individual parameter files from becoming really long.
- Network files ('/S:' to save the network file when originally generated and '/L:' to load the previously generated network file in subsequent simulations). These store networks formed by households, schools, workplaces, hospitals, or other networks. In ebolasim, we model extended family networks (i.e. households that are linked together) and hospitals. 
- Air travel files ('/A:'). These can be used to specify air travel networks within or between countries but we have not used then in ebolasim to date.
- School files ('/s:' - note the small s :-)). If exact data on the location and capacity of schools is available, this can be used to generate school networks instead of the heuristic algorithm that is used by default. As for the air travel files, we have not used these in ebolasim to date.
- R0 scaling ('/R:'). This is used to scale transmission up or down relative to the baseline value specified in the parameter file.
- Command line initialisation ('/CLP#:' where # ranges from 1 to 14). This allows us to vary parameters from the command line without having to generate a huge number of parameter files.

#### Density files

The first time the simulation is run, a population density text file is used as the input. Each line of the text file contains the following information for each cell:

`1atitude, longitude, population, country, admin code`

The name of the .txt file is put after the '/D:' and in addition, a '/M:' flag can be specified to a filename to save the population density information in a .bin binary file. Once the binary file has been generated, this can be used as the input as it loads faster.

#### Network files

The simulation code generates a network if households and/or places (i.e. extended family networks and/or hospital) are used in the model. This can be time intensive, so it is best to save the network to a .bin file using the '/S:' flag the first time this is done. On subsequent model run, the saved network file can be loaded using the '/L:' flag.

### Output files

The standard output file for each realisation of a simulation is aggregated over all spatial units and provides counts of the S, L, I, R, D compartments over time as well as the incidence of other variables of interest, such as infections, cases, detected cases, hospitalistions etc.

Additional files can also be output, including:
- summary files which average the metrics over all non-extinct and extinct realisations separately
- outputs disaggregated by administrative units
- outputs disaggregated by age group
- outputs disaggregated by healthcare worker (HCW) status
- outputs disaggregated by country
- infection trees
as well as a few more (infection types, interventions etc.) that can be turned on or off in the pre-parameter file but predate ebolasim and have not been used with or updated for ebolasim. 

## What the model does

This is a brief overview of the key functions called in `main()`. I will try and add further details in separate doc pages soon - both on function and on editing.

### Read in parameter files

After scanning in the command line parameters and setting up the threads (for parallelisation), the next function that is called is `ReadParams`. The purpose of this function is fairly self-explanatory as it reads in all of the parameters specified in both files. If one of the parameter values listed in the parameter files is `#n`, where n is a number between 1 and 14, the code will assign the corresponding `CLP#n` that was specified in the command line.

Within the `ReadParams` file, the code looks for names of parameters in the parameter and pre-parameter files, and if found, assigns it to the relevant variable in the global parameter struct.

**More details here soon!**

### Set up the model

`SetupModel` performs several key functions that only need to be performed once:

1. Read in the density file, determine the bounding box for the simulation space and the number of cells
2. Set up the population
3. Allocate memory for outputs
4. Initialise the spatial kernel (`Init
5. Assign the population to places (or load in a network file)
6. Stratify the population within places
7. Allocate key workers
8. Calibrate the infectiousness in the system and values of beta (i.e. infectiousness) in households, places and community transmission

#### Set up the population

#### Assign the population to places

#### Calibrate the infectiousness of the system

### Within a loop over the number of realisations specified in the parameter file...

### Initialise the model

`InitModel` (re)initialises parameters that change over the course of a single realisation and need to reset before the next realisation begins. This includes bookkeeping variables for counting cases, deaths, vaccine doses etc., and intervention capacities, such as ETU beds available per admin unit and number of cases that can be contact traced per day. It also resets variables associated with individuals (time of infection, outcome, time of outcome etc), households, etc.

**More details here soon!**

### Run the model

`RunModel` runs the main simulation loop. The model records samples at the frequency specified in the parameter file, but we pretty much always have more than one timestep per sample.

For each sample, the code:

1.	Records a sample
2. 	Updates some intervention parameters
3. 	For every timestep within a sample:
	-	Checks for new importation/infection seeding events
	-	Calls `InfectionSweep` (**super important function for transmission events**)
	- 	Calls `IncubRecoverySweep` (transfers individuals from infected to infectious, and from infectious to outcome)
	- 	Depending on interventions used: 
		-	Calls `HospitalSweep` (to admit/discharge individuals from hosptials and ETUs, if available)
		-	Calls `ContactTracingSweep` (to add/remove individuals from contact tracing lists)
		-	Calls `VaccSweep` (to vaccinated individuals in the queue if there are enough doses available)
	-	Calls `TreatSweep` - note that this was originally more of a legacy function in ebolasim as when we first started developing it. At the time, cell-based interventions wer not particularly relevant and it made more sense to link contact tracing and ring vaccination to individuals rather than geographical cells. This has slightly changed over time as we now also include geographically targeted vaccination, and could be updated in future.
	- Calls `TravelReturnSweep` - this (currently) isn't used in ebolasim
4. Calls  `RecordInfTypes` - this can be used to track where infections were made, but needs more documentation

**More details here soon!**

### Save model outputs

`SaveResults` writes all specified results to file, such as the aggregated case numbers, spatially-aggregated case numbers, age-aggregated case numbers.

After all realisations have run, `SaveSummaryResults` calculates the mean and variance of output parameters for extinct and non-extinct realisations, and saves mean results to file.

**More details here soon!**