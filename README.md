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
> latitude, longitude, population, country, admin code
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

### 