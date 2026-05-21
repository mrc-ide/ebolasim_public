# Parameters and parameter files

`ebolasim` uses a large number of input parameters to specify the synthetic population, transmission dynamics, interventions, and simulation parameters. These are then read into and stored in a global struct.

## The parameter struct

The simulation parameters are stored in a `PARAM` struct. This struct is defined in `SpatialSim.h` and is declared as `param P` before `main()`. As a global struct, all functions can access the parameters by referencing `P`, e.g. `P.DoAdunits`.

## Parameter files

To avoid individual parameter files becoming extremely long, we split them up into pre-parameter and parameter files. Pre-parameter files typically contain variables related to the population, disease transmission parameters and simulation parameters (i.e number of runs, size of a timestep, files to output) which don't usually change between scenarios. Parameter files typically contain intervention-specific parameters which do vary between scenarios.

The structure of the pre-parameter and parameter files is:

> \[Parameter name\]

> parameter_value

## Linking the parameter file and struct - `ReadParams()`

The `ReadParams()` function reads in the parameter files and stores the parameters in the struct `P`.

`ReadParams()` itself calls two key functions: `GetInputParameter()`, `GetInputParameter2()`.

`GetInputParameter()` is a wrapper for `GetInputParameter2()`, which in turn calls a third function called `GetInputParameter3()`. `GetInputParameter3()` actually does the string matching and assignment (using a pointer). The key difference between `GetInputParameter()` and `GetInputParameter2()` is that `GetInputParameter()` will throw an error if the parameter name isn't found, whereas `GetInputParameter2()` can be used to set default values if the parameter names isn't found in the parameter file. This can be useful, but can lead to errors and therefore there are some parameters that **must** be specified in the parameter file.

## Adding new parameters

To add a new parameter you need to:

1. Declare it (with type and size - if necessary) in the `PARAM` struct in `SpatialSim.h`.
2. Add it to the `ReadParams()` function making sure to use `GetInputParameter()` or `GetInputParameter2()` as appropriate. If using `GetInputParameter2()`, make sure to add the default value. `ReadParams()` is a long function due to the number of possible parameters, so it is important to think about where to add new parameters, particularly as they may be dependent on other parameters. For example, parameters related to admin units shouldn't be assigned until we have specified that we are subdividing the population into admin units.
3. Add it to the parameter file, being careful to type it exactly as it is in `ReadParams()`. Not that I've ever made that mistake....

### A few other things to think about:

1. If you need to add more command line parameters (`/CLPn:` etc.), this has to be added at the very beginning of `main()` (when all command line parameters are scanned in) and in `GetInputParameter3()` (just copy the template). For two-digit numbers, i.e. CLP1x, the code has to be placed before the corresponding one-digit number, i.e. CLP1; otherwise the `else if` statement will match everything to CLP1.
2. If the parameter needs to be reset between realisations, this will need to be added to `InitModel()`.

