# NorCPM1  


## About

Content: CMIP6 DCPP version of the Norwegian Climate Prodiction Model 

Version: 1.0.0 

Updated: 2021.03.08 (Ingo.Bethke@uib.no) 

## Installation and porting 

### System requirements 

* 350 GB of free space for input data 
* linux machine with minimum 32 cores
* compiler suite (preferably intel)
* mpi and netcdf libraries
* netcdf climate operators (NCO)

### Installing initial conditions and other input data 

External forcing and other input data required to rerun NorCPM1 CMIP6 DCPP 
simulations are obtained from https://doi.org/10.11582/2021.00013

Copy and untar files to a location of your choice accessible to compute nodes 
on your HPC facility. This will create a top-level folder inputdata and place 
files in respective subfolders.

### Porting of assimilation code

The assimilation code is located in NORCPMROOT/assim/norcpm1, where NORCPMROOT 
is NorCPM's installation directory.

Files that need to be modified are:
* assim/norcpm1/shared/make.inc.betzy
* assim/norcpm1/micom_init/Makefile.config.betzy
* setup/norcpm1/settings/setmach.sh (see Section "Porting of setup scripts")

To test compilation of i1 and i2 assimilation code, change directory to 
setup/norcpm1 and then run 

    scripts/build_assimcode.sh settings/norcpm1_assim-i1.sh 
    scripts/build_assimcode.sh settings/norcpm1_assim-i2.sh 

### Porting of model code

NorCPM's Earth system component needs to be ported in 
NORCPMROOT/model/norcpm1/scripts/ccsm_utils/Machines

#### Machine name and characteristics - config_machines.xml

Define the characteristics of your HPC system in config_machines.xml.

Edit/duplication the section for the machine betzy:

         <machine MACH="betzy"
         DESC="NTNU, Trondheim, 128 pes/node, batch system is SLURM"
         EXEROOT="/cluster/work/users/$CCSMUSER/noresm/$CASE"
         OBJROOT="$EXEROOT"
         INCROOT="$EXEROOT/lib/include" 
         LOGDIR="/cluster/work/users/$CCSMUSER/archive/logs"
         DIN_LOC_ROOT_CSMDATA="/cluster/projects/nn9039k/inputdata/subset/norcpm1-cmip6"
         DIN_LOC_ROOT_CLMQIAN="UNSET"
         DOUT_S_ROOT="/cluster/work/users/$CCSMUSER/archive/$CASE"
         DOUT_L_HTAR="FALSE"
         DOUT_L_MSROOT="UNSET"
         DOUT_L_MSHOST="UNSET"
         CCSM_BASELINE="UNSET"
         CCSM_CPRNC="UNSET"
         OS="Linux"
         BATCHQUERY="UNSET"
         BATCHSUBMIT="sbatch" 
         GMAKE_J="8" 
         MAX_TASKS_PER_NODE="128"
         MPISERIAL_SUPPORT="FALSE" />


Set MACH to an acronym/name tag of your choice of your HPC system. Update 
DESC correspondingly.

Set EXEROOT to the root location where the model should create the build and 
run-directories. The model will replace $CCSMUSER with your unix user and $CASE 
with the specific case name (=simulation/experiment name).

Set DIN_LOC_ROOT_CSMDATA to the location where you installed the forcing and 
boundary condition data (this path should end with “/inputdata”).

Set DOUT_S_ROOT to the root location for the short-term archiving (normally on 
the work disk area of the HPC system).

Optionally, set DOUT_L_MSROOT to the root location for the long-term archiving 
(normally a location that is not subject to automatic deletion). 
Leave unchanged if you don't plan to use automatic long-term archiving (e.g., 
if you want to move the data manually).

Optionally, set DOUT_L_MSHOST to the name of the (remote)-server for long-term 
archiving. 
Leave unchanged if you don't plan to use automatic long-term archiving.

Set BATCHSUBMIT to the submit command on your HPC system.

Set GMAKE_J to the number of make instances run in parallel when building the 
system. Set to 1 if licence or memory issues occur.

Set MAX_TASKS_PER_NODE to the maximum number of MPI tasks that can run on a 
single node on your HPC system (usually the same as the number of cores per 
node).

#### CPU configurations in config_pes.xml

Define NorCPM1 CPU-configurations for your HPC system in config_pes.xml.

We recommend to start with replacing all instances of betzy with your choice 
for MACH (i.e. the name tag of your machine). If you reuse the label betzy for
your machine then likely no changes are needed to the cpu configuration file. 

#### Building and runtime options in env_machopts.$MACH

Copy env_machopts.betzy to env_machopts.$MACH, where $MACH should be replaced 
with the name of your machine.

Edit settings as necessary. The settings should make compilers, libraries, 
queuing system commands ect. available during building and model execution. 
You will likely have to remove/replace all module specifications, while the 
runtime environment settings should work for most systems.

#### Compiler and linker options 

Copy Macros.betzy to Macros.$MACH if your compiler is intel. If your compiler 
is pgi, then instead use Macros.hexagon as template.

Following lines need to be customised:

    FC            := FORTRAN COMPILER COMMAND
    CC            := C COMPILER COMMAND
    NETCDF_PATH   := ROOT PATH FOR NETCDF LIBRARY 
    MPI_ROOT      := ROOT PATH FOR MPI LIBRARY 
    INC_NETCDF    := $(NETCDF_PATH)/include
    LIB_NETCDF    := $(NETCDF_PATH)/lib
    MOD_NETCDF    := $(NETCDF_PATH)/include
    INC_MPI       := $(MPI_ROOT)/include
    LIB_MPI       := $(MPI_ROOT)/lib

### Porting of setup scripts

In setup/norcpm1/settings/setmach.sh, following should be specified
* MIN_NODES: minimum allowed number of nodes per job 
* TASKS_PER_NODE: number of tasks per node
* WORK: location user-work directory  
* INPUTDATA: location of the installed model input data  
* modules and environmental variables for MPI, NetCDF, NCO, etc
* stack size limit (e.g. "ulimit -s unlimited")  

In case your machine does not use the SLURM job environment then also 
setup/norcpm1/scripts/submit_experiment.sh need to be modified. 

## Setting up and running experiments 

Scripts and setting files for creating and running NorCPM1 experiments are 
available in NORCPMROOT/setup/norcpm1. It is recommended (though not required) 
to execute the scripts from NORCPMROOT/setup/norcpm1. 

### Setting up an experiment

To set up a new experiment, run 

    ./setup_experiment.sh <path to experiment settings file> [VAR1=value1 VAR2=value2 ...]

where `<path to experiment settings file>` points to one of the experimental 
settings files in the settings sub-directory. 

For example 

    ./setup_experiment.sh settings/norcpm1_historical.sh ENSSIZE=15
    
will set up a no-assimilation historical experiment with 15 simulation members.

Optional VAR=value arguments can be specified to override the defaults from the 
settings file. The settings can also be modified directly but this would 
complicate git-updates to future versions of NorCPM.

By default, the experiment name is derived from the name of the settings file 
but can be set explitely changed by specifying `CASE_PREFIX=<experiment name>`. 

The setup script will prepare individual configuration directories for the 
simulation members in `NORCPMROOT/cases/norcpm1/$CASE_PREFIX` and individual 
run directories in `WORK/noresm/$CASE_PREFIX`.  

### Running the experiment 



### Miscellaneous 

#### Customizing simulation output etc 

-------------------------------------------------------------------------------

The directory of the first simulation member serves as a template for all other
simulation members.

If you want to apply any changes of the code or diagnostic output configuration 
to simulation members then do the following steps. 

Prepare the experiment as usual, for example with  

    ./setup_experiment.sh settings/norcpm1_historical.sh

Change to the configuration directory of the first simulation member in 
`NORCPMROOT/cases/norcpm1/$CASE_PREFIX`.

Make modifications to the code (by placing alternative code in SourceMods sub-
directory and then rebuild) and/or modifications to the output configuration 
(edit bld.nml.csh-files in Buildconf directory).

Rerun the setup_experiment.sh script but this time with argument SKIP_CASE1=1.
For example, 

    ./setup_experiment.sh settings/norcpm1_historical.sh SKIP_CASE1=1 

The SKIP_CASE1 argument will force the script skip the configuration of the 
first simulation members and then clone the settings of the existing first 
simulation members into the other simulation members. 

#### Reference for settings variables  

## Resources

