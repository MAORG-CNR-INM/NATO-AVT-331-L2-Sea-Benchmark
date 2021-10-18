- This is the L2 sea benchmark problem developed by CNR-INM.
- The tool provides the calm-water resistance of the DTMB-5415 (model scale) at Fr=0.28 by
a linear potential flow solver.
- The tool is implemented in fortran and works under linux on windows linux subsystem by GNU or Intel fortran compilers.
- External libraries (lapack and blas for GNU or mkl for Intel) are needed.
- OpenMP implementation is also available but not mandatory.

- To compile the code, 
    -- add the following line at the end of the .bashrc file
        "ulimit -s unlimited"
    -- open the makefile and
        -> for Intel compiler, uncomment the following compiling options
            FC       = ifort
            FCOPT    = -O5
            FCOPTOMP = -qopenmp
            LIB      = -mkl
        -> for GNU compiler, uncomment the following compiling options
            FC       = gfortran   
            FCOPT    = -O3
            FCOPTOMP = -fopenmp
            LIB      = -llapack -lblas
    -- save and close the makefile and finally type "make" in a linux shell

- The potential flow solver can be run at even keel or with 2DoF. For the bechmark pourpouse the simulation can be performed at even and keel and the example file are already set up.
- Up to 7 fidelity levels have been defined.

- Only two text files (SBDF.nml and variables.inp) need to be edited to run the code.
    --> SBDF.nml contains all the namelists
    ----> to select the fidelity levels, the parameters to be edited in the MAIN_PARAMETERS namelist is
        igrid           = 7				! Grid/fidelity level, 1=highest	(int)   
    --> variables.inp contains the design variables
        it is a column text files of 14 lines, each of them represents the design variable value and have to be within -1 and 1

- Run the code in the directory with the input files (see the /example/DTMB-5415 folder) executing the binary file in the /bin folder.
- A CPU000 folder will be created with all the input and output files. The objective function value, along with the costraints can be found in the objective.out file

References are available in the /doc folder
- For the design-space and optimization problem definitions, constraints, and solver refer to
-- Serani, A., Stern, F., Campana, E. F., & Diez, M. (2021). Hull-form stochastic optimization via computational-cost reduction methods. Engineering with Computers, 1-25.
-- Liuzzi, G., Lucidi, S., Rinaldi, F., Pellegrini, R., Serani, A., & Diez, M. (2021). Derivative-Free Line-Search Algorithm for Multi-Fidelity Optimization. In AIAA Scitech 2021 Forum (p. 1237).
-- Serani, A., Diez, M., Wackers, J., Visonneau, M., & Stern, F. (2019). Stochastic shape optimization via design-space augmented dimensionality reduction and rans computations. In AIAA SciTech 2019 Forum (p. 2218).
