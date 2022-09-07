!>\file mainpage.f90
!! This file is just to create the documentation page. No code here.
!
!>
!! @mainpage Using DSEM
!! 
!! @section Sections Sections
!! -# @ref Structure
!! -# @ref Case
!! -# @ref Running
!! 		-# @ref Serial
!! 		-# @ref Parallel
!! -# @ref Samples
!! 		-# @ref RCEV
!!		-# @ref reacRC
!!		-# @ref Sod
!! -# @ref Visualization
!!
!! @section Structure The Code Structure
!! The Discontinuous Spectral Element Method (DSEM) code is a high-order flow solver that solves the three-dimensional compressible Navier-Stokes
!! equations and is capable of simulating turbulent, supersonic, and reactive flows. DSEM uses unstructured grids and hence, handles complex 
!! geometries. The code is written in Fortran 90.
!!
!! The code consists of four folders:
!! -# <b>src:</b> contains all source code files.
!! -# <b>runFolder:</b> contains files with adjustable parameters to run the code.
!! -# <b>caseFiles:</b> contains files related to each specific case.
!! -# <b>meshing:</b> contains DSEM mesh generator.
!!
!! @section Case Setting Up a Case
!! -# <b>Make a grid.</b> DSEM uses three-dimensional hexahedral elements. You can directly make the grid using the DSEM mesh generator, or generate a grid with any other tool and convert it to the DSEM for- mat. The mesh file should be placed in the runFolder. Boundary conditions and the polynomial order of the elements are also defined in the mesh file. The number indicating each type of boundary condition is defined in the“runFolder/backstep surfaces” file. There are six types of boundary conditions:
!!      -# <b>inflow:</b> Inflow
!!      -# <b>outflow:</b> Dirichlet outflow boundary condition
!!      -# <b>walladiab:</b> Adiabatic wall
!!      -# <b>wallisoth:</b> Iso-thermall wall (set temperature in the bump file)
!!      -# <b>symmetry:</b> Symmetry plane
!!      -# <b>transmiss:</b> Supersonic boundary condition
!! Setting the periodic boundary condition is explained next.
!! -# <b>Set periodic boundary conditions.</b> If you have any periodic bound- ary condition, you need to set it in the problem file (caseFiles/Back- stepProblem.f90). You can activate the periodic boundary condition in each direction by setting the corresponding variable to 1.
!!      -# z direction: 
!! @code zperiodic = 1 @endcode
!!      -# y direction: 
!! @code yperiodic = 0 @endcode
!!      -# x direction: 
!! @code chicka = 0 @endcode
!! -# <b>Add source terms.</b> If you have any source term, you can add them in the subroutine compute source in the problem file. The variable source is a 5-dimension array in which you can put the source terms for the 5 Navier-Stokes equations. You also need to activate the source terms in the bump file (explained later).
!! -# <b>Set the size of the problem.</b> The size file (caseFiles/size.f90) needs to be updated by the size of the problem. The variable ngp is the number of elements, and nx, ny, nz, and nmax are the polynomial orders in each direction and the maximum of the order in all three directions.
!! -# <b>Initialize the solution.</b> The initial condition of the flow is set in the subroutine exactsol in the problem file. Variables rho, rhou, rhov, rhow, rhoe are ρ, ρu, ρv, ρw, ρe, respectively.
!!
!! @section Running Running the Code
!! DSEM may be either run in serial or parallel. However, it is always necessary to run the code in serial to generate a restart file and <a href="http://glaros.dtc.umn.edu/gkhome/metis/metis/overview">METIS</a> file. The <a href="http://glaros.dtc.umn.edu/gkhome/metis/metis/overview">METIS</a> file is used for partitioning the grid across many cores on a supercomputer. 
!! 
!! @subsection Serial Running the Code in Serial
!! In order to run the code, first, you need to set the running parameters in the bump file (runFolder/bump.in). Note that you can comment out a line by adding “!” at the beginning of the line.
!! -# <b>Name your case and set your mesh file:</b>
!! @code
!! run file name = shear
!! mesh file = Backstep.mesh
!! @endcode
!! \n
!! -# <b>Set the final time and/or the maximum number of time steps.</b>
!! @code
!! final time = 100.d0
!! maximum number of steps = 1000
!! @endcode
!! \n
!! -# <b>Set the CFL number.</b>
!! @code
!! use a cfl = 0.9d0
!! @endcode
!! \n
!! -# <b>Specify if you want a plot file at the end of the simulation or movie plot files during the simulation by including the corresponding line.</b>
!! @code
!! write endplot file
!! write movie file
!! @endcode
!! \n 
!! -# <b> Set the frequency of the movie files (in terms of time steps), if you have chosen to generate movie files.</b>
!! @code output frequency = 200 @endcode
!! \n
!! -# <b>Choose constant time step size (time accurate computation) or vari- able time step size (steady-state computation) by including the cor- responding line.</b>
!! @code
!! steady-state computation 
!! time accurate computation
!! @endcode
!! \n
!! -# <b>Activate the source term (if you have any) by replacing “none” by any other word.</b>
!! @code
!! source terms = none
!! source terms = time_dependant
!! source terms = fmdf
!! @endcode
!! \n
!! -# <b>Set the Prandtl number, Reynolds number, wall temperature (if using isothermal wall), and outlet pressure (if using transmissive boundary condition).</b>
!! @code
!! prandtl number = 0.72d0
!! reynolds number = 10000.d0 
!! wall temperature = 0.25d0 
!! outlet pressure = 17.568d0
!! @endcode
!! \n
!! -# <b>After you set the running parameters in the bump file, you need to compile the code.</b> 
!!      - You can compile the code by executing the following command in the run folder. Note that the address of the case folder is specified in the file “compileIt” and the compiler is defined in the Makefile in the case folder.
!! @code{.sh} $ ./compileIt @endcode
!!      - After compiling the code, an executable called “3DNS” should be generated. You may now run the code by executing this file.
!! @code{.sh} $ ./3DNS @endcode
!!      - All the plot files and the restart file will be generated in the runFolder.
!!
!! @subsection Parallel Running the Code in Parallel
!! In order to run the code in parallel, the code needs to read the initial condition from a restart file. You can make a restart file by running the code in serial for 0 time step.
!!
!! Once you have a restart file to read from, running the code in parallel is similar to running the code in serial, with a few differences:
!! 
!! -# Instead of naming the case and specifying the mesh file in the bump file, you specify the restart file, from which you read the initial condition:
!! @code
!! run file name = shear.rst 
!! mesh file = shear.rst
!! @endcode
!! \n 
!! -# You will need a METIS file. The METIS file determines which elements each core is responsible for. It is a list of numbers, with each number in a separate line. The number in the <em>n</em>th line is the label of the core which is responsible for the element <em>n</em>.
!!      - METIS is a third party software library available here: <a href="http://glaros.dtc.umn.edu/gkhome/metis/metis/overview">METIS - Serial Graph Partitioning and Fill-reducing Matrix Ordering</a>
!!      - You can generate the metis file using “gpmetis”. After you run the code in serial for 0 time step, the code generates a connectivity list. The file is named “yourcasename.graph”. For instance, you can generate a metis file for 400 cores by executing the following command:
!! @code{.sh} $ gpmetis yourcasename.graph 400 @endcode
!!      - You can run the code in parallel using the following command:
!! @code{.sh} $ mpirun -n 400 3DNS @endcode
!!
!! @section Samples Running Sample Cases
!! The following section will outline how to run the sample cases provided in the run folder and case folders for the code.
!! 
!! @subsection RCEV Supersonic Ramp Cavity
!! The supersonic ramp cavity case requires the user to first create a restart file. Once the restart file is created, it is then run with high turbulence viscosity (\f$ C_s=0.45 \f$) for about 10 flow-through times to establish a fully developed flow. Then, it can be run with lower turbulent viscosity (\f$ C_s=0.15 \f$).
!!
!! To do this, we follow these steps:
!! -# Specify the mesh file and case name in the bump.in file and comment out the lines for reading from restart file.
!! @code
!! ! ----------0 timesteps------------
!! run file name = shear
!! mesh file     = 78165_topwall.mesh
!! 
!! ! --------cluster(rst file)--------
!! ! run file name = shear.rst
!! ! mesh file     = shear.rst
!! @endcode
!! -# We then setup the case to run for 'zero timesteps'.
!! @code
!! maximum number of steps = 0
!! write endplot file
!! save restart file at end
!! @endcode
!! -# We then compile and run the code for zero timesteps by executing the 3DNS binary.
!! @code{.sh}
!! $ ./compileIt
!! $ ./3DNS
!! @endcode
!! -# Now the METIS and restart files are generated. The code may be run in parallel. We will run the code for 10 flow-through times with \f$ C_s=0.45 \f$.  Make the following modifications to the bump.in file.
!! @code
!! ! ----------0 timesteps------------
!! ! run file name = shear
!! ! mesh file     = 78165_topwall.mesh
!! 
!! ! --------cluster(rst file)--------
!! run file name = shear.rst
!! mesh file     = shear.rst
!! @endcode
!! @code
!! maximum number of steps = 22640
!! output frequency = 4528
!! write movie file
!! @endcode
!! @code
!! smagorinsk      ! use smaogorinsly model
!! rhoSens         ! apply rho sensor
!! shockSens       ! apply shock sensor
!! Cs = 0.45d0     ! smaogorinsky coefficient
!! numax = 60.d0   ! maximum turbulent [kinematic] viscosity
!! @endcode
!! -# We then generate the METIS file for 96 cores, if it does not already exist in the run folder.
!! @code{.sh} $ gpmetis shear.graph 96 @endcode
!! -# We then execute the code to run in parallel on 96 cores using MPIRUN.
!! @code{.sh} $ mpirun -n 96 3DNS @endcode
!! -# After the code has terminated with no faults, we then return to the bump.in file to reduce the Smagorinsky coefficient and update the restart file's name from which the code starts the simulation.
!! @code
!! Cs = 0.15d0     ! smaogorinsky coefficient
!! @endcode
!! @code
!! run file name = shear+.rst
!! mesh file     = shear+.rst
!! @endcode
!! -# We also need to rename the METIS file so that the name is compatible with the restart file's name. Then, we can run the code for another 10 flow-through times.
!! @code{.sh}
!! $ mv shear.graph.part.96 shear+.graph.part.96
!! $ mpirun -n 96 3DNS
!! @endcode
!! -# Once the code terminates you will have the resulting plot file of the simulation.
!! @image html "../../runFolder/2.Bfs/SampleResults/u.jpg" "Ramp Cavity U-Velocity" width=100%
!! @image html "../../runFolder/2.Bfs/SampleResults/Nu_t.jpg" "Ramp Cavity Turbulent Viscosity" width=100%
!! @image html "../../runFolder/2.Bfs/SampleResults/EV.jpg" "Ramp Cavity Entropy Viscosity" width=100%
!! 
!! @subsection reacRC Reacting Subsonic Ramp Cavity
!! The reacting subsonic ramp cavity case requires the user to first create a restart file. Once the restart file is created, it is then run for cold flow to establish the flow field prior to fuel injection. Once the flow field is established, we initialize the Monte Carlo particles with initial conditions interpolated from the Eulerian flow field. It is then that fuel injection begins and we run the reacting case with the Filtered Mass Density Function Method. 
!!
!! To do this, we follow these steps:
!! -# Specify the mesh file, case name, and particle file in the bump.in file.
!! @code
!! run file name = rcInjector
!! mesh file     = backstep.mesh
!! geometry file = backstep_surfaces
!! droplet file  = shear.part
!! @endcode
!! -# We then setup the case to run for 'zero timesteps' and ensure that all FMDF related options are turned off.
!! @code
!! maximum number of steps = 0
!! @endcode
!! -# We then run the code for zero timesteps by executing the 3DNS binary.
!! @code{.sh}
!! $ ./3DNS
!! @endcode
!! -# Now the METIS and restart files are generated. The code may be run in parallel. We will run the cold flow in parallel for a duration of 200 time units.  Make the following modifications to the bump.in file.
!! @code
!! run file name = rcInjector.rst
!! mesh file     = rcInjector.rst
!! final time    = 200.0d0
!! maximum number of steps = 10000000
!! output frequency = 100
!! read metis file
!! smagorinsk
!! @endcode
!! -# We then execute the code to run in parallel on 8 cores using MPIRUN
!! @code{.sh} $ mpirun -n 8 3DNS @endcode
!! -# After the code has terminated with no faults, we then return to the bump.in file to turn on FMDF and run the code from the new restart file.
!! @code
!! run file name = rcInjector+.rst
!! mesh file     = rcInjector+.rst
!! source terms = fmdf
!! fmdf
!! restfldr
!! drops
!! @endcode
!! -# Once the code terminates you will have the resulting plot file of the simulation.
!! @image html "../../runFolder/fmdfTests/rc4812ele1injector/SampleResults/mach.png" "Ramp Cavity Mach Number" width=100%
!! @image html "../../runFolder/fmdfTests/rc4812ele1injector/SampleResults/temp.png" "Ramp Cavity Temperature" width=100%
!! @image html "../../runFolder/fmdfTests/rc4812ele1injector/SampleResults/Fuel.png" "Ramp Cavity Fuel" width=100%
!! @image html "../../runFolder/fmdfTests/rc4812ele1injector/SampleResults/reaction.png" "Ramp Cavity Reaction" width=100%
!! 
!! @subsection Sod Reacting Sod Problem
!! The reacting sod case requires the user to first create a restart file. Once the restart file is created, we initialize the Monte Carlo particles with initial conditions interpolated from the Eulerian flow field, there is no need to run the cold flow before starting reaction like the previous case.
!!
!! To do this, we follow these steps:
!! -# Specify the mesh file, case name, and particle file in the bump.in file.
!! @code
!! run file name = sod
!! mesh file     = sod200.mesh
!! geometry file = backstep_surfaces
!! droplet file  = shearnew.part
!! @endcode
!! -# We then setup the case to run for 'zero timesteps' and ensure that all FMDF related options are turned off.
!! @code
!! maximum number of steps = 0
!! @endcode
!! -# We then run the code for zero timesteps by executing the 3DNS binary.
!! @code{.sh}
!! $ ./3DNS
!! @endcode
!! -# Now the METIS and restart files are generated. It is not recommended to run this case in parallel, as it is a computationally quick case, and parallel communication would outweigh computational time.  Make the following modifications to the bump.in file. After the code has terminated with no faults, we then return to the bump.in file to turn on FMDF and run the code from the new restart file.
!! @code
!! run file name = sod.rst
!! mesh file     = sod.rst
!! final time    = 0.20d0
!! maximum number of steps = 10000000
!! endplot
!! source terms = fmdf
!! shock
!! fmdf
!! restfldr
!! drops
!! @endcode
!! -# Once the code terminates you will have the resulting plot file of the simulation.
!! @image html "../../runFolder/fmdfTests/Sod/SampleResults/Non-reacting/Comparison.png" "Non-reacting Comparison" width=100%
!! @image html "../../runFolder/fmdfTests/Sod/SampleResults/Reacting/Comparison.png" "Reacting Comparison" width=100%
!!
!! @section Visualization Visualization
!! To visualize the results you may use any of the following programs:
!! -# Free Software
!! 		-# <a href="https://wci.llnl.gov/simulation/computer-codes/visit/">VisIt Visualization Software</a>
!! 		-# <a href="http://www.paraview.org/">ParaView</a>
!! -# Licensed Software
!!		-# <a href="http://www.tecplot.com/products/tecplot-360/">TecPlot 360</a>
!!
!! We strongly recommend using VisIt as it has been tested extensively with DSEM, however these and any other visualization packages that can read TecPlot ASCII and binary PLT formats should work.