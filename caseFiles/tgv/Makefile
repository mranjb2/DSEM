ARCH    = $(shell uname)
PROG    = 3DNS

HDRDIR  = ./

FC      = mpif90
LD      = mpif90
src     = ../../src

FFLAGS	= -O3 -m64 -mcmodel=medium
LDFLAGS	= -O3 -m64 -mcmodel=medium

FFTW    = /home/mranjb2/fftw/include/
FFTWL   = -lfftw3

.SUFFIXES: .f90 .f95

.f90.o:
	$(FC) -I $(FFTW) $(FFLAGS) -c $< $(FFTWL)

.f95.o:
	$(FC) $(FFLAGS) -c $<

OBJS    = size.o\
     $(src)/paraFMDF/initializationsc2.o\
     $(src)/paraFMDF/NS_data.o\
     $(src)/paraFMDF/particleRoutines/randomGenerationLibrary.o\
     $(src)/paraFMDF/particleRoutines/random.o\
     BackstepProblem-stoch.o\
     $(src)/paraFMDF/input_lib.o\
     $(src)/Parabolic3D/NS/NS/PhysLib.o\
     $(src)/Parabolic3D/NS/NS/stac3m/dtcalc.o\
     $(src)/turbulence/staticSmago/staticSmag.o\
     $(src)/paraFMDF/Plotting/WrtPlotFiles_scalar_turb_3d.o\
     $(src)/paraFMDF/Plotting/scatterPlot.o\
     $(src)/3DDriver_parfmdf.o\
     $(src)/paraFMDF/fmdfDriver2.o\
     $(src)/Grid_Routines/stac3m_geomLib.o\
     $(src)/Grid_Routines/refine_geometry.o\
     $(src)/Grid_Routines/element_library.o\
     $(src)/Grid_Routines/face_operations.o\
     $(src)/Grid_Routines/GeomLib3D.o\
     $(src)/Grid_Routines/read_mesh.o\
     $(src)/Grid_Routines/set_domain_geometry.o\
     $(src)/Grid_Routines/set_mortars.o\
     $(src)/Parabolic3D/stac3m/restart.o\
     $(src)/Parabolic3D/stac3m/fluxesH_newa.o\
     $(src)/Parabolic3D/stac3m/Vis_Fluxes_new.o\
     $(src)/Parabolic3D/stac3m/SolnDerivsc.o\
     $(src)/Parabolic3D/stac3m/send_q.o\
     $(src)/Parabolic3D/mortar_operations.o\
     $(src)/Parabolic3D/stac3m/send_flux.o\
     $(src)/Parabolic3D/NS/BCs/stac3m/wallbc_inv.o\
     $(src)/Parabolic3D/NS/BCs/stac3m/wallbc_isoth2.o\
     $(src)/Parabolic3D/NS/BCs/stac3m/wall_slide.o\
     $(src)/Parabolic3D/NS/BCs/stac3m/inflow_specAll.o\
     $(src)/Parabolic3D/NS/BCs/stac3m/outflow_specAll.o\
     $(src)/Parabolic3D/NS/BCs/stac3m/periodm_specAll.o\
     $(src)/Parabolic3D/NS/BCs/stac3m/periods_specAll.o\
     $(src)/Parabolic3D/NS/BCs/stac3m/ViscousBC_slide.o\
     $(src)/Parabolic3D/NS/BCs/stac3m/transmissive.o\
     $(src)/Independentroutines/baselib.o\
     $(src)/Independentroutines/io_lib.o\
     $(src)/Independentroutines/legendre_routines.o\
     $(src)/Independentroutines/Fft.o\
     $(src)/stac3m/updatea.o\
     $(src)/Spectral/mortarLib.o\
     $(src)/Spectral/matmula.o\
     $(src)/Spectral/stac3m/FDerivLibb.o\
     $(src)/Spectral/stac3m/FInterpLib.o\
     $(src)/Spectral/stac3m/FInterpLiba.o\
     $(src)/Spectral/stac3m/matrices.o\
     $(src)/Statistics/dongru/all_statsc.o\
     $(src)/Statistics/dongru/pre_stats_backstep_Jovic.o\
     $(src)/paraFMDF/particleRoutines/WrtPartFile-wVis.o\
     $(src)/Particles/part_bc.o\
     $(src)/paraFMDF/particleRoutines/part_integrate_adbash.o\
     $(src)/paraFMDF/particleRoutines/part_interp_vis.o\
     $(src)/paraFMDF/particleOperations.o\
     $(src)/paraFMDF/particleRoutines/part_map.o\
     $(src)/paraFMDF/particleRoutines/part_find.o\
     $(src)/Particles/readpart.o\
     $(src)/paraFMDF/sourcing.o\
     $(src)/paraFMDF/compressibility.o\
     $(src)/paraFMDF/MPI/part_mpi.o\
     $(src)/Statistics/dongru/sendrecv_dom.o\
     $(src)/MPI/gatherbcast_dt.o\
     $(src)/Statistics/dongru/sendreceive_Q.o\
     $(src)/MPI/grid_partition.o\
     $(src)/MPI/preprocessing.o\
     $(src)/MPI/metis_file.o\
     $(src)/EV_routines/Main/Driver.o\
     $(src)/tstep.o\
     $(src)/EV_routines/Main/Prolong_gtol.o\
     $(src)/EV_routines/Main/RHS.o\
     $(src)/EV_routines/Main/Ent_Fluxes.o\
     $(src)/EV_routines/Main/EViscosity_GL.o\
     $(src)/EV_routines/Main/Entropy_value.o\
     $(src)/EV_routines/Patching/mortar_operations_ent_RI.o\
     $(src)/EV_routines/Patching/send_flux_ent.o\
     $(src)/EV_routines/Filter/viscosity_smoothing.o\
     $(src)/Energydissipation/energydissipation.o\


default:$(OBJS)
	$(LD) -I $(FFTW)  $(LDFLAGS) -o $(PROG) *.o $(FFTWL)
	rm -rf *.mod *.o

clean :
	rm -r *.o *.mod
