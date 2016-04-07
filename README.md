#Flexinv

Flexinv is a Fortran/Python toolbox to create global, adaptive resolution transversely isotropic tomographic models of the entire mantle, via joint inversion of surface-wave dispersion and body-wave traveltimes, in the high-frequency ray approximation. 

Some components of the surface wave parts of this code go back to Woodhouse (1981) and Dziewonski & Anderson (1981). The body-wave routines go back to Gu (2005). Crustal corrections are based on CRUST2.0 by Laske et al. (2002). Models are parameterized in curvilinear hexahedrons, whose size is adapted to local ray coverage, following Schäfer et al. (2011).

Please note that this README file is not exhaustive and does not cover all intriciacies of hidden in the code, but should, at least, provide you a means to get started

## Installing Flexinv

After downloading the code with

```bash
git clone http://github.com/auerl/flexinv.git .
```

it can be compiled with any modern compiler. We have thoroughly tested our code with the 2014 and 2015 versions of the Intel compiler collection on Ubuntu 12.04 and 14.04, but any UNIXoid OS should work just as well. Build the code with

```bash
cd flexinv
./build.sh
```

## Running Flexinv

### Preliminaries

### Chosing a dataset

### Parameter selection

### Surface-wave matrices

### Body-wave matrices 

### Multi-scale meshing

### Inversion

Please note that the legacy Fortran routines sometimes comprise hard-coded 



===================================================================
May 7th, 2014, Ludwig:
===================================================================
Slight changes to the code by Ludwig in May 2014:

I continued cleaning up the code and made it a bit more flexibel, i.e.
one can now just freely pick the inversion parameters with a new flag
--inv_pars '1 2 3 4' where 1=vph, 2=vpv, 3=vsh, 4=vsv. You should also
be able to rotate the parameters, e.g. '3 4 2', to invert for vsh,
vsv and vp -> vph sensitivity of surface waves seems pretty limited
so this might be a reasonable choice!

The code now handles up to 4 parameters. Also I generalized the style
the damping is applied: You have now two 'damping schemes' controlling
how much roughness and difference damping is applied to each parameter
and at each depth.

Furthermore I merged this version of the code with the one where I
added support for the new crustal model CRUST1.0. Just call setup.py
with --crust10 to compute the new kernels





===================================================================
March 5th, 2014:
===================================================================
Hi Andrea,

here is the current version of our software for joint inversion of 
body wave traveltimes and surface wave. I prepared everything in such 
a way that it should be more or less straightforward to embark on the
problem of Joint inversion of P and S.

To get started I would suggest you to start with the following tomography,
which should be doable on my workstation:

- Inversion for v_sh and v_sv only
- Fundamental modes and direct S-Phases only
- 28 layers through the whole mantle
- 5°x5° lateral block size everywhere
- Fundamental modes of Ekström 
- Direct S-phases from Ritsema

If that works, you can proceed with the Joint P-S tomography.

Here is the steps that you need to do:

1.) Compute localized surface wave kernels for the crustal model CRUST2.0
*** THIS STEP CAN BE SKIPPED SINCE I ALREADY PRECOMPUTED THE KERNELS FOR YOU ***

COMMAND:
./setup.py --suwa_matrix --sw_data_id etl.ac --eqincr 5.0 --layers l28

--suwa_matrix tells setup.py to prepare everything to compute surface wave matrices
--layers specifies the vertical parameterization, l28 is an identifier which is 
associated with a certain layering. See 'inparam', where the actual layer depths
for the id 'l28' is defined. You can easily change it.
--eqincr is the equatorial increment in degrees, i.e. lateral block size
--sw_data_id is the id of the surface waves dataset to be used, again see 'inparam'
to check which surface wave dataset is meant with that. etl stands for the ekström
dataset .ac means it has been corrected for azimuthal anisotropy beforehand

The software is prepared to automatically check whether kernels for a particular
choice of vertical paramterization already exist, or not. What you do is basically
to call setup.py with certain command line arguments, telling it to compute surface
wave matrices. In case no kernels are found in the folder __kernels__ it will compute
them (actually it will prepare a bunch of shell scripts, which you can call sub-
sequently in parallel to compute the kernels). You can easily verify whether the
kernel database is complete: Just check if there is 16200 unique files called
L_??????_000_kernels in the folder c2-av-28/c2-28. The av indicates that the kernels
are parameterized in voigt average and anisotropy kernels. The 28 means 28 layers.
c2 means they are CRUST2.0 kernels. In the kernels files (ascii) there should be
28 lines for each mode.


2.) Compute the surface wave matrices

When you read 1) you will notice that the call to compute kernels and matrices is
actually the same

COMMAND:
./setup.py --suwa_matrix --sw_data_id etl.ac --eqincr 5.0 --layers l28

Navigate to the folder __matrices__/etl.ac-28-5.0 .. the name of the folder should
be self-explanatory. In the file you find a bunch of bash scripts, and one script
called submit.sh, which is meant to be run on a cluster with bsub(lsf) queing 
system. Since you want to run the thing on a normal workstation without that, please
just remove the 'bsub -R lustre -W35:59 < ' infront of the path and add a & behind
each .sh to run it in the background. After every 6th line add a new line with wait
such that it runs only 6 of the scripts at once (i only have 8 processors on my
machine). Sorry... you can of course change the setup.py to do all of that.

BEFORE:
bsub -R lustre -W35:59 < /home/andrea/flexinv2/__matrices__/etl.ac-28-5.0/sw.000.sh

AFTER:
/home/andrea/flexinv2/__matrices__/etl.ac-28-5.0/sw.000.sh &
/home/andrea/flexinv2/__matrices__/etl.ac-28-5.0/sw.001.sh &
/home/andrea/flexinv2/__matrices__/etl.ac-28-5.0/sw.002.sh &
wait
/home/andrea/flexinv2/__matrices__/etl.ac-28-5.0/sw.003.sh &
/home/andrea/flexinv2/__matrices__/etl.ac-28-5.0/sw.004.sh &
...

Then just type ./submit.sh to compute your matrix. This will involve a call to the
programm matrix_sw_vx which is in flexinv2/mat/suwa/, and needs all the 1D kernel
profiles. Computing the matrix will take quite a while and the outcome consumes
quite a bit of disk space.

Look at the sw.???.log files to see what happens. What you will get is sw
matrices, ready to be read by the LSQR based inversion code, in the Yale
CSR sparse matrix format (values,indices,pointers plus the rhs).


3.) Compute the body wave matrices

Now we compute the body wave matrices. Contrary to the surface wave matrices we
do not need to precompute the kernels, since crustal corrections for body wave
traveltime delays are computed in a very different fashion and our sensitivites
dont depend on the crust, so much.

COMMAND:
./setup.py --bowa_matrix --bw_data_id rit --eqincr 5.0 --layers l28

Since everything is a bit quicker for body waves i didn't take the effort to
parallelize this part of the code. The computation of the body wave matrices
just starts immediately. Do not worry about occasional (or sometimes not so
occasional) error messages. Mostly this just means that certain paths could
not be used due to a weird geometry.

What you will get is a folder rit-28-5.0 in __matrices__ containing matrices
in the same format as the surface wave matrices.


4.) Inverting different submatrices together

Finally we want to combine all our precomputed matrices in a certain way (i.e.
using a specific weighting scheme) and invert them together. Generally you
need to upweight the body wave portion of the data due to their lower RMS.
Initially it is okay to use a rather rough weighting, which can be refined 
later on. 

COMMAND
./setup.py --inversion --inv_id 'Test1' --inv_dirs S /home/andrea/flexinv2/__matrices__/etl.ac-28-5.0 1.0 B /home/andrea/flexinv2/__matrices__/rit-28-5.0 5.0 --layers l28 --eqincr 5.0 --dmp_id dfilt --weight_id ws1

Okay, what is happening here: --inversion tells setup.py to setup the 
inversion scripts. --inv_id is how the result folder will be named, --inv_dirs
is followed by a list of folders which are supposed to be crawled for usable matrices.
The Structure is always S/B to tell him whether it is body or surface waves, then
the absolute path to the folder, and then an manual additional weighting. Just use 1.0
for surface waves and 5.0 for body waves for now.You can make this as long as you want.
--dmp_id is the id of the damping schedule it defines how exactly vertical and horizontal 
roughness damping is tuned (it can be different in different layers). dfilt is a pretty 
simple damping schedule.
--weight_id is the weighting schedule, just use ws1. Like the damping schedule it is
defined in the file inparam.

You can then navigate to __inversions__/Test1.* in which, again, a bunch of scripts
for different damping parameters should have been created (pc.*.sh). Normally, you
would call them via the submit-script submit.sh, but you can also do chmod 755 pc*.sh
and call them individually. Their stdout will be written to a log file. Modify the
first line in the pc*sh scripts to just see the output on the command line.

5.) Postprocessing and plotting

Convert between vsh/vsv and voigt/xi by typing

postproc.py --convert

And plot your results with the command

comp_rg_vs.sh dvoigt.???.pcn layers.in 5.0

Just create your own copy of /work/ana/netcdfdrawmap/comp_rg_vs.sh and modify what
is in it, to make different plots.

Good luck!

:-)


