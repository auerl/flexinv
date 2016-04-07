#Flexinv

```fortran
  	          _|_|_|_|  _|        _|_|_|_|  _|      _|  _|_|_|  _|      _|  _|      _|  
	          _|        _|        _|          _|  _|      _|    _|_|    _|  _|      _|  
	          _|_|_|    _|        _|_|_|        _|        _|    _|  _|  _|  _|      _|  
	          _|        _|        _|          _|  _|      _|    _|    _|_|    _|  _|    
	          _|        _|_|_|_|  _|_|_|_|  _|      _|  _|_|_|  _|      _|      _|   
```

Flexinv is a Fortran/Python toolbox to create global, adaptive resolution transversely isotropic tomographic models of the entire mantle, via joint inversion of surface-wave dispersion and body-wave traveltimes, in the high-frequency ray approximation. See _Auer et al. (2014)_ and _Boschi (2009)_ for a complete description of the algorithm.

Some components of the surface wave parts of this code go back to Woodhouse _(1981)_ and _Dziewonski & Anderson (1981)_. The body-wave routines go back to _Gu (2005)_. Crustal corrections are based on CRUST2.0 by _Laske et al. (2002)_. Models are parameterized in curvilinear hexahedrons, whose size is adapted to local ray coverage, following _Schäfer et al. (2011)_.

Please note that this README file is not exhaustive and does not cover all intriciacies of hidden in the code, but should, at least, provide you a means to get started

## Installing Flexinv

### Requirements

The Fortran core of this software is self contained and only requires a recent Fortran compiler. For the python based wrapper scripts, we recommend to use the [Anaconda](https://www.continuum.io/downloads) Python distribution. All required libraries should be contained therein.

### Installation

After obtaining the code with

```bash
git clone http://github.com/auerl/flexinv.git .
```

it can be compiled with any modern compiler. We have thoroughly tested our code with the 2014 and 2015 versions of the Intel compiler collection on Ubuntu 12.04 and 14.04, but any UNIXoid OS should be fine. Build the code with

```bash
cd flexinv
./build.sh
```

## Running Flexinv

Setting up a new global tomographic imaging problem, requires reasoning about a myriad of parameters, and only a quick overview, leading approximately to _Model A_ from _(Auer et al. 2014), shall be given here. Everyone who aims addressing a specific problems is strongly encourage to get in touch with me.

### Preliminaries

While the core of flexinv is written in Fortran, we provide a python wrapper for convenience, which can be via handled via command line arguments and/or a simple parameter file `inparam` in the main program folder. First of all, you should replace all references to `/path/to/flexinv` with your actual installation path in in `inparam`. Also chose folders where to store matrices and inversion results in the `[folders]` category of `inparam`. 

### Parameter selection
First we need to chose lateral and vertical parameterization, as well as a set of physical inversion parameters. Flexinv currently supports orthogonal, curvilinear hexahedral basis functions or _voxels_, with adapted azimuthal increments at the poles to approximately maintain a uniform volume over the Earth's sphere. Since, for this test problem we focus on S-waves, we chose to invert for transversely isotropic shear wavespeeds v<sub>SH</sub> and v<sub>SV</sub> (Auer et al. 2014), at 28 layers of variable thickness from 0 to 2891 km depth and an equatorial increment of 5.0°.

### Surface-wave matrices
The module `mat/suwa` solves forward problem for surface waves, relating frequency dependent phase anomalies with wavespeed perturbations and assembles surface wave submatrices in the CSR (compressed-spares-row) Yale-type binary format. Using the

```bash
./setup.py --suwa_matrix --sw_data_id sw.fm.ac --eqincr 5.0 --inv_pars '3 4' --layers l28
```

The argument `--suwa_matrix` tells setup.py to prepare a series of submit scripts for surface wave submatrices, split up by frequency and overtone number. As you see, the desired physical parameters v<sub>SH</sub> and v<sub>SV</sub> have been passed through the `--inv_pars '3 4'`. Inverting for v<sub>PH</sub>, v<sub>PV</sub>, v<sub>SH</sub> and v<sub>SV</sub> would result in `--inv_pars '1 2 3 4'`. The argument after `--layers` specifies the desired vertical parameterization (see `inparam` to clarify which exact parameterization the _l28_ corresponds to), while `--eqincr 5.0` tells setup.py to use a five-degree lateral mesh size and `--sw_data_id` specifies the ID of the data set to be used (see inparam, to figure out which ID corresponds to which data set).

The script will create a folder to store matrices, and submit scripts for the LSF job schedulers at HPC facilities. Please adapt the code whenever you wish to submit to different schedulers such as SLURM. Submit all jobs to the scheduler by 
```bash
sh /path/to/flexinv/out/matrices/sw.fm.ac-28-p24-5.0/submit.sh
```

Make sure that enough disk space is available, since surface wave matrices may occupy up to a couple hundres of GB of disk space.

surface wave kernels ...

### Body-wave matrices 
The module `mat/bowa` solves the forward problem for body waves, relating cross-correlation based or picked traveltime delays with wavespeed perturbations, and assembles body-wave submatrices. Command line arguments are the same as for surface waves.

```bash
./setup.py --bowa_matrix --bw_data_id rit --eqincr 5.0 --inv_pars '3 4' --layers l28
```
Note that computation of sensitivity kernels as well as assembly of matrices for body waves is much more efficient than for surface waves. Hence, we chose to do it on the fly, directly on the interactive (login) nodes.

### Inversion

Finally, setup.py provides functionality to combine stored matrices via a specific weighting scheme into a text file, representative of an inversion 'schedule'. These schedule files can be read, e.g. with the module `inv`, which is a simple, single-core implementation of the LSQR algorithm of _Paige & Saunders (1982)_. We strongly recommend to use our new parallel tomography solver [PETScinv](../petscinv), which is based on the excellent [PETSc library](https://github.com/petsc/petsc). Create an inversion schedule by entering a command similar to

```bash
./setup.py --inversion --inv_id 'inv1' --inv_dirs S /path/to/flexinv/out/matrices/sw.fm.ac-p34-28-5.0 1.0 B /path/to/flexinv/out/matrices/bws-p34-28-5.0 .0 5.0 --layers l28 --eqincr 5.0 --ddmp_id dd28eq --rdmp_id dd28eq --weight_id ws1
```

Generally you
need to upweight the body wave portion of the data due to their lower RMS.
Initially it is okay to use a rather rough weighting, which can be refined 
later on. 

The argument `--inversion` tells setup.py to activate the inversion major mode, the `--inv_id` defines how the output folder will be titled, the list of folders after `--inv_dirs` are automatically crawled for pre-computed submatrices. The _S/B_ defines, whether the folder contains body or surface wave matrices, and the real value following the folder path is an additional weighting factor given to the matrices stored in that folder. The `--ddmp_id`,`--rdmp_id` and `--weight_id` define theID's of the anisotropy and roughness damping schemes, respectively (see inparam for details).


### Multi-scale meshing

... To be updated ... 

![Adaptive mesh](https://cloud.githubusercontent.com/assets/5452469/14347307/b81e88c4-fcb6-11e5-9b63-959849e71e6c.png =300x) 

--project
--proj_dirs
--adaptive
--inv_pars '1 2 3 4'

### Plotting and analysis

Plotting and analysis of models can be performed using our [pytomo](../pytomo) library.
