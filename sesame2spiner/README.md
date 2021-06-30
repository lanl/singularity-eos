# sesame2spiner

Self-contained tool that uses `EOSPAC` to read `SESAME` tables and convert to Spiner HDF format.

## Download

```bash
git clone --recursive git@gitlab.lanl.gov:jonahm/sesame2spiner.git
```

## Dependencies

`sesame2spiner` requires the following header-only libraries, which are provided via `git` submodules:
- [Spiner](https://gitlab.lanl.gov/jonahm/spiner)
- [Catch2](https://github.com/catchorg/Catch2)
- [nlohmann/json](https://github.com/nlohmann/json)

It also requires the following pre-installed libraries
- [EOSPAC](https://laws.lanl.gov/projects/data/eos/eospacReleases.php)
- [HDF5](https://portal.hdfgroup.org/display/support)

## Build and install

`sesame2spiner` uses a simple Makefile. Cloning and then typing `make`
will build the binary and run the tests.

You may need to change include libraries. In particular, you may need
to set the following flags:
- `EOSPAC_INCLUDE`, the `-I` flag for EOSPAC headers
- `EOSPAC_LIB`, the `-L` flag for EOSPAC 
- `EOSPAC_LINK`, the `-l` flag for EOSPAC. In particular, you may want to set the `rpath`.
- `CC`, the compiler. You likely want this to be the path to `h5c++`
  for `hdf5-serial`. You can also change it to, say, `icc` and simply
  point the appropriate include and link flags `INCLUDE_FLAGS`,
  `LFLAGS` to your `hdf5` installation.
  
Currently there is no `make install`. For now, please just copy the
binary to where you want it.

## Example usage

```bash
jonahm@sn-fey1:$ ./sesame2spiner -h
Usage: ./sesame2spiner[-p] [-h] <parameter file>

	 <parameter file>: input file in json format
	-p: print metadata associated with materials in parameter file
	-h: print this message

Example JSON file:

{
  "savename" : "materials.sp5",
  "materials" : [
    { // only matid is required. All others override defaults.
      "matid"  : 5030,
      "name"   : "air",
      "rhomin" : 1e-2, // g/cc
      "rhomax" : 10,
      "Tmin"   : 252,  // kelvin
      "Tmax"   : 1e4,
      "siemin" : 1e12, // erg
      "siemax" : 1e16
    },
    {
      "matid" : 4272,
      /* These shrink logarithm of bounds
         by a fraction of the total interval <= 1.
       */
      "shrinklRhoBounds" : 0.15,
      "shrinklTBounds" : 0.15,
      "shrinkleBounds" : 0.5
    }
  ]
}



jonahm@sn-fey1:$ ./sesame2spiner -p examples/air_and_steel.json 
sesame2spiner                            
-----------------------------------------
Author: Jonah Miller (jonahm@lanl.gov)   
-----------------------------------------

Saving to file materials.sp5
Processing 2 materials...
...5030
MATID: 5030
	name: dry air
	exchange coefficient = 0
	mean atomic mass     = 14.803
	solid bulk modulus   = 0
	normal density       = 0.001293
	[rho min, rho max]   = [1e-07, 15]
	[T min, T max]       = [175.235, 3.4815e+08]
	[sie min, sie max]   = [8.402e+08, 2.484e+16]
	num rho                = 21
	num T                  = 31
	Comments:
101: material. dry air (z=7.37296, a=14.80304)/source. h. c. graboske/date. dec 81/refs. ucid-16901/comp. n2 (0.7809), o2 (0.2195), ar (0.0096)/codes. see ucid-16901 /Classification. Unclassified /                                               

...4272
MATID: 4272
	name: stainless steel 347
	exchange coefficient = 0
	mean atomic mass     = 55.3678
	solid bulk modulus   = 0
	normal density       = 7.91
	[rho min, rho max]   = [0, 39550]
	[T min, T max]       = [0, 1.16048e+09]
	[sie min, sie max]   = [-8.22568e+08, 7.08239e+16]
	num rho                = 110
	num T                  = 79
	Comments:
101: material. stainless steel 347 (z= 25.80 a= 55.37) /source. j. c. boettger /date. nov 01 /refs. none /comp. weight % : fe 70, cr 19, ni 11 /codes. grizzly /Classification. Unclassified /                                                      
102: mixed tfd used for the electronic table and high compression cold curve. jdjnuc used for nuclear table with igrun=7, gamref=1.66, dgamma=-0.05, and tmelt=1986 kelvin. the cold curve was obtained with the chug model with cmat=2.0, faclj=0.7, ecohkc=99.4, and us=4.62+1.42up.                                              

All done. Cleaning up.
jonahm@sn-fey1:$ h5ls materials.sp5 
4272                     Group
5030                     Group
air                      Soft Link {5030}
stainless\ steel\ 347    Soft Link {4272}
```

## Copyright

© 2021. Triad National Security, LLC. All rights reserved.  This
program was produced under U.S. Government contract 89233218CNA000001
for Los Alamos National Laboratory (LANL), which is operated by Triad
National Security, LLC for the U.S.  Department of Energy/National
Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of
Energy/National Nuclear Security Administration. The Government is
granted for itself and others acting on its behalf a nonexclusive,
paid-up, irrevocable worldwide license in this material to reproduce,
prepare derivative works, distribute copies to the public, perform
publicly and display publicly, and to permit others to do so.
