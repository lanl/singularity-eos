## Build and install

`sesame2spiner` is built with ``singularity-eos`` 
if ``-DSINGULARITY_BUILD_SESAME2SPINER=ON`` is specified at the configure stage.

## Example usage

```bash
jonahm@sn-fey1:$ ./sesame2spiner -p examples/air.dat examples/steel.dat
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

© 2021-2022. Triad National Security, LLC. All rights reserved.  This
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
