
# geomorph

Geomorph is a grid conversion tool for popular geodynamics and seismic
tomography data. The geomorph grid structure is based on an icosahedron, and
has the same native grid-structure as the TERRA and MantleVis applications.

## Building
```
git clone https://github.com/sglim2/geomorph
cd geomorph/src
make ARCH=GNU
```
If compiling on OSX used ```make ARCH=OSX``` instead.


## Input/Output grids
GeoMorph can handle various input grid structures: TERRA_CC, TERRA_CV, MVIS, MITP, FILT, GYPSUMS, GYPSUMP.

GeoMorph output is limited to grids structure identical to its native format: MVIS, TERRA_CC, and TERRA_CV


 * MVIS: MantleVis format
 * TERRA_CV: TERRA Convection model
 * TERRA_CC: TERRA Circulation model
 * MITP: Based on the MIT-P08 data format
 * FILT: FILT data-format
 * GYPSUMS: GypSum S-wave model
 * GYPSUMP: GypSum P-wabe model

## Examples

You may want to obtain some example data
```
git clone https://bitbucket.org/sglim2/geomorph-example-data
cd geomorph-example-data
```
Then any of the following commands should work:

 * Exmaple 1: Convert MITP-MVIS using optimised-searching nearest-neighbour
```
geomorph --mt 256 --nt 16 --nd 10 --suffix 01 --infile MITP08.txt --outfile mvis001 --intype mitp --outtype mvis --interp nearest2
```
 * Exmaple 2: Convert FILT-MVIS using optimised-search neareast-neighbour
```
geomorph --mt 32 --nt 16 --nd 10 --suffix 01 --infile Filt --outfile mvis001 --intype filt --filtinstart 50 --filtinend 2850 --filtinnum 57 --outtype mvis --interp nearest2
```
 * Exmaple 3: Convert MITP-TERRA_CV using optimised-search neareast-neighbour
```
geomorph --mt 32 --nt 8 --nd 5 --suffix 01 --infile MITP08.txt --outfile c001 --intype mitp --outtype terra_cv --interp nearest2
```
 * Exmaple 4: Convert MITP-TERRA_CC using optimised-search neareast-neighbour
```
geomorph --mt 32 --nt 8 --nd 5 --suffix 01 --infile MITP08.txt --outfile c002 --intype mitp --outtype terra_cc --interp nearest2
```
 *  Exmaple 5: Convert MITP-MVIS using linear interpolation
```
geomorph --mt 16 --nt 8 --nd 10 --suffix 01 --mtin 32 --ntin 8 --ndin 10 --suffixin 01 --cmbin 0.54940 --infile mvis002 --outfile mvis002 --intype mvis --outtype mvis --interp linear
```
 * Exmaple 6: Convert TERRA_CV-MVIS using linear interpolation
```
geomorph --mt 32 --nt 8 --nd 10 --suffix 00 --mtin 32 --ntin 8 --ndin 10 --suffixin 01 --cmbin 0.54940 --infile c002 --outfile mvis001 --intype terra_cv --outtype mvis --interp linear
```
 *  Exmaple 7: Scale mt=32 to mt=64 TERRA_CV-TERRA_CV using linear intepolation
```
geomorph --mt 64 --nt 16 --nd 10 --suffix 00 --mtin 32 --ntin 8 --ndin 10 --suffixin 01 --cmbin 0.54940 --infile c002 --outfile c003 --intype terra_cv --outtype terra_cv --interp linear
```
 * Example 8: Convert GYPSUMP-MVIS using brute-force nearest-neighbour
```
geomorph --mt 64 --nt 16 --nd 10 --suffix 01 --gypsuminnum 22 --gypsumlatloninfile Grid.LatsLons.1deg.txt --gypsumdepthinfile Grid.dpths.txt --infile Grid.GyPSuM.P --intype gypsump --outfile mvis001 --outtype mvis --interp nearest --cmbin 2775.0
```


