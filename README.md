# Reconstruction

Reconstruction software for the GRAND data.
Initial soft from TREND recons.

Software is based on the PORT library for function minimisation.
Execution is carried out through command:
```bash
bin/recons [RunId] 0
```
Input files are:
- coord_antennas.txt with format [AntennaId x(WE) y(SN) z(alt asl] where x,y, z are in meters
- R[RunId]_coinctable.txt with format [UnixTime AntId EvtId CoincId TrigTime x x x x x x x ]
TrigTime is given in nanoseconds wrt earliest trigger in coinc.

Output files are:
- -R[RunId]_planerecons.txt with format [CoincId UnixTime Mult Theta dTheta Phi dPhi Chi2 Signif]

Theta and Phi are given in GRAND conventions.
- -R[RunId]_sphrecons.txt with format [CoincId UnixTime Mult x y z t@source Chi2 Signif] 

where (x,y,z) are given in TREND coordinates.
To be done (04/12/2017): change code to write (x,y,z) in GRAND referential.

# Build

Requires gfortran. Then:
```bash
make clean && make
```

# License

The GRAND reconstruction package is under the **GNU LGPLv3** license. See the
provided [LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files.
