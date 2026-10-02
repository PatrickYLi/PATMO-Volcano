![](https://bytebucket.org/tgrassi/planatmo/raw/bf879231524ee15ec7f8205fee61af593133a8ca/icon_small.png?token=d5bb9ae2a8489f7f4b2e71a1bb9e75b8ba0bb12d)

# Welcome to PATMO

PATMO is a flexible code aimed at modelling 1D planetary atmospheres including (photo)chemistry, molecular/eddy diffusion, and multi-frequency radiative transfer.  
It is partially based on [KROME](http://www.kromepackage.org/), as it produces optimized problem-dependent Fortran code by using a Python pre-processor.

### Getting Started 

For the maintained volcanic case, see
[Pinatubo 1991 inputs and outputs](tests/volcano_pinatubo_1991/README.md).
Generate with `./compile.sh volcano_pinatubo_1991`, then run `make` in `build`.
`./test_volcano` is the standalone optical pre-run; `./test` runs background
spin-up followed by volcanic chemistry. Input generation does not refresh old
simulation outputs.

Check the [Wiki](https://sites.google.com/sophia-atmochem-lab.org/patmo-user-guide) (Currently maintained by Patrick and is always under construction)

### Credits  
PATMO is developed and mantained by [Tommaso Grassi](http://starplan.dk/users/tommaso)  
Niels Bohr Insitute, Starplan Center, University of Copenhagen   
Contributors: E.Simoncini, A.Chiavassa, J.Ramsey, N.Vaytet, A.Popovas, T.Haugbølle, S.Bovino.
