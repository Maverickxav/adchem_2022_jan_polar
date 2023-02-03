This directory contains the molecular cluster dynamics plugin ClusterIn, and an example of coupling the plugin to an aerosol dynamics model.

Input files for the H<sub>2</sub>SO<sub>4</sub>-NH<sub>3</sub> test system (Besel et al., *J. Phys. Chem. A* 124, 5931–5943, 2020) have been pre-generated, and an example of generating the files for this or other set-ups is given by run_perl.sh.
Example properties for the aerosol distribution are here obtained by subroutines in acdc_aerosol_parameters.f90; when the plugin is coupled to an aerosol model, the properties come from the host model.

Compile and run run_clusterin_example.f90 to see how the plugin works:
```console
$ make
$ ./run
```

Detailed information can be found in the technical manual.


Addition to readme by Carlton.
This version has 2 data set for testing. One H2SO4-NH3 and one H2SO4-DMA. After certain modifications to the copy_cluster_chem.sh it takes the bases used in different chemistries and updates it in the run_perl script. But the user should still manually change the file names (dipole/cluster energies/inp_neutral.inp).
This should atleast generate the *equation and *system files.