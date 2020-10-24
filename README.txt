Code for Extracting the hydrophobic defects or lipid packing defects on a membrane monolayer surface

Author: Madhusmita Tripathy, MBU IISc Bangalore, madhu.cfl@gmail.com

This code identifies the lattice grid points that constitute the hydrophobic defects in a lipid membrane.

As input, one has to povide the box dimensions, the selection of lipid atoms from one leaflet and the reference atoms over the trajectory. In simulations, where the membrane CoM is not fixed and evolves over the trajectory, one has to provide the box centers over the trajectory as well. 

The output is the defect grid points over the trajectory. For checks, one can also output the lipid selection and referece atom data.

Here, we analyze a lipid membrane with 5 kinds of atoms, and the radius information is hardcoded, so as most of the other parameters. One can easily extend this to multiple atom kinds by using parameter files.

The input files for the atom selections and the reference atom input data have 4 fields: atoms type(based on which the radius is assigned) and the 3 spatial coordinates. The atom selection is chosen from the bilayer as the lipid atoms that lie above a specific z_cutoff. One can use VMD or similar tools to generate these input data files. Or, one can also modify the code to handle other file types. However, these 4 are the only required fields.


Cite: M. Tripathy, Subashini T., A. Srivastava, JCTC 2020
