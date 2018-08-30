# DSID 
This project includes the implementation of Differential Stress Induced Damage (DSID) model and its application for engineering problems based on our paper: 

Jin, W., Xu, H., Arson, C., & Busetti, S. (2017). Computational model coupling mode II discrete fracture propagation with continuum damage zone evolution. International Journal for Numerical and Analytical Methods in Geomechanics, 41(2), 223-250.

Please kindly cite above paper if you used any of the functions or algorithms listed in this Github repository, thank you


Matlab code:

DSID_CP.m -> Matlab implementation of DSID model at Gauss Point, cutting plane method (return mapping) is used for iteration;

DSID_direct_iteration.m -> Matlab implementation of DSID model at Gauss Point, direct iteration is used;

Kachanov.m -> Calculate the damage effective modulus with crack inside a plate


UMAT

UMAT_DSID_CP_3D.for/UMAT_DSID_CP_plain_strain.for -> Abaqus UMAT implementation of DSID model using cutting plane algorithm for 3 dimentiona/plane strain cases


Abaqus Input

input files used for the above mentioned paper


