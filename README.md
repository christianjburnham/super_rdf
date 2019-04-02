# super_rdf
Code for calculating the radial distribution function using a supercell approach

Compilation:

gfortran -O3 -ffree-form -ffast-math super_rdf.f -o super_rdf
(requires common block file super_rdf.common)


Input files:
   super_rdf.control // control file
   input.xyz // xyz input file. Name can be set in control file.

super_rdf.control:
   line1: name of input file (string)
   line2: max distance (real) number of bins (integer)
       for constructing rdf histogram
   line3: type_a (int),type_b (int)
       calculates rdf between pairs of type a b, where type_a,type_b = 1 (O) or 2 (H)

input.xyz:
   first line: number of atoms
   second line: 'ANGLE' (keyword), amag, bmag, cmag, alpha, beta, gamma
      where amag,bmag,cmag,alpha,beta,gamma are the unit cell vectors in Angstroms
   next natomns lines: atomic_name ("O" or "H"), x,y,z (Angstroms)
      assume an order O1,H1a,H1b,O2,H2a,H2b,etc.

output:
    Will print the radial distribution function to file rdf.dat

*******
NOTES
*******

Code is intended to illustrate supercell method for summing over ij interactions in periodic boundary conditions where the
sum includes all interactions inside a cut-off sphere that can be of arbitrary radius; i.e. not restricted to fit inside
one unit cell.

The code, as is, is set up to calculate the radial distribution functions for water molecule input only. However, it
should be possible to modify to handle other systems.
