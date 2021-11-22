# Polarization-from-BEC
This document provides codes that may be useful for estimation of spontaneous polarization using Born-effective charges. 
The program reads the following input files and calculates the spontaneous polarization.

Input files

	polar: atomic positions of the polar phase (in fractional coordinates)
	non-polar: atomic positions of the non-polar phase (in fractional coordinates)
	born: Born-effective charge tensor in the same order of elements as in POSCAR
	coord: lattice parameter matrix

To run the code, gfortran is required. The spontaneous polarization is calculated in μC/cm^2.
