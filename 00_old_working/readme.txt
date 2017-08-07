****** The usage of the buckling_2D code ******

Author: Ziyi Zhu
Contact: wazzytrent@gmail.com
Last_updated: March 6th, 2017
***********************************************

Prerequisite:
	module load icc

To Compile:
	make                   -> DLVO version
  make srfrpl            -> srfrpl version
  make p_srfrpl          -> srfrpl parallel version

To Delete the compiled file:
	make clean

To Run the compiled file:
	b2 input output.dat

The sample input file and job submission scripts have been included
in the input and execute folder respectively.
