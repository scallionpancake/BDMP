****** The usage of the buckling_2D code ******

Author: Ziyi Zhu
Contact: wazzytrent@gmail.com
Last_updated: March 6th, 2017
***********************************************
UPDATED:
Steve Kuei
kuei.steve@rice.edu
080717

***********************************************
to be run on davinci.rice.edu

Prerequisite:
	module load icc

To Compile:
	make                   -> openMP, 2D version, named paramag.out

To Delete the compiled file:
	make clean

To Run the compiled file:
	paramag input output.dat

(note 'input.in' will fail - remove extension)

execute
- slurm batch files
input
- .in files with different input parameters.
output
- saved outputs. 
