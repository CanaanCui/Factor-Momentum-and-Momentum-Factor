# Factor-Momentum-and-Momentum-Factor

Code Files:

1) FM_MainCode.do

This code replicates almost all of the tables and figures in the paper. Except for one test that requires the original CRSP file (taken from WRDS), all other input files necessary for running this code are provided.

When running this file, please ÔcdÕ (change the path) into the data directory Ñ the code assumes that the input files are located here, and it saves some temporary files into a subdirectory called Temp.
 
2) FM_ConstructKNSFactors.do

This code constructs the momentum-neutral factors and runs some tests based on those factors. The required additional input files are the monthly and daily CRSP files and the Kozak et al. file with characteristics, available on Serhiy KozakÕs website at https://sites.google.com/site/serhiykozak/data (address valid as of April 2022).

3) FM_ConstructOtherFormsOfMomentum.do

This file constructs monthly Industry-adjusted, Sharpe ratio, Novy-MarxÕs intermediate, and industry momentum factors. The required additional file is the monthly CRSP file.

4) FM_ConstructResidualMomentum.do

This code constructs residual momentum strategies. The required additional input files are the monthly and daily CRSP files 
