# FCS-ECal-pi0-reconstruction

STAR FCS ECal pi0 reconstruction analysis framework

Author: Xilin Liang

Email: xilin.liang@email.ucr.edu


Instruction:

The framework contain:

	1) StRoot : the maker for pi0 reconstruction

	2) runMuDst.C : run macro, include the pi0 reconstruction by highest energy pair

	5) fcsgaincorr.txt : gain correction file (in GitHub page, only the gain correction with 1.2 value, which only for the first iteration of the first pair of run for each period!). Generally, the values are varied by different towers. And for each iteration, you will need to change this file before each iteration
	
	6) submitScheduler: the folder for submitting jobs to Scheduler
	
	7) run22root / run23root / run24root: the folder for saving root files and analysis

	8) Gain corrections predict scripts (current for run 22 only): include: integralLumi.txt , led_files2/ , ledRatio_time_V4.C , ledruns.txt , run_lumin.txt


How to do the analysis on pi0 reconstruction:
	
	1) Make sure the files and folders are created and accessable

	2) Complete the one-time set up on the scripts (No needed anymore)
	
	3) Compile the library: (follow the command below at the main working folder (current folder))
		 
		 stardev
		 
		 cons

	4) Start with the first set of runs, you will use the fcsgaincorr.txt for the first iteration. Use the fcsgaincorr.txt from GitHub page for the first pair of run (for each pair, choose the first run number, for example: choose 23074017 for 23074017/23074018). 

	5) Then start the iteration (repeat several times until you get the good invariant mass close to pi0 invariant mass):
		
		5.1) Go to submitScheduler/ folder, follow the instruction in submitScheduler/README.md

		5.2) Go to run22root/ folder, follow the instruction in run22root/README.md

	6) Save the final gain correction file fcsgaincorr_{$day}_{$N}.txt and fcsgaincorrOffline_{$day}_{$N}.txt . Akio will need the calibrated gain correction files. 

	7) Predict the initial gain correction files for the other pairs (do NOT do this for run 24 high lumin pp runs until further notice):

		7.1) Create a file call "physicsrunlist.txt" and write down the first run number of each pair.

		7.2) root4star -b -q ledRatio_time_V4.C'({$run},"{$gaincorrfile}")'

		,where {$run} is the first run number of the first pair; {$gaincorrfile} is the final gain correction file (fcsgaincorr_{$day}_{$N}.txt) from step 6) for the first pair.

		7.3) Then you will get fcsgaincorr_{$run}_ini.txt for the predicted initial gain correction files of each pair. 

	8) Use the predicted gain correction files as the input file for each pair (for run 24 high lumin pp runs, use the final gain correction file from first set of runs). Follow step 5) and 6) to do the invariant mass interation.
	
	Note: rename the fcsgaincorr_{$run}_ini.txt to fcsgaincorr.txt before you submit job!


