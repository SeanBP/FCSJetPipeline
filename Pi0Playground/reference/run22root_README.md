Instruction on fit for each tower invariant mass

0) Make sure all the jobs are completed successfully.

1) Merge all the root files: 

	./hadd_job.sh

2) Put the final root file ({$run}/StFcsPi0invariantmassAll.root) under run22root/ folder, you need to rename the name so as they are not overwrite each other.

	Merge all the root files for 2 runs into one final root file. The final root file must name under: StFcsPi0invariantmass{$day}testAll_iteration{$N}.root  where {$day} can be run day of the run , {$N} is number of iterations. Same as other step.


3) Use the command below to fit for all tower invariant mass plot and get fit status: (Please keep the single and double quote symbols, and keep the variable that you set in step 1)

	root4star -b -q draw_all_towersmassfit2.C'("{$day}","{$N}")' 

4) Check the tower fit status plot for south and north ECal (statusN_iteration{$N}.root, and statusS_iteration{$N}.root) 

	root -l statusN_iteration{$N}.root 

	statusN->Draw("colz")

	root -l statusS_iteration{$N}.root 

	statusS->Draw("colz")

	Status 5: good invariant mass peak : distance between invariant mass peak from fit and 𝜋^0 peak less than 0.0135 GeV, 10% difference

	Status 4: distance between invariant mass peak from fit and 𝜋^0 peak 0.0135 GeV to 0.045 GeV. (not too far away, really need gain correction)

	Status 3: We can find invariant mass peak beyond region of status 4&5 define. If it has good invariant mass peak fitting around 𝜋^0 invariant mass, we still consider to calculate the gain correction

	Status 2:bad fit result: 𝑠𝑖𝑔𝑚𝑎>0.07 𝑜𝑟 𝑚𝑎𝑠𝑠<0.08 𝑜𝑟 𝑚𝑎𝑠𝑠>0.29 [𝐺𝑒𝑉]

	Status 1: entries < 150.     
	
	Status 0 (white): no entry.   

	The ideal case for a good fit at pi0 invariant mass (0.135 GeV) will give you status 5  

	The towers close to beam (column < 5 ; 10 < row < 20) might have status 2 or 3 due to the radiation damage

	
5) You also need to check and save the invariant mass plot:

	root -l StFcsPi0invariantmass{$day}testAll_iteration{$N}.root

	h1_inv_mass_cluster->Draw()

Then you can see the invariant mass plot.

6) Put the gain correction files from calculation to the main working folder

	cp fcsgaincorr_{$day}_{$N}.txt ../fcsgaincorr.txt

7) One iteration is completed, you can back to the main working folder

