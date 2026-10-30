
const float peakmin=0.1;
const float peakmax=0.2;
const float tailmax=0.4;
double xexp(double* x, double* p, int flag=0) {
	if(x[0]<p[1]) return 0.0;
	return p[0]/p[2]/p[2]*pow(x[0]-p[1],p[3])*exp(-(x[0]-p[1])/p[2]);
}

double gaus(double* x, double* p, int flag=0) {
	return p[0]*exp(-0.5*pow((x[0]-p[1])/p[2],2));
}

double gxexp(double* x, double* p, int flag=0) {
	return xexp(x,p) + gaus(x,&p[4]);
}

//using Akio's fitting function
void draw_all_towersmassfit2(const char* day, const char* iteration)
{
	//TFile* File1005=TFile::Open(Form("StFcsPi0Result%s_All.root",date));
	TFile* File1005=TFile::Open(Form("StFcsPi0invariantmass%stestAll_iteration%s.root",day,iteration));

	TCanvas *c5[36];
	TCanvas *c6[36];

	TCanvas* c1=new TCanvas("c1","c1",22000,22000);
	TCanvas* c2=new TCanvas("c2","c2",22000,22000);
	TCanvas* c3=new TCanvas("c3","c3",800,800);
	TCanvas* c4=new TCanvas("c4","c4",800,800);
	c5[1]=new TCanvas("c5","c5",6000,6000);
	c6[1]=new TCanvas("c6","c6",6000,6000);
	FILE *tower_input = fopen("fcsdet.txt","r");
	//FILE *depchid_input = fopen("depch_id.txt","r");


	double peakmass = 0.135;
	int fitentries = 150;

	c1.Divide(22,34);
	c2.Divide(22,34);
	c5[1].Divide(5,5);
	c6[1].Divide(5,5);
	float mergin=0.1;
	int ns,id,row,col,ehp,dep,ch;
	Double_t oldgaincorr;
	int depchtoid[2][24][32];
	for (int i=0; i<24 ; i=i+1)
	{
		for (int j=0; j<32; j=j+1)
		{
			depchtoid[0][i][j]=-1;
			depchtoid[1][i][j]=-1;
		}
	}
	
	Double_t NorthTowerPeak[748],SouthTowerPeak[748];
	Double_t NorthGainCorr[748],SouthGainCorr[748];
	Double_t NorthGainCorrRatio[748],SouthGainCorrRatio[748];
	Double_t oldNorthGainCorr[748],oldSouthGainCorr[748];
	TH2I* statusN = new TH2I("statusN","North Ecal tower invariant mass plot status;column;row",22,1,23,34,-34,0);
	TH2I* statusS = new TH2I("statusS","South Ecal tower invariant mass plot status;column;row",22,-22,0,34,-34,0);
	TH2D* chi2ndfN = new TH2D("chi2ndfN","North Ecal tower fit chi2/ndf;column;row",22,1,23,34,-34,0);
	TH2D* chi2ndfS = new TH2D("chi2ndfS","South Ecal tower fit chi2/ndf;column;row",22,-22,0,34,-34,0);
	TH2D* h2_gaincorrectionN = new TH2D("h2_gaincorrectionN","North Ecal tower gain correction;column;row",22,1,23,34,-34,0);
	TH2D* h2_gaincorrectionS = new TH2D("h2_gaincorrectionS","South Ecal tower gain correction;column;row",22,-22,0,34,-34,0);
	TH2D* h2_masspeakN = new TH2D("h2_masspeakN","North Ecal tower Gaussian peak fit mass;column;row",22,1,23,34,-34,0);
	TH2D* h2_masspeakS = new TH2D("h2_masspeakS","South Ecal tower Gaussian peak fit mass;column;row",22,-22,0,34,-34,0);
	TH2D* h2_Northgaincorrration = new TH2D("h2_Northgaincorrration","North Ecal tower gain correction ratio;column;row",22,1,23,34,-34,0);
	TH2D* h2_Southgaincorrration = new TH2D("h2_Southgaincorrration","South Ecal tower gain correction ratio;column;row",22,-22,0,34,-34,0);
	TH1D* h1_Northoldgaincorr = new TH1D("h1_Northoldgaincorr","old gain correction for North towers",748,0,748); 		
	TH1D* h1_Northnewgaincorr = new TH1D("h1_Northnewgaincorr","new gain correction for North towers",748,0,748); 		
	TH1D* h1_Southoldgaincorr = new TH1D("h1_Southoldgaincorr","old gain correction for South towers",748,0,748); 		
	TH1D* h1_Southnewgaincorr = new TH1D("h1_Southnewgaincorr","new gain correction for South towers",748,0,748); 		
	TH1D* h1_Northgaincorrdist = new TH1D("h1_Northgaincorrdist","gain correction for North towers",80,0.9,1.7); 		
	TH1D* h1_Southgaincorrdist = new TH1D("h1_Southgaincorrdist","gain correction for South towers",80,0.9,1.7); 		
	TH1D* h1_Northgaincorrdist_used = new TH1D("h1_Northgaincorrdist_used","gain correction used for North towers",80,0.9,1.7); 		
	TH1D* h1_Southgaincorrdist_used = new TH1D("h1_Southgaincorrdist_used","gain correction used for South towers",80,0.9,1.7); 		
	TH1D* h1_Northgaincorrdist_status3 = new TH1D("h1_Northgaincorrdist_status3","gain correction for North towers (status3)",80,0.9,1.7); 		
	TH1D* h1_Southgaincorrdist_status3 = new TH1D("h1_Southgaincorrdist_status3","gain correction for South towers (status3)",80,0.9,1.7); 		
	TH1D* h1_Allgaincorrdist = new TH1D("h1_Allgaincorrdist","gain correction for All towers",80,0.9,1.7); 		
	
	ifstream depchid_input;
	depchid_input.open("depch_id.txt");
	ifstream inputGainCorr;
	inputGainCorr.open("../fcsgaincorr.txt");

	ofstream outputGainCorr;
	outputGainCorr.open(Form("fcsgaincorr_%s_%s.txt",day,iteration));
	ofstream outputGainCorrOffline;
	outputGainCorrOffline.open(Form("fcsgaincorrOffline_%s_%s.txt",day,iteration));
	ofstream outputGainCorrRatio;
	outputGainCorrRatio.open(Form("fcsgaincorrratio_%s_%s.txt",day,iteration));

	for (int i=0; i<(748*2);i=i+1)
	{

		depchid_input >> ehp >> ns >> dep >> ch >> id;
		depchtoid[ns][dep][ch] = id;
		inputGainCorr >> ehp >> ns >> dep >> ch >> oldgaincorr;
		if (ns == 0) {oldNorthGainCorr[id]=oldgaincorr;} else {oldSouthGainCorr[id]=oldgaincorr;}

	}


	for (int i=0; i<748; i=i+1)
	{

		cout<<"tower:"<<i<<endl;
		char title_hist[100];
		fscanf(tower_input,"%i %i %i %i", &ns, &id, &row, &col);
		sprintf(title_hist,"mass_by_Ntower_%i",i);
		TH1F* hist=(TH1F*)File1005->Get(title_hist);
		hist->SetTitle(Form("invariant mass sorted by North %i tower (row: %i , col: %i)",id,row,col));
		TVirtualPad *pad=c1->cd(i+1);
		pad->SetRightMargin(mergin); 
		pad->SetLeftMargin(mergin);
   		pad->SetTopMargin(mergin);   
		pad->SetBottomMargin(mergin);
		TVirtualPad *pad2=c5[row]->cd(col);
		pad2->SetRightMargin(mergin); 
		pad2->SetLeftMargin(mergin);
   		pad2->SetTopMargin(mergin);   
		pad2->SetBottomMargin(mergin);

		hist->SetLineWidth(1);
		hist->Draw("same");
		TLine* line=new TLine(0.135,0,0.135,(hist->GetMaximum()));
                line.SetLineColor(2);
		line->SetLineWidth(2);
		Int_t allentries;
		Int_t bkgtail;
		int status = 5;
		allentries=hist->GetEntries();
			cout<<"tower:"<<i<<endl;
		if (allentries < fitentries) {NorthTowerPeak[i]=peakmass;if (allentries == 0) {status = 0;} else {status = 1;NorthTowerPeak[i]=peakmass;}}
		else
		{
			Int_t max_bin=hist->GetMaximumBin();
			double max = hist->GetMaximum();
			Double_t max_pos = (max_bin-0.5)*0.008;
			if ((max_pos > 0.30)||(max_pos < 0.07)) max_pos = peakmass;

  			double peak = hist->GetBinCenter(hist->GetMaximumBin());
			//TH1* rebin = hist->Rebin(2,"newname");   //or 3
  			//double peak = rebin->GetBinCenter(rebin->GetMaximumBin());
  			if(peak<0.05 || peak>0.30) peak=peakmass;
  			printf("max=%f peak=%f\n",max,peak);

  			TF1* ftx = new TF1("xexp",xexp,0,tailmax,4); 
  			ftx->SetLineColor(kGreen);
  			ftx->SetParName(0,"A");     ftx->SetParameter(0,max/2);
  			ftx->SetParName(1,"x0");    ftx->FixParameter(1,0.0);
  			ftx->SetParName(2,"tau");   ftx->SetParameter(2,0.1);  ftx->SetParLimits(2,0.01,0.25);
  			//ftx->SetParName(2,"tau");   ftx->SetParameter(2,0.2); 
  			ftx->SetParName(3,"p");     ftx->SetParameter(3,1.0);
  			TFitResultPtr rest = hist->Fit("xexp","b","",0,tailmax);

  			TF1* f = new TF1("gxexp",gxexp,0,tailmax,7); 
  			f->SetParName(0,"A");     f->FixParameter(0,ftx->GetParameter(0));		 
  			f->SetParName(1,"x0");    f->FixParameter(1,ftx->GetParameter(1));		 
  			f->SetParName(2,"tau");   f->FixParameter(2,ftx->GetParameter(2));		 
  			f->SetParName(3,"p");     f->FixParameter(3,ftx->GetParameter(3));		 
  			f->SetParName(4,"peak");  f->FixParameter(4,max/2); f->SetParLimits(4,0.0,max);
  			f->SetParName(5,"mean");  f->FixParameter(5,peak);				 
  			f->SetParName(6,"sigma"); f->SetParameter(6,0.04); f->SetParLimits(6,0.005,0.1);
  			TFitResultPtr resg = hist->Fit("gxexp","","",0,tailmax);
			TF1* f2 = new TF1("gxexp2",gxexp,0,tailmax,7);
  			f2->SetLineColor(kRed);
  			f2->SetParameter(0,ftx->GetParameter(0));		 
  			f2->SetParameter(1,ftx->GetParameter(1));		 
  			f2->SetParameter(2,ftx->GetParameter(2)); f2->SetParLimits(2,0.01,0.25); 
  			f2->SetParameter(3,ftx->GetParameter(3));		 
  			f2->SetParameter(4,max/2); f2->SetParLimits(4,0.0,max);
  			f2->SetParName(5,"mean");f2->SetParameter(5,peak);				 
  			f2->SetParName(6,"sigma");f2->SetParameter(6,f->GetParameter(6)); f2->SetParLimits(6,0.002,0.1);
  			TFitResultPtr resg = hist->Fit("gxexp2","","",0,tailmax);
			f2->SetLineWidth(2);
			f2->Draw("same");
  			
  			TF1* bg = new TF1("bg",xexp,0,tailmax,4);
  			bg->SetLineColor(kBlue);
  			bg->SetParameter(0,f2->GetParameter(0));
  			bg->SetParameter(1,f2->GetParameter(1));
  			bg->SetParameter(2,f2->GetParameter(2));
  			bg->SetParameter(3,f2->GetParameter(3));
			bg->SetLineWidth(2);
			bg->Draw("same");

  			TF1* g = new TF1("gaus",gaus,0,tailmax,3);
  			g->SetLineColor(kMagenta);
  			g->SetParameter(0,f2->GetParameter(4));
  			g->SetParameter(1,f2->GetParameter(5));
  			g->SetParameter(2,f2->GetParameter(6));
			g->SetLineWidth(2);
			g->Draw("same");

			NorthTowerPeak[i]=f2->GetParameter(5);
			Double_t final_sigma=f2->GetParameter(6);
			Double_t final_chi2 = f2->GetChisquare();
			Double_t final_ndf = f2->GetNDF();
			Double_t final_chi2ndf = 1.0*final_chi2/final_ndf;
			chi2ndfN->Fill(col,-row,final_chi2ndf);

			if ((final_sigma > 0.07)||(NorthTowerPeak[i]<0.06)||(NorthTowerPeak[i]>0.29)) {status = 2;NorthTowerPeak[i]=peakmass;}
			else
			{
				if ((max_bin < 13)||(max_bin > 21)) {status = 3;}
				else
				{if ((NorthTowerPeak[i] < 0.1215)||(NorthTowerPeak[i]>0.1485)){status=4;} else {status = 5;}}
			}

		}
		statusN->Fill(col,-row,status);
		cout<<"Peak:"<<NorthTowerPeak[i]<<endl;
		NorthGainCorr[i] = oldNorthGainCorr[i]*(0.135 / NorthTowerPeak[i]);
		NorthGainCorrRatio[i]=0.135 / NorthTowerPeak[i];
		cout<<"North tower gain corr:"<<NorthGainCorr[i]<<endl;
		h2_gaincorrectionN->Fill(col,-row,NorthGainCorr[i]);
		if (status >= 2){h2_masspeakN->Fill(col,-row,NorthTowerPeak[i]);}
		if (status==5) { h1_Northoldgaincorr->Fill(i,oldNorthGainCorr[i]); h1_Northnewgaincorr->Fill(i,NorthGainCorr[i]);h1_Northgaincorrdist->Fill(NorthGainCorr[i]);h1_Northgaincorrdist_used->Fill(oldNorthGainCorr[i]);}
		if (status==3) { h1_Northgaincorrdist_status3->Fill(NorthGainCorr[i]);}
		h1_Allgaincorrdist->Fill(NorthGainCorr[i]);
		h2_Northgaincorrration->Fill(col,-row,NorthGainCorrRatio[i]);
		outputGainCorrOffline << ns <<"    " <<i<<"     "<<NorthGainCorr[i]<<endl;
		outputGainCorrRatio << ns <<"    " <<i<<"     "<<NorthGainCorr[i]<<"      "<<status<<endl;
		if (col==22)
		{

			if (row ==1) {c5[row]->SaveAs(Form("NorthTowerResult_iteration%s.pdf[",iteration));}
			if (row ==34) {c5[row]->SaveAs(Form("NorthTowerResult_iteration%s.pdf]",iteration));}
			if (!((row ==34)||(row ==1))) {c5[row]->SaveAs(Form("NorthTowerResult_iteration%s.pdf",iteration));}
			c5[row+1]=new TCanvas("c5","c5",6000,6000);
			c5[row+1].Divide(5,5);

		}

	}	
	c1->SaveAs("mass_by_tower_Ncluster.jpg");
	statusN->SaveAs(Form("statusN_iteration%s.root",iteration));
	h2_gaincorrectionN->SaveAs("gaincorrectionN.root");
	h2_masspeakN->SaveAs("masspeakN.root");
	h2_Northgaincorrration->SaveAs("NorthGainCorrectionRatio.root");
	chi2ndfN->SaveAs("Chi2NDFN.root");

	for (int i=0; i<748; i=i+1)
	{

		cout<<"tower:"<<i<<endl;
		char title_hist[100];
		fscanf(tower_input,"%i %i %i %i", &ns, &id, &row, &col);
		sprintf(title_hist,"mass_by_Stower_%i",i);
//		sprintf(title_hist,"mass_point_by_tower_%i",i);
		TH1F* hist2=(TH1F*)File1005->Get(title_hist);
		//cout<<hist->GetBinContent(1)<<endl;
		hist2->SetTitle(Form("invariant mass sorted by South %i tower (row: %i , col: %i)",id,row,col));
		TVirtualPad *pad=c2->cd(i+1);
		pad->SetRightMargin(mergin); 
		pad->SetLeftMargin(mergin);
   		pad->SetTopMargin(mergin);   
		pad->SetBottomMargin(mergin);
		TVirtualPad *pad2=c6[row]->cd(col);
		pad2->SetRightMargin(mergin); 
		pad2->SetLeftMargin(mergin);
   		pad2->SetTopMargin(mergin);   
		pad2->SetBottomMargin(mergin);

		hist2->SetLineWidth(2);
		hist2->Draw("same");
		TLine* line=new TLine(0.135,0,0.135,(hist2->GetMaximum()));
                line.SetLineColor(2);
		line->SetLineWidth(2);
		Int_t allentries;
		Int_t bkgtail;
		int status = 5;
		allentries=hist2->GetEntries();
		if (allentries < fitentries) {SouthTowerPeak[i]=peakmass;if (allentries == 0) {status = 0;} else {status = 1;SouthTowerPeak[i]=peakmass;}}
		else
		{
			Int_t max_bin=hist2->GetMaximumBin();
			double max = hist2->GetMaximum();
			Double_t max_pos = (max_bin-0.5)*0.008;
			if ((max_pos > 0.29)||(max_pos < 0.07)) max_pos = peakmass;

  			double peak = hist2->GetBinCenter(hist2->GetMaximumBin());
			//TH1* rebin = hist2->Rebin(2,"newname");   //or 3
  			//double peak = rebin->GetBinCenter(rebin->GetMaximumBin());
  			if(peak<0.05 || peak>0.29) peak=peakmass;
  			printf("max=%f peak=%f\n",max,peak);

  			TF1* ftx = new TF1("xexp",xexp,0,tailmax,4); 
  			ftx->SetLineColor(kGreen);
  			ftx->SetParName(0,"A");     ftx->SetParameter(0,max/2);
  			ftx->SetParName(1,"x0");    ftx->FixParameter(1,0.0);
  			ftx->SetParName(2,"tau");   ftx->SetParameter(2,0.1);  ftx->SetParLimits(2,0.01,0.25);
  			//ftx->SetParName(2,"tau");   ftx->SetParameter(2,0.2); 
  			ftx->SetParName(3,"p");     ftx->SetParameter(3,1.0);
  			TFitResultPtr rest = hist2->Fit("xexp","b","",0,tailmax);

  			TF1* f = new TF1("gxexp",gxexp,0,tailmax,7); 
  			f->SetParName(0,"A");     f->FixParameter(0,ftx->GetParameter(0));		 
  			f->SetParName(1,"x0");    f->FixParameter(1,ftx->GetParameter(1));		 
  			f->SetParName(2,"tau");   f->FixParameter(2,ftx->GetParameter(2));		 
  			f->SetParName(3,"p");     f->FixParameter(3,ftx->GetParameter(3));		 
  			f->SetParName(4,"peak");  f->FixParameter(4,max/2); f->SetParLimits(4,0.0,max);
  			f->SetParName(5,"mean");  f->FixParameter(5,peak);				 
  			f->SetParName(6,"sigma"); f->SetParameter(6,0.04); f->SetParLimits(6,0.005,0.1);
  			TFitResultPtr resg = hist2->Fit("gxexp","","",0,tailmax);
			TF1* f2 = new TF1("gxexp2",gxexp,0,tailmax,7);
  			f2->SetLineColor(kRed);
  			f2->SetParameter(0,ftx->GetParameter(0));		 
  			f2->SetParameter(1,ftx->GetParameter(1));		 
  			f2->SetParameter(2,ftx->GetParameter(2)); f2->SetParLimits(2,0.01,0.25); 
  			f2->SetParameter(3,ftx->GetParameter(3));		 
  			f2->SetParameter(4,max/2); f2->SetParLimits(4,0.0,max);
  			f2->SetParName(5,"mean");f2->SetParameter(5,peak);				 
  			f2->SetParName(6,"sigma");f2->SetParameter(6,f->GetParameter(6)); f2->SetParLimits(6,0.002,0.1);
  			TFitResultPtr resg = hist2->Fit("gxexp2","","",0,tailmax);
			f2->SetLineWidth(2);
			f2->Draw("same");
  			
  			TF1* bg = new TF1("bg",xexp,0,tailmax,4);
  			bg->SetLineColor(kBlue);
  			bg->SetParameter(0,f2->GetParameter(0));
  			bg->SetParameter(1,f2->GetParameter(1));
  			bg->SetParameter(2,f2->GetParameter(2));
  			bg->SetParameter(3,f2->GetParameter(3));
			bg->SetLineWidth(2);
			bg->Draw("same");

  			TF1* g = new TF1("gaus",gaus,0,tailmax,3);
  			g->SetLineColor(kMagenta);
  			g->SetParameter(0,f2->GetParameter(4));
  			g->SetParameter(1,f2->GetParameter(5));
  			g->SetParameter(2,f2->GetParameter(6));
			g->SetLineWidth(2);
			g->Draw("same");

			SouthTowerPeak[i]=f2->GetParameter(5);
			Double_t final_sigma=f2->GetParameter(6);
			Double_t final_chi2 = f2->GetChisquare();
			Double_t final_ndf = f2->GetNDF();
			Double_t final_chi2ndf = 1.0*final_chi2/final_ndf;
			chi2ndfS->Fill(-col,-row,final_chi2ndf);

			if ((final_sigma > 0.07)||(SouthTowerPeak[i]<0.06)||(SouthTowerPeak[i]>0.29)) {status = 2;SouthTowerPeak[i]=peakmass;}
			else
			{
				if ((max_bin < 13)||(max_bin > 21)) {status = 3;}
				else
				{if ((SouthTowerPeak[i] < 0.1215)||(SouthTowerPeak[i]>0.1485)){status=4;} else {status = 5;}}
			}

		}
		//if ((final_sigma > 0.1)||(SouthTowerPeak[i]<0.09)||(SouthTowerPeak[i]>0.18)) {SouthTowerPeak[i]=peakmass;}
		cout<<"Peak:"<<SouthTowerPeak[i]<<endl;
		SouthGainCorr[i] = oldSouthGainCorr[i]*(0.135 / SouthTowerPeak[i]);
		SouthGainCorrRatio[i]=0.135 / SouthTowerPeak[i];
		cout<<"South tower gain corr:"<<SouthGainCorr[i]<<endl;
		h2_gaincorrectionS->Fill(-col,-row,SouthGainCorr[i]);
		statusS->Fill(-col,-row,status);
		if (status >= 2){h2_masspeakS->Fill(-col,-row,SouthTowerPeak[i]);}
		if (status==5) { h1_Southoldgaincorr->Fill(i,oldSouthGainCorr[i]); h1_Southnewgaincorr->Fill(i,SouthGainCorr[i]);h1_Southgaincorrdist->Fill(SouthGainCorr[i]);h1_Southgaincorrdist_used->Fill(oldSouthGainCorr[i]);}
		if (status==3) { h1_Southgaincorrdist_status3->Fill(SouthGainCorr[i]);}
		h1_Allgaincorrdist->Fill(SouthGainCorr[i]);
		h2_Southgaincorrration->Fill(-col,-row,SouthGainCorrRatio[i]);
		outputGainCorrOffline << ns <<"    " <<i<<"     "<<SouthGainCorr[i]<<endl;
		outputGainCorrRatio << ns <<"    " <<i<<"     "<<SouthGainCorr[i]<<"      "<<status<<endl;
		if (col==22)
		{

			if (row ==1) {c6[1]->SaveAs(Form("SouthTowerResult_iteration%s.pdf[",iteration));}
			if (row ==34) {c6[34]->SaveAs(Form("SouthTowerResult_iteration%s.pdf]",iteration));}
			if (!((row ==34)||(row ==1))) {c6[row]->SaveAs(Form("SouthTowerResult_iteration%s.pdf",iteration));}
			c6[row+1]=new TCanvas("c6","c6",6000,6000);
			c6[row+1].Divide(5,5);

		}

	}	
	c2->SaveAs("mass_by_tower_Scluster.jpg");
	statusS->SaveAs(Form("statusS_iteration%s.root",iteration));
	h2_gaincorrectionS->SaveAs("gaincorrectionS.root");
	h2_masspeakS->SaveAs("masspeakS.root");
	h2_Southgaincorrration->SaveAs("SouthGainCorrectionRatio.root");
	chi2ndfS->SaveAs("Chi2NDFS.root");

	for (ns=0; ns <2; ns=ns+1)
	{

		for (dep=0;dep<24;dep=dep+1)
		{

			for (ch=0;ch<32;ch=ch+1)
			{
				if (depchtoid[ns][dep][ch] == -1) continue;
				if (ns == 0)
				{
					outputGainCorr <<" 0  "<< ns <<"    "<< dep << "    "<<ch<<"     "<<NorthGainCorr[depchtoid[ns][dep][ch]]<<endl;
				}
				else
				{
					outputGainCorr <<" 0  "<< ns <<"    "<< dep << "    "<<ch<<"     "<<SouthGainCorr[depchtoid[ns][dep][ch]]<<endl;
				}

			}

		}

	}

	c3.cd();
	h1_Northoldgaincorr->SetMarkerStyle(20);
	h1_Northoldgaincorr->SetMarkerColor(2);
	h1_Northoldgaincorr->Draw("P");
	h1_Northnewgaincorr->SetMarkerStyle(20);
	h1_Northnewgaincorr->SetMarkerColor(3);
	h1_Northnewgaincorr->Draw("P same");
	c3->SaveAs("Northtowergaincorr.root");
	c4.cd();
	h1_Southoldgaincorr->SetMarkerStyle(20);
	h1_Southoldgaincorr->SetMarkerColor(2);
	h1_Southoldgaincorr->Draw("P");
	h1_Southnewgaincorr->SetMarkerStyle(20);
	h1_Southnewgaincorr->SetMarkerColor(3);
	h1_Southnewgaincorr->Draw("P same");
	c4->SaveAs("Southtowergaincorr.root");

	h1_Northgaincorrdist->SaveAs("Northgaincorrdist.root");
	h1_Southgaincorrdist->SaveAs("Southgaincorrdist.root");
	h1_Northgaincorrdist_used->SaveAs("Northgaincorrdist_used.root");
	h1_Southgaincorrdist_used->SaveAs("Southgaincorrdist_used.root");
	h1_Northgaincorrdist_status3->SaveAs("Northgaincorrdist_status3.root");
	h1_Southgaincorrdist_status3->SaveAs("Southgaincorrdist_status3.root");
	h1_Allgaincorrdist->SaveAs("Allgaincorrdist.root");

	fclose(tower_input);
	depchid_input.close();
	inputGainCorr.close();
	outputGainCorr.close();
	outputGainCorrOffline.close();


}
