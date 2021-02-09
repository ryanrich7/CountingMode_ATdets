
//Author -- Ryan Richards
////Looks at Qsq,Lab Phi, Cos(Phi),Sin(Phi) and theta,phi_tg for Main, A_T detectors

//adc_up - upstream detector
//adc_lo - downstream detector
//Ebeam in GeV - Run dependent...
void ATRunL(int run, double th0, double adc_up, double adc_lo, double Ebeam, TString targ){ 

   gStyle->SetOptStat(1002211);


   TChain *T = new TChain("T");
  // T->Add(Form("/lustre19/expphy/volatile/halla/parity/crex_optics_rootfiles/prexLHRS_%d_-1*.root",run));
   T->Add(Form("../prex2files/prexLHRS_%d_-1*.root",run));

//   T->Add(Form("/lustre19/expphy/volatile/halla/parity/crex_optics_rootfiles/prexLHRS_%d_-1_1.root",run));
//   T->Add(Form("/lustre19/expphy/volatile/halla/parity/crex_optics_rootfiles/prexLHRS_%d_-1_2.root",run));


   double d2r = TMath::Pi()/180; 
   double r2d = 1/d2r;

   //th0 is for data. I will reconstruct lab phi for this analysis
   double cth0 = TMath::Cos(th0*d2r); 
   double sth0  = TMath::Sin(th0*d2r);  
   double tanth0 = sth0/cth0;
   double secth0 = 1/cth0;



   //[0] - up main det, [1]- down main det
   TH1F *Phi[2],*Cos[2], *Sin[2];
   TH1F *Qsq[2];

   TCanvas *c_q, *c_ph, *c_cos, *c_sin, *c_thtg, *c_phtg;


   int color[2] = {2,3};
   char det[2][20] = {"Upstream","Downstream"};

   for(int i = 0; i < 2; i++){ 
 
      Qsq[i] = new TH1F(Form("Qsq[%i]",i),Form("LHRS Q^{2} for %s Main Detector, Run %d, %s",det[i],run,targ.Data()),150,0,0.015);//0.07 for CREX
      //Relative to horizontal 
      Phi[i] = new TH1F(Form("Phi[%i]",i),Form("LHRS #Phi_{H} for %s Main Detector, Run %d, %s",det[i],run, targ.Data()),150,-50,50);
     //Relative to vertical
     Cos[i] = new TH1F(Form("Cos[%i]",i),Form("LHRS Cos(#Phi_{V}) for %s Main Detector, Run %d, %s",det[i],run, targ.Data()),150,-0.8,0.8);
     Sin[i] = new TH1F(Form("Sin[%i]",i),Form("LHRS Sin(#Phi_{V}) for %s Main Detector, Run %d, %s",det[i],run, targ.Data()),150,0.6,1.05);

     Qsq[i]->SetLineColor(color[i]); Phi[i]->SetLineColor(color[i]); Cos[i]->SetLineColor(color[i]);
     Sin[i]->SetLineColor(color[i]); 


     Qsq[i]->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
     Sin[i]->GetXaxis()->SetTitle("Sin(#Phi_{V})");
     Cos[i]->GetXaxis()->SetTitle("Cos(#Phi_{V})");
     Phi[i]->GetXaxis()->SetTitle("#Phi_{H} (deg)");



    }

    //Relevant data from root tree
    double thisu1, thisv1, thisu2, thisv2;
    double thisPdat,thisdPdat, thisThdat, thisPhdat,thisEvt;
    double thisCosAngdat, thisuADC, thisQsqdat, thisCosPhi, thisSinPhi, thisPhi;
    double thisupADC, thisloADC;


    //A_T in
//    double thisinCosAng,thisinADC,thisinQsq,thisinCosPhi, thisinSinPhi, thisinPhi;
    //A_T out 
//    double thisoutCosAng,thisoutADC,thisoutQsq,thisoutCosPhi,thisoutSinPhi,thisoutPhi;

   
 
    vector <double> Qsqup, Qsqdown;
    vector <double> Phiup,Phidown;
    vector <double> Cosphiup,Cosphidown;
    vector <double> Sinphiup, Sinphidown;


    T->SetBranchAddress("L.vdc.u1.nclust",&thisu1); T->SetBranchAddress("L.vdc.v1.nclust",&thisv1);
    T->SetBranchAddress("L.vdc.u2.nclust",&thisu2); T->SetBranchAddress("L.vdc.v2.nclust",&thisv2);
    T->SetBranchAddress("L.gold.th",&thisThdat); T->SetBranchAddress("L.gold.ph",&thisPhdat);
    T->SetBranchAddress("L.gold.p",&thisPdat); T->SetBranchAddress("P.upQadcL",&thisupADC);
    T->SetBranchAddress("L.gold.dp",&thisdPdat); T->SetBranchAddress("P.evtypebits",&thisEvt);
    T->SetBranchAddress("P.loQadcL",&thisloADC);


   long n = T->GetEntries();
   
   //Looping over tree
   for(long j = 0; j < n; j++){
      T->GetEntry(j);
     
      int thisevent = (int) thisEvt;
      if(thisu1 == 1 && thisv1 == 1 && thisu2 == 1 && thisv2 == 1 && abs(thisThdat)<0.08 && abs(thisPhdat)<0.05 && thisPdat > Ebeam*0.96 && thisPdat < Ebeam*1.002 && ((thisevent&2)==2) ) {
 
      //Will deal with ADC cuts separately  
      //Qsq
      thisCosAngdat = (cth0 - thisPhdat*sth0)/(TMath::Sqrt(1+thisThdat*thisThdat+thisPhdat*thisPhdat));
      thisQsqdat = 2*thisPdat*Ebeam*(1-thisCosAngdat);  
      //CosPhi relative to vertical
      thisCosPhi = thisThdat/(TMath::Sqrt(thisThdat*thisThdat+sth0*sth0+2*tanth0*thisPhdat + thisPhdat*thisPhdat*secth0*secth0));     
      thisSinPhi = TMath::Sqrt(1-thisCosPhi*thisCosPhi);     
      //Phi relative in horizontal
      thisPhi = 90 - r2d*TMath::ACos(thisCosPhi);

      if(thisupADC > adc_up){ 
         Qsqup.push_back(thisQsqdat);
         Cosphiup.push_back(thisCosPhi);
         Sinphiup.push_back(thisSinPhi);
         Phiup.push_back(thisPhi);
      } if(thisloADC > adc_lo){
         Qsqdown.push_back(thisQsqdat);
         Cosphidown.push_back(thisCosPhi);
         Sinphidown.push_back(thisSinPhi);
         Phidown.push_back(thisPhi);
      }



     }


  }


  for(int k = 0;  k < Qsqup.size(); k++){     
    Qsq[0]->Fill(Qsqup[k]);Cos[0]->Fill(Cosphiup[k]); Sin[0]->Fill(Sinphiup[k]); Phi[0]->Fill(Phiup[k]);   
  }
  for(int l = 0; l < Qsqdown.size(); l++){
    Qsq[1]->Fill(Qsqdown[l]); Cos[1]->Fill(Cosphidown[l]); Sin[1]->Fill(Sinphidown[l]); Phi[1]->Fill(Phidown[l]); 
  }


   c_q = new TCanvas();
   c_q->Divide(2,1);
   for(int i = 0; i < 2; i++){
   c_q->cd(i+1);
   Qsq[i]->Draw("HIST");
   }
//   c_q->SaveAs(Form("qsqLHRS_Run%d.pdf",run));


   c_cos = new TCanvas();
   c_cos->Divide(2,1);
   for(int i = 0; i < 2; i++){
   c_cos->cd(i+1);
   Cos[i]->Draw("HIST");
   }
  // c_cos->SaveAs(Form("cosLHRS_Run%d.pdf",run));
  
   c_sin = new TCanvas();
   c_sin->Divide(2,1);
   for(int i = 0; i < 2; i++){ 
   c_sin->cd(i+1);
   Sin[i]->Draw("HIST");
   }
  // c_sin->SaveAs(Form("sinLHRS_Run%d.pdf",run));

   c_ph = new TCanvas();
   c_ph->Divide(2,1); 
   for(int i = 0; i < 2; i++){
   c_ph->cd(i+1);
   Phi[i]->Draw("HIST");
  }
 //  c_ph->SaveAs(Form("phiLHRS_Run%d.pdf",run));


   TCanvas *ctest = new TCanvas();
   ctest->Divide(2,2);
   ctest->cd(1);
   Qsq[0]->Draw("HIST");
   ctest->cd(2);
   Cos[0]->Draw("HIST");
   ctest->cd(3);
   Sin[0]->Draw("HIST");
   ctest->cd(4);
   Phi[0]->Draw("HIST");
   ctest->SaveAs(Form("KineLHRSupdet_Run%d_%s.pdf",run,targ.Data()));

   TCanvas *ctest1 = new TCanvas();
   ctest1->Divide(2,2);
   ctest1->cd(1);
   Qsq[1]->Draw("HIST");
   ctest1->cd(2);
   Cos[1]->Draw("HIST");
   ctest1->cd(3);
   Sin[1]->Draw("HIST");
   ctest1->cd(4);
   Phi[1]->Draw("HIST");
   ctest1->SaveAs(Form("KineLHRSdowndet_Run%d_%s.pdf",run,targ.Data())); 


    



}
