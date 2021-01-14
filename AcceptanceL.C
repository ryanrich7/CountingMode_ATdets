
//Author -- Ryan Richards
////Looks at Qsq,Lab Phi, Cos(Phi),Sin(Phi) and theta,phi_tg for Main, A_T detectors

//adc_up - main detector
//adc_in - A_T in, adc_out - A_T out
//Ebeam in GeV - Run dependent...
void AcceptanceL(int run, double th0, double adc_up, double adc_in, double adc_out, double Ebeam){ 

   gStyle->SetOptStat(1002211);


   TChain *T = new TChain("T");
   T->Add(Form("/lustre19/expphy/volatile/halla/parity/crex_optics_rootfiles/prexLHRS_%d_-1.root",run));

   double d2r = TMath::Pi()/180; 
   double r2d = 1/d2r;

   //th0 is for data. I will reconstruct lab phi for this analysis
   double cth0 = TMath::Cos(th0*d2r); 
   double sth0  = TMath::Sin(th0*d2r);  
   double tanth0 = sth0/cth0;
   double secth0 = 1/cth0;



   //[0] - main det, [1]- A_T in, [2] - A_T out
   TH1F *Phi[3],*Cos[3], *Sin[3];
   TH1F *ThTg[3],*PhTg[3];
   TH1F *Qsq[3];

   TCanvas *c_q, *c_ph, *c_cos, *c_sin, *c_thtg, *c_phtg;


   int color[3] = {2,3,4};
   char det[3][20] = {"Main","ATin","ATout"};

   for(int i = 0; i < 3; i++){ 
 
      Qsq[i] = new TH1F(Form("Qsq[%i]",i),Form("LHRS Q^{2} for %s Detector, Run %d",det[i],run),150,0,0.07);
      //Relative to horizontal 
      Phi[i] = new TH1F(Form("Phi[%i]",i),Form("LHRS #Phi_{H} for %s Detector, Run %d",det[i],run),150,-50,50);
     //Relative to vertical
     Cos[i] = new TH1F(Form("Cos[%i]",i),Form("LHRS Cos(#Phi_{V}) for %s Detector, Run %d",det[i],run),150,-0.8,0.8);
     Sin[i] = new TH1F(Form("Sin[%i]",i),Form("LHRS Sin(#Ph_{V}) for %s Detector, Run %d",det[i],run),150,0.6,1.05);
     //Thtg and Phtg
     ThTg[i] = new TH1F(Form("ThTg[%i]",i),Form("LHRS #theta_{tg} for %s Detector, Run %d",det[i],run),150,-0.06,0.06);
     PhTg[i] = new TH1F(Form("PhTg[%i]",i),Form("LHRS #phi_{tg} for %s Detector, Run %d",det[i],run),150,-0.03,0.03);

     Qsq[i]->SetLineColor(color[i]); Phi[i]->SetLineColor(color[i]); Cos[i]->SetLineColor(color[i]);
     Sin[i]->SetLineColor(color[i]); ThTg[i]->SetLineColor(color[i]); PhTg[i]->SetLineColor(color[i]);

    }

    //Relevant data from root tree
    double thisu1, thisv1, thisu2, thisv2;
    double thisPdat,thisdPdat, thisThdat, thisPhdat,thisEvt;
    double thisCosAngdat, thisuADC, thisQsqdat, thisCosPhi, thisSinPhi, thisPhi;
    double thisupADC,thisinADC,thisoutADC;


    //A_T in
//    double thisinCosAng,thisinADC,thisinQsq,thisinCosPhi, thisinSinPhi, thisinPhi;
    //A_T out 
//    double thisoutCosAng,thisoutADC,thisoutQsq,thisoutCosPhi,thisoutSinPhi,thisoutPhi;

   
 
    vector <double> Qsqmain, Qsqatin, Qsqatout;
    vector <double> Phimain, Phiatin, Phiatout;
    vector <double> Cosphimain, Cosphiatin, Cosphiatout;
    vector <double> Sinphimain, Sinphiatin, Sinphiatout;
    vector <double> thtgmain, thtgatin, thtgatout;
    vector <double> phtgmain, phtgatin, phtgatout;


    T->SetBranchAddress("L.vdc.u1.nclust",&thisu1); T->SetBranchAddress("L.vdc.v1.nclust",&thisv1);
    T->SetBranchAddress("L.vdc.u2.nclust",&thisu2); T->SetBranchAddress("L.vdc.v2.nclust",&thisv2);
    T->SetBranchAddress("L.gold.th",&thisThdat); T->SetBranchAddress("L.gold.ph",&thisPhdat);
    T->SetBranchAddress("L.gold.p",&thisPdat); T->SetBranchAddress("P.upQadcL",&thisupADC);
    T->SetBranchAddress("L.gold.dp",&thisdPdat); T->SetBranchAddress("P.evtypebits",&thisEvt);
    T->SetBranchAddress("P.atlQadcL",&thisinADC); T->SetBranchAddress("P.atrQadcL",&thisoutADC);

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
         Qsqmain.push_back(thisQsqdat);
         Cosphimain.push_back(thisCosPhi);
         Sinphimain.push_back(thisSinPhi);
         Phimain.push_back(thisPhi);
         thtgmain.push_back(thisThdat);
         phtgmain.push_back(thisPhdat);
      } if(thisinADC > adc_in){
         Qsqatin.push_back(thisQsqdat);
         Cosphiatin.push_back(thisCosPhi);
         Sinphiatin.push_back(thisSinPhi);
         Phiatin.push_back(thisPhi);
         thtgatin.push_back(thisThdat);
         phtgatin.push_back(thisPhdat);
      } if(thisoutADC > adc_out){     
         Qsqatout.push_back(thisQsqdat);
         Cosphiatout.push_back(thisCosPhi);
         Sinphiatout.push_back(thisSinPhi);
         Phiatout.push_back(thisPhi);
         thtgatout.push_back(thisThdat);
         phtgatout.push_back(thisPhdat);
      }



     }


  }


  for(int k = 0;  k < Qsqmain.size(); k++){     
    Qsq[0]->Fill(Qsqmain[k]);Cos[0]->Fill(Cosphimain[k]); Sin[0]->Fill(Sinphimain[k]); ThTg[0]->Fill(thtgmain[k]); PhTg[0]->Fill(phtgmain[k]);
    Phi[0]->Fill(Phimain[k]);   
  }
  for(int l = 0; l < Qsqatin.size(); l++){
    Qsq[1]->Fill(Qsqatin[l]); Cos[1]->Fill(Cosphiatin[l]); Sin[1]->Fill(Sinphiatin[l]); ThTg[1]->Fill(thtgatin[l]); PhTg[1]->Fill(phtgatin[l]);
    Phi[1]->Fill(Phiatin[l]); 
  }
  for(int m = 0; m < Qsqatout.size(); m++){
    Qsq[2]->Fill(Qsqatout[m]); Cos[2]->Fill(Cosphiatout[m]); Sin[2]->Fill(Sinphiatout[m]); ThTg[2]->Fill(thtgatout[m]); PhTg[2]->Fill(phtgatout[m]);
    Phi[2]->Fill(Phiatout[m]); 
   }




   c_q = new TCanvas();
   c_q->Divide(2,2);
   for(int i = 0; i < 3; i++){
   c_q->cd(i+1);
   Qsq[i]->Draw("HIST");
   }
   c_q->SaveAs(Form("qsqLHRS_Run%d.pdf",run));


   c_cos = new TCanvas();
   c_cos->Divide(2,2);
   for(int i = 0; i < 3; i++){
   c_cos->cd(i+1);
   Cos[i]->Draw("HIST");
   }
   c_cos->SaveAs(Form("cosLHRS_Run%d.pdf",run));
  
   c_sin = new TCanvas();
   c_sin->Divide(2,2);
   for(int i = 0; i < 3; i++){ 
   c_sin->cd(i+1);
   Sin[i]->Draw("HIST");
   }
   c_sin->SaveAs(Form("sinLHRS_Run%d.pdf",run));

   c_ph = new TCanvas();
   c_ph->Divide(2,2); 
   for(int i = 0; i < 3; i++){
   c_ph->cd(i+1);
   Phi[i]->Draw("HIST");
   }
   c_ph->SaveAs(Form("phiLHRS_Run%d.pdf",run));

   c_thtg = new TCanvas();
   c_thtg->Divide(2,2);
   for(int i = 0; i < 3; i++){
   c_thtg->cd(i+1);
   ThTg[i]->Draw("HIST");
   }
   c_thtg->SaveAs(Form("thtgLHRS_Run%d.pdf",run));

   c_phtg = new TCanvas();
   c_phtg->Divide(2,2);
   for(int i = 0; i < 3; i++){
   c_phtg->cd(i+1);
   PhTg[i]->Draw("HIST");
   }
   c_phtg->SaveAs(Form("phtgLHRS_Run%d.pdf",run));




    


    



}
