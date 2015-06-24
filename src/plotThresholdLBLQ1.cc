#include "CommonHead.h"
#include "CommonFunc.h"

using namespace std;


int MarkerStyleWheel(int idx){
  switch(idx){
  case 0: return 20;

  case 1: return 21;
  case 2: return 25;
    
  case 3: return 22;
  case 4: return 26;
    
  case 5: return 23;
  case 6: return 32;
    
  default: return idx+1;
  }
}

void exportMaskMap(TH2 *map, TString outputName){
  ofstream fout(outputName);
  int nbinx=map->GetNbinsX();
  int nbiny=map->GetNbinsY();

  for(int ibiny=1;ibiny<=nbiny;ibiny++){
    for(int ibinx=1;ibinx<=nbinx;ibinx++){
      fout<<(map->GetBinContent(ibinx,ibiny)>0.5?1:0)<<" ";
    }
    fout<<endl;
  }
  fout.close();
}

int main(int argc, char** argv){
  if(argc<3){
    cout<<"Usage: "<<argv[0]<<" <jobname> <input base dir>"<<endl;
    return 0;
  }
  
  TString jobname=argv[1];
  TString inputBaseDir=argv[2];
  
  TString outputDir="fig/threshold/onStave_LBLQ1/"+jobname;  

  const int nstep=16;
  const int nchip=2;
    
  TString thresholdHistNames[nchip]={"threshdist_Mod_101", "threshdist_Mod_102"};
  TString FE[nchip]={"FE1", "FE2"};

  double AltFines[nstep]={0};

  for(int istep=0;istep<nstep;istep++){
    // if(istep<35) AltFines[istep]=65+istep;
    // else 
    AltFines[istep]=65+istep*5;
  }

  int startIdx=1;

  TString inputThresholdFileNames[nchip][nstep];


  for(int istep=0;istep<nstep;istep++){
    for(int ichip=0;ichip<nchip;ichip++){
      inputThresholdFileNames[ichip][istep]=inputBaseDir+Form("/THRESHOLD_%s_AltFine%.0f/analysis.root", FE[ichip].Data(), AltFines[istep]);
    }
  }

  double thresholds[nchip][nstep]={{0}}, thresRMS[nchip][nstep]={{0}}, zeros[nstep]={0};

  system("mkdir -vp "+outputDir);

  TH1 *hthresholds[nchip][nstep]={{NULL}};

  for(int ichip=0;ichip<nchip;ichip++){
    for(int istep=0; istep<nstep; istep++){
      TFile *ft=TFile::Open(inputThresholdFileNames[ichip][istep], "r");
      if(!ft) continue;
      TH1D *ht=(TH1D*)ft->Get(thresholdHistNames[ichip]);
      hthresholds[ichip][istep]=(TH1*)ht->Clone();
      // ht->SetBinContent(1,0);	// Remove potential pile-up at 0
      TF1 *func=hthresholds[ichip][istep]->GetFunction("gauss");
      thresholds[ichip][istep]=func->GetParameter(1);
      thresRMS[ichip][istep]=func->GetParameter(2);
    }
  }

  CommonFunc::SetAtlasStyle();
  // Plot threshold vs. Altfine
  vector<TString> pavetext;
  pavetext.push_back("LBLQ1 on-stavelet test");
  pavetext.push_back("40 Mbps");

  TCanvas *c=new TCanvas("ThresholdVsAltFine","",800,600);
  TH1D *hframe=new TH1D("hframe","",1000, 60, 140);
  hframe->SetLineWidth(0);
  hframe->GetXaxis()->SetTitle("AltFine");
  hframe->GetYaxis()->SetTitle("Threshold / e");
  hframe->SetMinimum(800);
  hframe->SetMaximum(3800);
  hframe->Draw();
  
  c->SetGrid();

  double x1, y1, x2, y2;
  CommonFunc::GetX1Y1X2Y2(c, x1, y1, x2, y2);

  TLegend *leg=CommonFunc::FastLegend(x1+0.03,y2-0.3,x1+0.4,y2-0.02,0.05);

  TGraphErrors *gr[nchip];

  for(int ichip=0;ichip<nchip;ichip++){
    int nstep_ths=nstep;
    double AltFine_ths[nstep]={0};
    double mean_ths[nstep]={0}, rms_ths[nstep]={0};
	
    int idx=0;
    for(int istep=0;istep<nstep; istep++){
      if(thresholds[ichip][istep]!=0){
	AltFine_ths[idx]=AltFines[istep];
	mean_ths[idx]=thresholds[ichip][istep];
	rms_ths[idx]=thresRMS[ichip][istep];
	idx++;
      }
      else{
	nstep_ths--;
      }
    }
      
    gr[ichip]=new TGraphErrors(nstep_ths, AltFine_ths, mean_ths, zeros, rms_ths);
    gr[ichip]->SetMarkerColor(CommonFunc::ColorWheel(ichip+1));
    gr[ichip]->SetLineColor(CommonFunc::ColorWheel(ichip+1));
    gr[ichip]->SetMarkerStyle(MarkerStyleWheel(ichip+1));
    gr[ichip]->Draw("PL");
	
    leg->AddEntry(gr[ichip],FE[ichip],"lep");

    ofstream fout(outputDir+"/"+FE[ichip]+".txt");
    TSpline5 *s=new TSpline5("grs",gr[ichip]);
    for(int altfine=65;altfine<140;altfine++){
      fout<<altfine<<"\t"<<s->Eval(altfine)<<endl;
    }
    s->SetLineWidth(2);
    s->SetLineColor(CommonFunc::ColorWheel(ichip+1));
    s->Draw("same");
    fout.close();
  }

  leg->Draw();
      
  TPaveText* text=CommonFunc::CreatePaveText(x2-0.4, y1+0.02, x2-0.03, y1+0.2, pavetext, 0.04);
  text->Draw("same");

  c->RedrawAxis();
  TString output1DName=outputDir+"/ThresholdVsAltFine";

  CommonFunc::PrintCanvas(c,output1DName);

  {
    for(int ichip=0;ichip<nchip;ichip++){
      vector<TString> pavetext;
      pavetext.push_back("LBLQ1 on-stavelet test");
      pavetext.push_back(FE[ichip]+" 40 Mbps");
      
      ofstream fout(outputDir+"/summary_thresholds_"+FE[ichip]+".txt");
      TCanvas *c=new TCanvas("ThresholdHists_"+FE[ichip],"",800,600);
      TH1D *hframe=new TH1D("hframe","",1000, 800, 3800);
      hframe->SetLineWidth(0);
      hframe->GetXaxis()->SetTitle("Threshold / e");
      hframe->GetYaxis()->SetTitle("Pixels");
      hframe->SetMinimum(0);
      hframe->SetMaximum(5000);
      hframe->Draw();
  
      c->SetGrid();

      double x1, y1, x2, y2;
      CommonFunc::GetX1Y1X2Y2(c, x1, y1, x2, y2);

      TPaveText* text=CommonFunc::CreatePaveText(x1+0.03, y2-0.5, x1+0.4, y2-0.3, pavetext, 0.04);
      text->Draw("same");

      TLegend *leg=CommonFunc::FastLegend(x1+0.03,y2-0.25,x1+0.75,y2-0.02,0.035);
      leg->SetNColumns(4);
      TGraphErrors *gr[nchip];
      int index=1;
      for(int istep=0;istep<nstep;istep++){
  	if(int(AltFines[istep])%5!=0) continue;
  	hthresholds[ichip][istep]->SetLineWidth(3);
  	hthresholds[ichip][istep]->SetMarkerColor(CommonFunc::ColorWheel(index));
  	hthresholds[ichip][istep]->SetLineColor(CommonFunc::ColorWheel(index));
  	hthresholds[ichip][istep]->SetMarkerStyle(MarkerStyleWheel(index));
  	//hthresholds[ichip][istep]->Draw("hist,same");
  	hthresholds[ichip][istep]->Draw("E,same");
  	TF1 *func=hthresholds[ichip][istep]->GetFunction("gauss");
  	cout<<func->GetParameter(0)<<" "<<func->GetParameter(1)<<" "<<func->GetParameter(2)<<endl;
  	fout<<Form("%.0f\t%.2f\t%.2f", AltFines[istep], func->GetParameter(1), func->GetParameter(2))<<endl;;
  	func->SetLineColor(CommonFunc::ColorWheel(index));
  	leg->AddEntry(hthresholds[ichip][istep],Form("AltFine=%.0f", AltFines[istep]),"lep");
  	index++;
      }
      fout.close();
      leg->Draw();
      TString output1DName=outputDir+"/ThresholdHists_"+FE[ichip];

      CommonFunc::PrintCanvas(c,output1DName);
    }
  }
}
