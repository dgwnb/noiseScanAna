#include "CommonHead.h"
#include "CommonFunc.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"

using namespace std;
using namespace CommonFunc;

const double epsilon=1e-6;
const int MAXSTEP=200;
const int MAXSAMPLE=10;

const int nGangPixRow=7;
  
const int nPixX=80, nPixY=336;

const double nPixel=26880;

const double noiseThreshold=1e-12;

const double burstThreshold=0.1;
  
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

map <TString, TGraph *>grMap;

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

void ANDMaskMap(double input[nPixX][nPixY], double standard[nPixX][nPixY]){

  for(int ibiny=0;ibiny<nPixY;ibiny++){
    for(int ibinx=0;ibinx<nPixX;ibinx++){
      double content=(input[ibinx][ibiny]>0.5)&&(standard[ibinx][ibiny]>0.5);
      input[ibinx][ibiny]=content;
    }
  }
}

void readMaskMap(TString inputFileName, double PixMask[nPixX][nPixY]){
  ifstream fin(inputFileName);
  assert(fin);
  cout<<"REGTEST: Using input mask file "<<inputFileName<<endl;
  for(int ibiny=0;ibiny<nPixY;ibiny++){
    for(int ibinx=0;ibinx<nPixX;ibinx++){
      fin>>PixMask[ibinx][ibiny];
    }
  }
  fin.close();
}

double BinarySearch(TF1 *s, double target, double xmin, double xmax, double tolerance=1e-3){
  if(xmin<s->GetXmin()){
    cout<<"Warning: input range is too wide. Use the range of the input spline: "
	<<s->GetXmin()<<" "<<xmax<<endl;
    xmin=s->GetXmin();
  }
  if(xmax>s->GetXmax()){
    cout<<"Warning: input range is too wide. Use the range of the input spline: "
	<<xmin<<" "<<s->GetXmax()<<endl;
    xmax=s->GetXmax();
  }
  if(xmin>xmax){
    cout<<"wrong input."<<endl;
    cout<<"xmin: "<<xmin<<", xmax: "<<xmax<<endl;
    abort();
  }
  // start binary search

  double begin = xmin;
  double end = xmax;
  double median=(end-begin)/2;
  double thisUL=s->Eval(begin);
  int counter=0;
  while(fabs(thisUL-target)>tolerance){
    median=(end+begin)/2;

    thisUL=s->Eval(median);
    if((thisUL-target)*(s->Eval(end)-target)<0) begin=median;
    else if((thisUL-target)*(s->Eval(begin)-target)<0) end=median;
    counter++;
  }
  return median;//   return median;
}

double BinarySearch(TSpline *s, double target, double xmin, double xmax, double tolerance=1e-3){
  if(xmin<s->GetXmin()){
    cout<<"Warning: input range is too wide. Use the range of the input spline: "
	<<s->GetXmin()<<" "<<xmax<<endl;
    xmin=s->GetXmin();
  }
  if(xmax>s->GetXmax()){
    cout<<"Warning: input range is too wide. Use the range of the input spline: "
	<<xmin<<" "<<s->GetXmax()<<endl;
    xmax=s->GetXmax();
  }
  if(xmin>xmax){
    cout<<"wrong input."<<endl;
    cout<<"xmin: "<<xmin<<", xmax: "<<xmax<<endl;
    abort();
  }
  // start binary search

  double begin = xmin;
  double end = xmax;
  double median=(end-begin)/2;
  double thisUL=s->Eval(begin);
  int counter=0;
  while(fabs(thisUL-target)>tolerance){
    median=(end+begin)/2;

    thisUL=s->Eval(median);
    if((thisUL-target)*(s->Eval(end)-target)<0) begin=median;
    else if((thisUL-target)*(s->Eval(begin)-target)<0) end=median;
    counter++;
  }
  return median;//   return median;
}

double translateAltFine(double altfine, TString FE){
  if(!grMap[FE]){
    // Create the TGraph as map
    cout<<"Create the TGraph as map"<<endl;
    ifstream fin("Input/threshold/onStave_LBLQ1/"+FE+".txt");
    assert(fin);
    double AltFine[MAXSTEP], Threshold[MAXSTEP];
    int idx=0;
    while(fin>>AltFine[idx]>>Threshold[idx]) idx++;
    grMap[FE]=new TGraph(idx, AltFine, Threshold);
    fin.close();
  }
  
  return grMap[FE]->Eval(altfine);
}

double translateThreshold(double threshold, TString FE){
  if(!grMap[FE]){
    // Create the TGraph as map
    cout<<"Create the TGraph as map"<<endl;
    ifstream fin("Input/threshold/onStave_LBLQ1/"+FE+".txt");
    double AltFine[MAXSTEP], Threshold[MAXSTEP];
    int idx=0;
    while(fin>>AltFine[idx]>>Threshold[idx]) idx++;
    grMap[FE]=new TGraph(idx, AltFine, Threshold);
    fin.close();
  }
  return BinarySearch(new TSpline5("spline", grMap[FE]), threshold, 65, 90, 1e-3);
}

bool checkNoiseBurst(TH2 *h, TString option="noGang"){
  int nbinx=h->GetNbinsX();
  int nbiny=h->GetNbinsY();
  int fired=0;
  for(int ibinx=1;ibinx<=nbinx;ibinx++){
    for(int ibiny=1;ibiny<=nbiny;ibiny++){
      if(option=="noGang"&&ibiny<=nGangPixRow) continue;
      if(h->GetBinContent(ibinx,ibiny)>0) fired++;
    }
  }

  double total=nbinx*nbiny;
  if(option=="noGang") total-=nGangPixRow*nPixX;
  return fired/total>burstThreshold;
}


double GausCDF_C(double *x, double *p){
  return ROOT::Math::gaussian_cdf_c(x[0], p[1], p[0])*p[2];
}

int main(int argc, char** argv){
  if(argc<3){
    cout<<"Usage: "<<argv[0]<<" <jobname> <input list>"<<endl;
    return 0;
  }

  const int nchip=2;
  
  TString jobname=argv[1];
  TString inputList=argv[2];
  
  TString option=argv[argc-1];

  TString outputDir="fig/noiseScan/onStave_LBLQ1/"+jobname;
  
  TString histNames[nchip]={"Mod_101_NoiseOccupancy","Mod_102_NoiseOccupancy"};
  TString nTriggerHistName[nchip]={"Mod_101_nEvents","Mod_102_nEvents"};

  TString FE[nchip]={"FE1", "FE2"};

  system("mkdir -vp "+outputDir);
  
  ifstream fin(inputList,ios::in);
  vector<TString> inputSet, inputArg;
  vector<TString> legtext,legopt, Opt;
  vector<int> style,color;
  vector<double> width;
    
  TString fileName, argument, legName, optName, styleName, colorName, widthName, optionName;
  while(fin>>fileName>>argument>>legName>>optName>>styleName>>colorName>>widthName>>optionName){
    legtext.push_back(legName);
    legopt.push_back(optName);
      
    inputSet.push_back(fileName);
    inputArg.push_back(argument);
    style.push_back(atoi(styleName));
    color.push_back(atoi(colorName));
    width.push_back(atof(widthName));

    Opt.push_back(optionName);
  }

  fin.close();

  int nsample=inputSet.size();
  vector<TString> inputlist=SplitString(inputArg[0],',');
  TString poiName=inputlist[0];
  double xmin=atof(inputlist[1]);
  double xmax=atof(inputlist[2]);

  for(int isam=0;isam<nsample;isam++){
    vector<TString> inputlist=SplitString(inputArg[isam],',');
    double start=atof(inputlist[1]);
    double end=atof(inputlist[2]);

    if(start<xmin) xmin=start;
    if(end>xmax) xmax=end;
  }

  double AltFines[MAXSAMPLE][MAXSTEP]={{0}};

  TCanvas *cmask[nchip]={NULL};
  TString maskOutputName[nchip];
  
  for(int ichip=0;ichip<nchip;ichip++){
    cmask[ichip]=new TCanvas("Mask"+FE[ichip],"",800,600);
    maskOutputName[ichip]=outputDir+"/Mask"+FE[ichip]+".gif";
    gSystem->Unlink(maskOutputName[ichip]);
  }

  int nstep[MAXSAMPLE]={0};
  double rangeMin[MAXSAMPLE][nchip]={{-1}}, rangeMax[MAXSAMPLE][nchip]={{-1}};
  double entries[MAXSAMPLE][nchip][MAXSTEP], uncert[MAXSAMPLE][nchip][MAXSTEP];
  double entries_noGang[MAXSAMPLE][nchip][MAXSTEP], uncert_noGang[MAXSAMPLE][nchip][MAXSTEP];
  double entries_noHot[MAXSAMPLE][nchip][MAXSTEP], uncert_noHot[MAXSAMPLE][nchip][MAXSTEP]; // This is the baseline figure of merit
  double entries_noEither[MAXSAMPLE][nchip][MAXSTEP], uncert_noEither[MAXSAMPLE][nchip][MAXSTEP];

  double Xmax[nchip]={-1};

  TString inputMaskFile[nchip];

  for(int ichip=0;ichip<nchip;ichip++) inputMaskFile[ichip]="Input/hotPixMask/onStave_LBLQ1/Mask"+FE[ichip]+".txt";
  TString inputStdMaskFile="Input/hotPixMask/onStave_LBLQ1/Mask_standard_with_ganged_pixel.txt";
  double StdPixMask[nPixX][nPixY];
  readMaskMap(inputStdMaskFile, StdPixMask);
  
  for(int isam=0;isam<nsample;isam++){
    cout<<"Processing "<<inputSet[isam]<<endl;
    if(Opt[isam].Contains("stdmask")){
      cout<<"\tREGTEST: Using standard mask"<<endl;
      for(int ichip=0;ichip<nchip;ichip++){
	inputMaskFile[ichip]="Input/hotPixMask/onStave_LBLQ1/Mask_standard_with_ganged_pixel.txt";
      }
    }
    vector<TString> inputlist=SplitString(inputArg[isam],',');
    poiName=inputlist[0];
    double start=atof(inputlist[1]);
    double end=atof(inputlist[2]);
    double binw=atof(inputlist[3]);
    nstep[isam]=int((end-start+epsilon)/binw)+1;

    for(int istep=0;istep<nstep[isam];istep++){
      AltFines[isam][istep]=start+istep*binw;
    }

    TString inputBaseDir=inputSet[isam];
    TString inputFileNames[nchip][MAXSTEP];
    
    for(int istep=0;istep<nstep[isam];istep++){
      for(int ichip=0;ichip<nchip;ichip++){
	inputFileNames[ichip][istep]=inputBaseDir+Form("/NOISESCAN_%s_AltFine%.0f/histos.root", FE[ichip].Data(), AltFines[isam][istep]);
      }
    }
    
    TString legText[nchip];
    for(int ichip=0;ichip<nchip;ichip++) legText[ichip]=FE[ichip]+legtext[isam];

    double PixMask[nchip][nPixX][nPixY]={{{0}}};

    double firedpix[nchip][MAXSTEP]={{0}}, firedpix_noGang[nchip][MAXSTEP]={{0}}, firedpix_noHot[nchip][MAXSTEP]={{0}}, firedpix_noEither[nchip][MAXSTEP]={{0}};

    TCanvas *csanim[nchip]={NULL};
  
    TString outputAnimFileFE[nchip];

    // Here we load the mask from the nominal scan to make fair comparison
    
    for(int ichip=0;ichip<nchip;ichip++){
      rangeMin[isam][ichip]=-1;
      rangeMax[isam][ichip]=-1;

      TFile *f[MAXSTEP]={NULL};
      TH2 *hin[MAXSTEP]={NULL};
      TH1 *hnevt[MAXSTEP]={NULL};
      
      for(int istep=0;istep<nstep[isam];istep++){
	f[istep]=TFile::Open(inputFileNames[ichip][istep],"r");
	hin[istep]=(TH2I*)f[istep]->Get(histNames[ichip]);
	hnevt[istep]=(TH1*)f[istep]->Get(nTriggerHistName[ichip]);
      }

      if(Opt[isam].Contains("refmask")||Opt[isam].Contains("stdmask")){
	readMaskMap(inputMaskFile[ichip], PixMask[ichip]);
      }
      else{
	cout<<"REGTEST: Use its own mask file"<<endl;
      }
      // First do a quick analysis: Figure out threshold etc.

      int threshold=10;
      int start_step=-1;
      for(int istep=0;istep<nstep[isam];istep++){
	if(checkNoiseBurst(hin[istep])) continue;
	else{
	  cout<<"REGTEST: noise burst stopped at "<<AltFines[isam][istep]<<endl;
	  start_step=istep;
	  break;
	}
      }

      for(int ibinx=1;ibinx<=nPixX;ibinx++){
	for(int ibiny=1;ibiny<=nPixY;ibiny++){
	  int counter=0;
	  bool previousHit=false;
	  for(int istep=start_step;istep<nstep[isam];istep++){
	    if(!previousHit && counter<threshold) counter=0; // Reset the counter if there is no previous hits
	    double content=hin[istep]->GetBinContent(ibinx,ibiny);
	    if(content>0){
	      // Prevent to mask noisy pixels which only appear at the last step
	      if(istep<nstep[isam]-threshold) previousHit=true;
	      counter++;
	    }
	    else{
	      previousHit=false;
	    }
	    if(counter>=threshold){
	      if(!previousHit){
		cerr<<"ERROR: How come!?"<<endl; getchar();
		cout<<counter<<endl;
	      }
	      break;
	    }
	  }
	  if(!Opt[isam].Contains("refmask")&&!Opt[isam].Contains("stdmask")){
	    if(previousHit) PixMask[ichip][ibinx-1][ibiny-1]=0;
	    else PixMask[ichip][ibinx-1][ibiny-1]=1;
	  }
	}
      }
      ANDMaskMap(PixMask[ichip],StdPixMask); // We use the standard mask from now on
      
      outputAnimFileFE[ichip]=outputDir+"/OccupancyAnimation_"+legText[ichip]+".gif";

      csanim[ichip]=new TCanvas("animation"+legText[ichip],"",800,600);
      gSystem->Unlink(outputAnimFileFE[ichip]);
      csanim[ichip]->SetRightMargin(0.15);

      for(int istep=0; istep<nstep[isam]; istep++){

	TH2D *h=new TH2D("norm_occ_mod_"+legText[ichip]+Form("_AltFine=%.0f", AltFines[isam][istep]), legText[ichip]+Form(", AltFine=%.0f", AltFines[isam][istep]), 
			 hin[istep]->GetNbinsX(), hin[istep]->GetXaxis()->GetXmin(), hin[istep]->GetXaxis()->GetXmax(),
			 hin[istep]->GetNbinsY(), hin[istep]->GetYaxis()->GetXmin(), hin[istep]->GetYaxis()->GetXmax());
	TH2D *h_noGang=new TH2D("norm_occ_mod_noGang_"+legText[ichip]+Form("_AltFine=%.0f", AltFines[isam][istep]), "", 
				hin[istep]->GetNbinsX(), hin[istep]->GetXaxis()->GetXmin(), hin[istep]->GetXaxis()->GetXmax(),
				hin[istep]->GetNbinsY(), hin[istep]->GetYaxis()->GetXmin(), hin[istep]->GetYaxis()->GetXmax());
	TH2D *h_noHot=new TH2D("norm_occ_mod_noHot_"+legText[ichip]+Form("_AltFine=%.0f", AltFines[isam][istep]), "", 
			       hin[istep]->GetNbinsX(), hin[istep]->GetXaxis()->GetXmin(), hin[istep]->GetXaxis()->GetXmax(),
			       hin[istep]->GetNbinsY(), hin[istep]->GetYaxis()->GetXmin(), hin[istep]->GetYaxis()->GetXmax());
	TH2D *h_noEither=new TH2D("norm_occ_mod_noEither_"+legText[ichip]+Form("_AltFine=%.0f", AltFines[isam][istep]), "", 
				  hin[istep]->GetNbinsX(), hin[istep]->GetXaxis()->GetXmin(), hin[istep]->GetXaxis()->GetXmax(),
				  hin[istep]->GetNbinsY(), hin[istep]->GetYaxis()->GetXmin(), hin[istep]->GetYaxis()->GetXmax());

	// Analyze the raw data: Though calibGui did it for us, the uncertainty info is not saved!
	int nevt=hnevt[istep]->GetBinContent(1);
	int nbinx=hin[istep]->GetNbinsX();
	int nbiny=hin[istep]->GetNbinsY();

	double occ_noGang=0, occ_noHot=0, occ_noEither=0;
      
	double nPixel_noGang=nPixel-nGangPixRow*nPixX, nPixel_noHot=nPixel, nPixel_noEither=nPixel-nGangPixRow*nPixX;

	for(int ibinx=1;ibinx<=nbinx;ibinx++){
	  for(int ibiny=1;ibiny<=nbiny;ibiny++){
	    double content=hin[istep]->GetBinContent(ibinx,ibiny);
	    h->SetBinContent(ibinx, ibiny, content/nevt);
	    h_noGang->SetBinContent(ibinx, ibiny, content/nevt);
	    h_noHot->SetBinContent(ibinx, ibiny, content/nevt);
	    h_noEither->SetBinContent(ibinx, ibiny, content/nevt);

	    // cout<<content<<endl; getchar();
	    if(content>0) firedpix[ichip][istep]++;

	    bool isHot=false, isGang=false;
	  
	    if(ibiny<=nGangPixRow){
	      isGang=true;
	      h_noGang->SetBinContent(ibinx, ibiny, 0);
	      h_noEither->SetBinContent(ibinx, ibiny, 0);
	    }
	    else{
	      occ_noGang+=content;
	      if(content>0) firedpix_noGang[ichip][istep]++;
	    }
	  
	    if(PixMask[ichip][ibinx-1][ibiny-1]==0){
	      isHot=true;
	      h_noHot->SetBinContent(ibinx, ibiny, 0);
	      nPixel_noHot-=1;
	      nPixel_noEither-=1;
	    }
	    else{
	      occ_noHot+=content;
	      if(content>0) firedpix_noHot[ichip][istep]++;
	    }

	    if(!isHot&&!isGang){
	      occ_noEither+=content;
	      if(content>0) firedpix_noEither[ichip][istep]++;
	    }
	  }
	}

	firedpix[ichip][istep]/=nPixel/100.;
	firedpix_noGang[ichip][istep]/=nPixel_noGang/100.;
	firedpix_noHot[ichip][istep]/=nPixel_noHot/100.;
	firedpix_noEither[ichip][istep]/=nPixel_noEither/100.;

	gStyle->SetOptStat(0);

	entries[isam][ichip][istep]=h->Integral()/nPixel; // 26880: number of total pixels. 600s: Noise scan time
	uncert[isam][ichip][istep]=sqrt(hin[istep]->Integral())/nevt/nPixel;

	entries_noGang[isam][ichip][istep]=occ_noGang/nevt/nPixel_noGang;
	uncert_noGang[isam][ichip][istep]=sqrt(occ_noGang)/nevt/nPixel_noGang;

	entries_noHot[isam][ichip][istep]=occ_noHot/nevt/nPixel_noHot;
	uncert_noHot[isam][ichip][istep]=sqrt(occ_noHot)/nevt/nPixel_noHot;

	entries_noEither[isam][ichip][istep]=occ_noEither/nevt/nPixel_noEither;
	uncert_noEither[isam][ichip][istep]=sqrt(occ_noEither)/nevt/nPixel_noEither;

	// TCanvas *c=new TCanvas(histNames[ichip],"",800,600);

	// Starting point: number of fired pixels < 1%
	if(rangeMin[isam][ichip]<0&&firedpix_noGang[ichip][istep]<burstThreshold*100){
	  rangeMin[isam][ichip]=AltFines[isam][istep];
	  cout<<"Starting point for "<<FE[ichip]<<" decided: "
	      <<rangeMin[isam][ichip]<<", fraction="<<firedpix_noGang[ichip][istep]
	      <<endl;
	}

	// Ending point: occupancy<10*threshold
	if(rangeMax[isam][ichip]<0&&entries_noHot[isam][ichip][istep]<10*noiseThreshold){
	  rangeMax[isam][ichip]=AltFines[isam][istep];
	  cout<<"Ending point for "<<FE[ichip]<<" decided: "
	      <<rangeMax[isam][ichip]<<", occupancy="<<entries_noHot[isam][ichip][istep]
	      <<endl;
	}

	if(entries_noHot[isam][ichip][istep]>Xmax[ichip]&&isam==0) Xmax[ichip]=entries_noHot[isam][ichip][istep];
	TString histTitle="OCCUPANCY "+FE[ichip]+Form(" AltFine %.0f", AltFines[isam][istep]);
	csanim[ichip]->cd();
	h->SetTitle(histTitle);
	h->Draw("colz");
	csanim[ichip]->Modified();
	csanim[ichip]->Update();
	if(!option.Contains("noanim")){
	  if(istep<nstep[isam]-1) csanim[ichip]->Print(outputAnimFileFE[ichip]+"+50");
	  else csanim[ichip]->Print(outputAnimFileFE[ichip]+"++");
	}

	if(option.Contains("verbose")) cout<<inputFileNames[ichip][istep]<<" "<<histNames[ichip]<<" "<<entries[isam][ichip][istep]<<endl;
	delete h;

      }

      // Plot the mask
      cmask[ichip]->cd();
      TH2D *hframe=new TH2D("Mask_"+legText[ichip],"Mask_"+legText[ichip],nPixX, 0, nPixX, nPixY, 0, nPixY);
      hframe->GetXaxis()->SetTitle("Column");
      hframe->GetYaxis()->SetTitle("Row");

      for(int ibinx=1;ibinx<=nPixX;ibinx++){
	for(int ibiny=1;ibiny<=nPixY;ibiny++){
	  hframe->SetBinContent(ibinx, ibiny, PixMask[ichip][ibinx-1][ibiny-1]);
	}
      }
      hframe->DrawClone("colz");
      // cmask[ichip]->RedrawAxis();
      cmask[ichip]->Modified();
      cmask[ichip]->Update();
      if(!option.Contains("noanim")){
	if(isam<nsample-1) cmask[ichip]->Print(maskOutputName[ichip]+"+50");
	else cmask[ichip]->Print(maskOutputName[ichip]+"++");
      }
      if(Opt[isam].Contains("exportmask")){
	exportMaskMap(hframe, outputDir+"/Mask"+FE[ichip]+".txt");

	inputMaskFile[ichip]=outputDir+"/Mask"+FE[ichip]+".txt";
	cout<<"REGTEST: Export output mask file "<<inputMaskFile[ichip]<<endl;
      }
    }
    
  }

  // Plot noise vs. Altfine
  CommonFunc::SetAtlasStyle();
  double zeros[MAXSTEP]={0};
    
  TGraph *gr[nchip][MAXSAMPLE], *gr_noGang[nchip][MAXSAMPLE], *gr_noHot[nchip][MAXSAMPLE], *gr_noEither[nchip][MAXSAMPLE];
    
  TF1 *func[nchip][MAXSAMPLE]={{NULL}};
  
  for(int ichip=0;ichip<nchip;ichip++){
    vector<TString> pavetext;
    pavetext.push_back("LBLQ1 on-stavelet test");
    pavetext.push_back(FE[ichip]+" 40 Mbps");
    TCanvas *c=new TCanvas("NoiseVsAltFine"+FE[ichip],"",800,600);
    double xmin_th=translateAltFine(xmin, FE[ichip]);
    double xmax_th=translateAltFine(xmax, FE[ichip]);
    TH1D *hframe=new TH1D("hframe","",1000, xmin_th, xmax_th);
      
    hframe->SetLineWidth(0);
    hframe->GetXaxis()->SetTitle("Threshold / e");
    hframe->GetYaxis()->SetTitle("Averaged noise occupancy");
      
    hframe->SetMaximum(Xmax[ichip]*1.5);
    hframe->SetMinimum(5e-14);
    hframe->Draw();
      
    CommonFunc::DrawConstantLine(c, noiseThreshold, xmin_th, xmax_th,2,4,2);
    c->SetGrid();
      
    double x1, y1, x2, y2;
    CommonFunc::GetX1Y1X2Y2(c, x1, y1, x2, y2);
      
    TLegend *leg=CommonFunc::FastLegend(x2-0.5,y2-0.35,x2-0.03,y2-0.03,0.03);
    if(!option.Contains("nofit")) leg->SetFillStyle(1001);

    
    for(int isam=0;isam<nsample;isam++){
      double threshold[MAXSTEP];
      for(int istep=0;istep<nstep[isam];istep++){
	threshold[istep]=translateAltFine(AltFines[isam][istep], FE[ichip]);
      }
    
      gr_noHot[ichip][isam]=new TGraphErrors(nstep[isam], threshold, entries_noHot[isam][ichip], zeros, uncert_noHot[isam][ichip]);
      
      gr_noHot[ichip][isam]->SetMarkerColor(color[isam]);
      gr_noHot[ichip][isam]->SetLineColor(color[isam]);
      gr_noHot[ichip][isam]->SetLineStyle(style[isam]);
      gr_noHot[ichip][isam]->SetLineWidth(width[isam]);
      gr_noHot[ichip][isam]->SetMarkerStyle(MarkerStyleWheel(3));
      gr_noHot[ichip][isam]->Draw("PL");

      double critical=-1;
      if(rangeMax[isam][ichip]-rangeMin[isam][ichip]<=1){
	rangeMax[isam][ichip]=xmax>(rangeMin[isam][ichip]+2)?(rangeMin[isam][ichip]+2):xmax;
	cerr<<"Warning: "<<FE[ichip]<<" at "<<legtext[isam]<<" does not have enough points to perform fit. Extend the range to include three points"<<endl;
	cout<<rangeMin[isam][ichip]<<"->"<<rangeMax[isam][ichip]<<endl;
      }
      if(!option.Contains("nofit")){
	gr_noHot[ichip][isam]->Fit("expo","","",translateAltFine(rangeMin[isam][ichip], FE[ichip]), translateAltFine(rangeMax[isam][ichip], FE[ichip]));
	func[ichip][isam]=gr_noHot[ichip][isam]->GetFunction("expo");
	func[ichip][isam]->SetName("Expo_"+FE[ichip]+Form("_isamation%d", isam));
	
	// func[ichip][isam]=new TF1("GausCDF_C_"+FE[ichip]+Form("_isamation%d", isam), GausCDF_C, translateAltFine(rangeMin[isam][ichip], FE[ichip]), translateAltFine(rangeMax[isam][ichip], FE[ichip]), 3);
	// func[ichip][isam]->SetParName(0,"Mean");
	// func[ichip][isam]->SetParameter(0,100);
	// func[ichip][isam]->SetParLimits(0,-10,10);
	
	// func[ichip][isam]->SetParName(1,"Sigma");
	// func[ichip][isam]->SetParameter(1,100);
	// func[ichip][isam]->SetParLimits(1,0,10000);
	
	// func[ichip][isam]->SetParName(2,"Constant");
	// func[ichip][isam]->SetParameter(2,1);
	// func[ichip][isam]->SetParLimits(2,0,2);
	
	gr_noHot[ichip][isam]->Fit(func[ichip][isam],"","",translateAltFine(rangeMin[isam][ichip], FE[ichip]), translateAltFine(rangeMax[isam][ichip], FE[ichip]));
	func[ichip][isam]->SetLineColor(color[isam]);
	func[ichip][isam]->SetLineWidth(width[isam]);
      
	func[ichip][isam]->DrawClone("same");
	
	func[ichip][isam]->SetRange(xmin_th, xmax_th);
	func[ichip][isam]->SetLineStyle(3);
	func[ichip][isam]->DrawClone("same");
	critical=BinarySearch(func[ichip][isam], noiseThreshold, xmin_th, xmax_th, noiseThreshold/1e3);
	leg->AddEntry(gr_noHot[ichip][isam],legtext[isam].ReplaceAll("_"," ")+Form(": %.0fe (%.0f)", critical, translateThreshold(critical, FE[ichip])),legopt[isam]);
      
      }
      else leg->AddEntry(gr_noHot[ichip][isam],legtext[isam].ReplaceAll("_"," "),legopt[isam]);

    }
    TPaveText* text=CommonFunc::CreatePaveText(x1+0.05, y1+0.05, x1+0.35, y1+0.3, pavetext, 0.04);
    text->Draw("same");
    leg->Draw();
    
    c->RedrawAxis();
    TString output1DName=outputDir+"/NoiseVsAltFine"+FE[ichip];

    CommonFunc::PrintCanvas(c,output1DName);

    c->SetLogy();

    CommonFunc::PrintCanvas(c,output1DName+"_logy");
  }
}
