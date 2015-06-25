#include "CommonHead.h"
#include "CommonFunc.h"

using namespace std;
using namespace CommonFunc;

void PrintTemperatureProfile (const char *fname_in, const char *fname_out) {
    //
    // Reads a file outputted by HP_MultiMeter_Connect.exe on atlas2 and
    // outputs a TGraph of the temperature readings versus time.
    //
    // inputs: fname_in:  path to *.txt input file written
    //                    by HP_MultiMeter_Connect.exe
    //         fname_out: desired path to *.root output file;
    //                    if file exists, error is printed to
    //                    cerr and the function returns without
    //                    printing
    //
    // Could easily be edited to save *.pdf or *.eps files as well.
    //
    ifstream data_file(fname_in);

    if (!data_file.good()) {
        cerr << "Error in <PrintTemperatureProfile>: Cannot open input file "
             << fname_in << endl;
        return;
    }

    if (data_file.peek() == ifstream::traits_type::eof()) {
        cout << "No events present in input file "
             << fname_in << endl;
        return;
    }

    vector<float> temperatures;
    temperatures.clear();

    vector<TDatime> datimes;
    datimes.clear();

    string reading;
    while (getline(data_file, reading)) {
        float temp_crt;
        int h, m, s;
	int day;
        // TODO: edit leading trash to fit spec - might not currently fit
        int n_reads = sscanf(reading.c_str(), "%*s %*f %*s Temperature: %f K \tTime: %*s %*s %d %d:%d:%d", &temp_crt, &day, &h, &m, &s);
	// cout<<n_reads<<" "<<temp_crt<<" "<<h<<" "<<m<<" "<<s<<endl;
        if (n_reads != 5) continue; // end of file (trailing whitespace)
	cout<<temp_crt<<endl;
        // slow, but easier than manually dealing with the offset of time_t
        TDatime *datime_crt = new TDatime(1995, 1, day, h, m, s);
	temp_crt-=273.15;	// Convert it to Celsius
        temperatures.push_back(temp_crt);
        datimes.push_back(*datime_crt);
    }

    int n_entries = datimes.size();

    float *temp_array = (float *) &temperatures[0];
    float time_vals[n_entries];
    for (int i = 0; i < n_entries; ++i){
      time_vals[i] = datimes[i].Convert();
      cout<<time_vals[i]<<endl;
    }

    TGraph *h_temp_vs_time = new TGraph(n_entries, time_vals, temp_array);
    h_temp_vs_time->SetTitle("chip temperature during testing");
    h_temp_vs_time->SetMinimum(20);
    h_temp_vs_time->SetMaximum(30);
    h_temp_vs_time->GetXaxis()->SetTimeDisplay(1);
    h_temp_vs_time->GetXaxis()->SetTimeFormat("%d %H:%M");

    h_temp_vs_time->GetXaxis()->SetTitle("time (day min:sec)");
    h_temp_vs_time->GetYaxis()->SetTitle("temperature (C)");
    h_temp_vs_time->GetXaxis()->SetNdivisions(505);
    TCanvas *c = new TCanvas("c", "temperature vs. time", 0, 0, 800, 600);
    c->SetBorderMode(0);

    h_temp_vs_time->Draw("AP");

    system("mkdir -vp fig/temperatureLog/");
    PrintCanvas(c,TString("fig/temperatureLog/")+fname_out);
    // TFile *out_file = new TFile(fname_out, "RECREATE");

    // if (!out_file->IsOpen()) {
    //     cerr << "Error in <PrintTemperatureProfile>: Cannot open output file "
    //          << fname_out << endl;
    //     return;
    // }

    // h_temp_vs_time->Write();
}

// enable running outside of ROOT if, for whatever reason, this is desirable
#ifndef __CINT__

void standalone_application (int argc, char** argv) {
    if (argc != 3) gSystem->Exit(1);
    PrintTemperatureProfile(argv[1], argv[2]);
    gSystem->Exit(0);
}

int main (int argc, char** argv) {
    if ((argc == 2 && !strcmp(argv[1], "-h")) || argc != 3) {
        cout << "Usage: ./PrintTemperatureProfile fname_in.txt fname_out.root"
             << endl;
        exit(0);
    }

    TApplication app("app", &argc, argv);
    standalone_application(app.Argc(), app.Argv());
    app.Run();
 
    return 0;
}

#endif
