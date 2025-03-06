#include <TRandom.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

void shov_filter(std::vector<double>& data){
    double mean = 0;
    double stddev = 0;
    for (double val : data) {
        mean += val;
    }
    mean /= data.size();
    for (double val : data) {
        stddev += pow(val - mean, 2);
    }
    stddev = sqrt(stddev / data.size());

    auto it = data.begin();
    while (it != data.end()) {
        double z = fabs((*it - mean) / stddev);
        double p = 2 * (1 - TMath::Freq(z));
        if (p < 1.0 / (2 * data.size())) {
            it = data.erase(it);
        } else {
            ++it;
        }
    }
}
void shov_check(){
    short n = 50;
    std::vector<double> x;

    double  mean = 0.0,
            meanSquare = 0.0;
    std::fstream file("chi2task_v4.dat");
    TH1D *h1 = new TH1D("h1", "Experimental distribution, no filters", 70, 80, 120);
    TH1D *h2 = new TH1D("h2", "Experimental distribution + shov", 70, 80, 120);

    TCanvas *c1 = new TCanvas("c1", "Experimental distributions", 1200, 400); 
    for (size_t i = 0; i < n; ++i){
        std::string line;
        getline(file, line);
        std::istringstream iss(line);
        double value;
        iss >> value;
        x.push_back(value);
        mean += value;
        meanSquare += value * value;
        h1->Fill(value);
    }
    file.close();

    sort(x.begin(), x.end());
    mean /= n;
    double var = 0.0;

    for (size_t i = 0; i < 70; ++i){
        double value = h1->GetBinCenter(i);
        double weight  = h1->GetBinContent(i);
        var += (value - mean) * (value - mean) * weight / 50.0;
    }
    
    std::cout << "Mean = " << mean << "; Var = " << var << std::endl;
    
    TF1 *gaus = new TF1("gaus", "gaus", 80, 120);
    gaus->SetParameters(1000, mean, TMath::Sqrt(var));
    h1->Fit("gaus", "Q");

    double chi2 = 0.0;
    int ndf = 0;
    for (size_t i = 0; i < h1->GetNbinsX(); ++i){
        double observed = h1->GetBinCenter(i);
        double expected = gaus->Eval(h1->GetBinCenter(i));
        if (expected){
            chi2 += (observed - expected) * (observed - expected) / expected;
            ndf++;
        }
    }

    ndf -=3;
    std::cout << "Chi2 = " << chi2 << std::endl;
    std::cout << "NDF = " << ndf << std::endl;
    std::cout << "Chi2/ndf = " << (double) chi2/ndf << std::endl;

    bool shovFlag = true;

    while (shovFlag){
        int n = x.size();
        shov_filter(x);
        if (n > x.size())
            shovFlag = true;
        else
            shovFlag = false;

    }

    mean = 0.0;
    for (size_t i = 0; i < x.size(); ++i){
        h2->Fill(x[i]);
        mean += x[i];
    }

    mean /= n;
    var = 0.0;

    for (size_t i = 0; i < h2->GetNbinsX(); ++i){
        double value = h2->GetBinCenter(i);
        double weight  = h2->GetBinContent(i);
        var += (value - mean) * (value - mean) * weight / x.size();
    }

    std::cout << "Mean = " << mean << "; Var = " << var << std::endl;

    TF1 *gaus1 = new TF1("gaus1", "gaus", 80, 120);

    gaus1->SetParameters(1000, mean, TMath::Sqrt(var));
    h2->Fit("gaus", "Q");

    chi2 = 0.0;
    ndf = 0;

    for (size_t i = 0; i < h2->GetNbinsX(); ++i){
        double observed = h2->GetBinCenter(i);
        double expected = gaus->Eval(h2->GetBinCenter(i));
        if (expected){
            chi2 += (observed - expected) * (observed - expected) / expected;
            ndf++;
        }
    }

    ndf -=3;
    std::cout << "Chi2 = " << chi2 << std::endl;
    std::cout << "NDF = " << ndf << std::endl;
    std::cout << "Chi2/ndf = " << (double) chi2/ndf << std::endl;

    c1->Divide(2, 1);
    c1->cd(1);
    h1->Draw();
    c1->cd(2);
    h2->Draw();
    c1->SaveAs("shov_check.pdf");


}