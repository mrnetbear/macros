#include <TRandom.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>

void weight_measurements(){
    TRandom *rand = new TRandom();
    const int nSamples = 10000;
    const double    var1 = 0.15,
                    var2 = 0.35,
                    var3 = 0.20;
    const double    mean = 0.5;

    TH1D *h1 = new TH1D("h1", "First Hist, #mu = 0.5, #sigma^{2} = 0.15 ", 100, -1.0, 2.0);
    TH1D *h2 = new TH1D("h2", "Second Hist, #mu = 0.5, #sigma^{2} = 0.35 ", 100, -1.0, 2.0);
    TH1D *h3 = new TH1D("h3", "Third Hist, #mu = 0.5, #sigma^{2} = 0.20 ", 100, -1.0, 2.0);
    TH1D *h4 = new TH1D("h4", "Weighted Hist", 100, -1.0, 2.0);

    TCanvas *c1 = new TCanvas("c1", "Gauss distribution with different accuracy", 1200, 1000);
    c1->Divide(2, 2);

    for (int i = 0; i < nSamples; ++i) {
        double x1 = rand->Gaus(mean, TMath::Sqrt(var1));
        double x2 = rand->Gaus(mean, TMath::Sqrt(var2));
        double x3 = rand->Gaus(mean, TMath::Sqrt(var3));

        double w1 = 1.0 / var1;
        double w2 = 1.0 / var2;
        double w3 = 1.0 / var3;
        double w = w1 + w2 + w3;
        double weighted_x = (w1 * x1 + w2 * x2 + w3 * x3) / w;

        h1->Fill(x1);
        h2->Fill(x2);
        h3->Fill(x3);
        h4->Fill(weighted_x);
    }
    c1->cd(1);
    h1->Draw();
    c1->cd(2);
    h2->Draw();
    c1->cd(3);
    h3->Draw();
    c1->cd(4);
    h4->Draw();

    c1->SaveAs("weight_measurements.pdf");
}

void weighted_mean_measurements(){
    TRandom *rand = new TRandom();
    const int nSamples = 10000;
    const double    var1 = 0.05,
                    var2 = 0.10;
    const double    mean1 = 0.49,
                    mean2 = 0.55;

    TH1D *h11 = new TH1D("h1", "First Hist, #mu = 0.49, #sigma^{2} = 0.05 ", 100, -1.0, 2.0);
    TH1D *h12 = new TH1D("h2", "Second Hist, #mu = 0.55, #sigma^{2} = 0.10 ", 100, -1.0, 2.0);
    TH1D *h13 = new TH1D("h3", "Weighted Mean Hist", 100, -1.0, 2.0);

    TCanvas *c2 = new TCanvas("c2", "Weighted mean distribution with different accuracy", 1200, 400);
    c2->Divide(3, 1);

    for (int i = 0; i < nSamples; ++i) {
        double x1 = rand->Gaus(mean1, TMath::Sqrt(var1));
        double x2 = rand->Gaus(mean2, TMath::Sqrt(var2));

        double w1 = 1.0 / var1;
        double w2 = 1.0 / var2;
        double w = w1 + w2;
        double weighted_mean = (w1 * x1 + w2 * x2) / w;

        h11->Fill(x1);
        h12->Fill(x2);
        h13->Fill(weighted_mean);
    }

    c2->cd(1);
    h11->Draw();
    c2->cd(2);
    h12->Draw();
    c2->cd(3);
    h13->Draw();
    c2->SaveAs("weighted_mean_measurements.pdf");
}

int main(){
    weight_measurements();
    weighted_mean_measurements();
    return 0; 
}