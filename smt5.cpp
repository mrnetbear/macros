#include <TRandom3.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>

void monte_carlo_acception(){

    TRandom3 *rand = new TRandom3(0);
    
    TH1D *hAcc = new TH1D("h1", "Monte Carlo acception method cos(2x)", 100, 0, TMath::Pi()/2);

    double n_samples = 1000000;

    for (size_t i = 0; i < n_samples; ++i){
        double x = rand->Uniform(0, TMath::Pi()/2);
        double y = rand->Uniform(0, 1);

        if (y < TMath::Cos(2*x))
            hAcc->Fill(x);
    }
    TCanvas *c1 = new TCanvas("c1", "Monte Carlo acception method", 1200, 400);
    hAcc->Draw();
    c1->SaveAs("monte_carlo_acceptance.pdf");

}