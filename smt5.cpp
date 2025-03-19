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
        double x = rand->Uniform(0, TMath::Pi()/4);
        double y = rand->Uniform(0, 1);

        if (y < TMath::Cos(2*x))
            hAcc->Fill(x);
    }
    TCanvas *c1 = new TCanvas("c1", "Monte Carlo acception method", 1200, 400);
    hAcc->Draw();
    c1->SaveAs("monte_carlo_acceptance.pdf");

}

void monte_carlo_reversed(){
    TRandom3 *rand = new TRandom3(0);
    
    TH1D *hRev = new TH1D("h2", "Monte Carlo reversed method cos(2x)", 100, 0, TMath::Pi()/2);

    double n_samples = 100000;

    for (size_t i = 0; i < n_samples; ++i){
        double r = rand->Uniform(0, 1);
        double x = - acos(r) / 2.0 + TMath::Pi()/4; 
        hRev->Fill(x);
    }
    TCanvas *c2 = new TCanvas("c2", "Monte Carlo acception method", 1200, 400);
    hRev->Draw();
    c2->SaveAs("monte_carlo_reverse.pdf");
}

int main(){
    monte_carlo_acception();
    monte_carlo_reversed();
    return 0;
}

