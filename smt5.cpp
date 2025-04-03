#include <TRandom3.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
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

    double n_samples = 1000000;

    for (size_t i = 0; i < n_samples; ++i){
        double r = rand->Uniform(0, 1);
        double x = - acos(r) / 2.0 + TMath::Pi()/4; 
        hRev->Fill(x);
    }
    TF1* frev = new TF1("frev", "cos(2*x)*31000", 0, TMath::Pi()/4);
    TCanvas *c2 = new TCanvas("c2", "Monte Carlo acception method", 1200, 400);
    hRev->Draw();
    frev->Draw("samel");
    c2->SaveAs("monte_carlo_reverse.pdf");
}

void monte_carlo_acception_new(){
    TRandom3 *rand = new TRandom3(0);

    size_t n_samples = 1000000;
    std::vector<double> x, y;
    x.reserve(n_samples);
    y.reserve(n_samples);

    size_t attempts = 0;
    const size_t max_attempts = n_samples * 10; // Limit total attempts

    while (x.size() < n_samples && attempts < max_attempts) {
        double x_t = rand->Uniform(0, TMath::Pi()/4);
        double y_t = rand->Uniform(0, 1);

        if (y_t < TMath::Cos(2*x_t)){
            x.push_back(x_t);
            y.push_back(y_t);
        }
        attempts++;
        if (attempts % 10000 == 0) {
            std::cout << "Attempts: " << attempts << ", Accepted: " << x.size() << std::endl;
        }
    }

    if (x.size() < n_samples) {
        std::cout << "Warning: Only generated " << x.size() << " samples out of " << n_samples << " requested." << std::endl;
    }

    std::cout << "Sampling complete. Attempts: " << attempts << ", Accepted: " << x.size() << std::endl;

    TGraph *gAcc = new TGraph(x.size(), x.data(), y.data());
    TCanvas *c3 = new TCanvas("c3", "Monte Carlo acceptance method", 1200, 400);
    TF1 *cos2x = new TF1("cos2x", "cos(2*x)", 0, TMath::Pi()/4);
    gAcc->Draw("AP");
    cos2x->Draw("same");
    c3->SaveAs("monte_carlo_acceptance_new.pdf");
}

int main(){
    monte_carlo_acception();
    monte_carlo_reversed();
    monte_carlo_acception_new();
    return 0;
}

