#include <TRandom.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>

void central_limit_theorem_check() {
    //modelation parameters
    const int nExperiments = 10000;
    const int nSamples = 1000;
    const double lambda = 1.0;

    //Hists declaration
    TH1D *hSum = new TH1D("hSum", "Distribution of Sums", 100, 0, 0);
    TH1D *hNorm = new TH1D("hNorm", "Normalized Distribution of Sums", 100, -5, 5);

    //Randomizer
    TRandom *rand = new TRandom();

    for (int i = 0; i < nExperiments; ++i) {
        double sum = 0.0;
        for (int j = 0; j < nSamples; ++j) {
            sum += rand->Exp(lambda); //exponential distribution generation
        }

        //mean and variance calculation
        double mean = nSamples / lambda;
        double variance = nSamples / (lambda * lambda);

        //Hists Filling
        hSum->Fill(sum);
        hNorm->Fill((sum - mean) / TMath::Sqrt(variance));
    }

    // Drawing
    TCanvas *c1 = new TCanvas("c1", "Central Limit Theorem", 800, 400);
    c1->Divide(2, 1);

    c1->cd(1);
    hSum->Draw();
    hSum->SetXTitle("Sum of Random Variables");
    hSum->SetYTitle("Counts");

    c1->cd(2);
    hNorm->Draw();
    hNorm->SetXTitle("Normalized Sum");
    hNorm->SetYTitle("Counts");

    // Сохранение в файл
    c1->SaveAs("central_limit_theorem_check.pdf");
}