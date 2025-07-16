#include <iostream>
#include <vector>
#include <cmath>
#include "TMath.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TF3.h"
#include "TF2.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1F.h"
#include <thread>   // for threads
#include <chrono>   // for expresssing duration
#include <atomic>   
#include <iomanip>
#include <functional>   // for std::ref()
using namespace std; 

// Глобальные переменные для хранения данных
std::vector<double> Ns = {40, 80, 40};    // Сигнал
std::vector<double> sigma_s = {0.1, 0.05, 0.1};  // Относительные погрешности сигнала
std::vector<double> Nb = {60, 50, 40};    // Фон
std::vector<double> sigma_b = {0.05, 0.1, 0.15}; // Относительные погрешности фона
std::vector<double> Nobs = {100, 110, 90}; // Наблюдаемые события

// Функция предсказания числа событий
double Npred(int bin, double mu, double theta_s, double theta_b) {
    return mu * Ns[bin] * (1 + sigma_s[bin] * theta_s) + Nb[bin] * (1 + sigma_b[bin] * theta_b);
}

// Двойной отрицательный логарифм правдоподобия
double doubleNLL(double mu, double theta_s, double theta_b, double theta_s0, double theta_b0) {
    double logL = 0.0;
    for (size_t i = 0; i < Nobs.size(); ++i) {
        double Ni = Npred(i, mu, theta_s, theta_b);
        logL += Nobs[i] * log(Ni) - Ni; // Poisson part
    }
    logL += -0.5 * (theta_s - theta_s0) * (theta_s - theta_s0); // Гауссово ограничение на theta_s
    logL += -0.5 * (theta_b - theta_b0) * (theta_b - theta_b0); // Гауссово ограничение на theta_b
    return -2 * logL;
}

void GenerateToyMC(double mu, double theta_s, double theta_b, std::vector<double>& toy_data) {
    TRandom3 rnd(0);
    toy_data.clear();
    for (size_t i = 0; i < Nobs.size(); ++i) {
        double Ni = Npred(i, mu, theta_s, theta_b);
        toy_data.push_back(rnd.Poisson(Ni));
    }
}

// Функция для анализа одного бина
void AnalyzeSingleBin(int bin, double &mu_hat, double &mu_68_low, double &mu_68_low_toy,  double &mu_68_high, double &mu_68_high_toy, 
                    double &mu_95_low, double &mu_95_low_toy, double &mu_95_high, double &mu_95_high_toy,
                    TGraph* &graph_pvalue, TGraph* &graph_teststat, TH1D* &hist_teststat_mc_bin, 
                    TGraph* &graph_pvalue_bin_toy, TGraph* &graph_teststat_bin_toy) {
    // Сохраняем оригинальные данные
    std::vector<double> Ns_orig = Ns;
    std::vector<double> sigma_s_orig = sigma_s;
    std::vector<double> Nb_orig = Nb;
    std::vector<double> sigma_b_orig = sigma_b;
    std::vector<double> Nobs_orig = Nobs;
    
    // Оставляем только выбранный бин
    Ns = {Ns_orig[bin]};
    sigma_s = {sigma_s_orig[bin]};
    Nb = {Nb_orig[bin]};
    sigma_b = {sigma_b_orig[bin]};
    Nobs = {Nobs_orig[bin]};
    
    // Находим оценку максимального правдоподобия для mu
    TF3 *f = new TF3("f", "doubleNLL(x, y, z, [0], [1])", 0, 5, -5, 5, -5, 5);
    f->SetParameter(0, 0);
    f->SetParameter(1, 0);
    double gmin[3] = {0, 0, 0};
    double global_min = f->GetMinimum(gmin);
    
    mu_hat = gmin[0];
    double theta_s_hat = gmin[1];
    double theta_b_hat = gmin[2];
    
    std::cout << "\nResults for bin " << bin << ":" << std::endl;
    std::cout << "mu_hat = " << mu_hat << std::endl;
    std::cout << "theta_s_hat = " << theta_s_hat << std::endl;
    std::cout << "theta_b_hat = " << theta_b_hat << std::endl;
    
    // Создаем графики
    graph_pvalue = new TGraph();
    graph_teststat = new TGraph();
    hist_teststat_mc_bin = new TH1D("hist_teststat_mc_bin", 
                                        "Test statistic distribution (MC) for bin, #mu = 1;q_{#mu}", 
                                        100, 0, 10);
    graph_pvalue_bin_toy = new TGraph();
    
    double mu_step = 0.05;
    double mu_test = 1.;
    bool    low68_found = false,
            high68_found = false,
            low95_found = false,
            high95_found = false;
    bool    low68_found_toy = false,
            high68_found_toy = false,
            low95_found_toy = false,
            high95_found_toy = false;
    for (double mu = 0; mu < 3; mu += mu_step) {
        // Фиксируем mu и находим минимум
        TF2 *f_fixed_mu = new TF2("f_fixed_mu", "doubleNLL([0], x, y, [1], [2])", -5, 5, -5, 5);
        f_fixed_mu->SetParameter(0, mu);
        f_fixed_mu->SetParameter(1, 0);
        f_fixed_mu->SetParameter(2, 0);
        
        double min_params[2] = {0, 0};
        double min_val = f_fixed_mu->GetMinimum(min_params);

        TRandom3* r_data = new TRandom3(1);
        TRandom3* r_thetas0 = new TRandom3(2);
        TRandom3* r_thetab0 = new TRandom3(3);

        int n_toys = 1e4;
        bool draw_flag = true;
        double probability = 0;

        double qmu = min_val - global_min;
        double pvalue = (qmu <= 0) ? 1.0 : TMath::Prob(qmu, 1);
        if (pvalue > 0.05 && !low95_found){
            low95_found = true;
            mu_95_low = mu;
        }
        if (pvalue < 0.05 && !high95_found && low95_found){
            high95_found = true;
            mu_95_high = mu;
        }

        if (pvalue > 0.32 && !low68_found){
            low68_found = true;
            mu_68_low = mu;
        }
        if (pvalue < 0.32 && !high68_found && low68_found){
            high68_found = true;
            mu_68_high = mu;
        }

        if (abs(mu - mu_test) < 1e-6) {
            n_toys = 1e5;
            for (int i = 0; i < n_toys; ++i) {
                // Генерируем псевдоэксперимент
                std::vector<double> toy_data;
                GenerateToyMC(mu, min_params[0], min_params[1], toy_data);
                double thetas_0 = r_thetas0->Gaus(min_params[0], 1);
                double thetab_0 = r_thetab0->Gaus(min_params[1], 1);
                
                // Сохраняем оригинальные данные
                std::vector<double> Nobs_orig_temp = Nobs;
                
                // Подставляем данные псевдоэксперимента
                Nobs = toy_data;
                

                // Вычисляем тестовую статистику для псевдоэксперимента
                TF3 *f_toy = new TF3("f_toy", "doubleNLL(x, y, z, [0], [1])", 0, 5, -5, 5, -5, 5);
                f_toy->SetParameter(0, thetas_0);
                f_toy->SetParameter(1, thetab_0);
                double gmin_toy[3] = {0, 0, 0};
                double global_min_toy = f_toy->GetMinimum(gmin_toy);
                
                TF2 *f_fixed_mu_toy = new TF2("f_fixed_mu_toy", "doubleNLL([0], x, y, [1], [2])", -5, 5, -5, 5);
                f_fixed_mu_toy->SetParameter(0, mu);
                f_fixed_mu_toy->SetParameter(1, thetas_0);
                f_fixed_mu_toy->SetParameter(2, thetab_0);
                double min_val_toy = f_fixed_mu_toy->GetMinimum(min_params);
                
                double qmu_toy = min_val_toy - global_min_toy;
                hist_teststat_mc_bin->Fill(qmu_toy);
                
                // Восстанавливаем данные
                Nobs = Nobs_orig_temp;
                if(qmu < qmu_toy)
                    probability++;
            }
            draw_flag = false;
        }
        if (draw_flag){
            for (int i = 0; i < n_toys; ++i) {
                // Генерируем псевдоэксперимент
                std::vector<double> toy_data;
                GenerateToyMC(mu, min_params[0], min_params[1], toy_data);
                double thetas_0 = r_thetas0->Gaus(min_params[0], 1);
                double thetab_0 = r_thetab0->Gaus(min_params[1], 1);
                
                // Сохраняем оригинальные данные
                std::vector<double> Nobs_orig_temp = Nobs;
                
                // Подставляем данные псевдоэксперимента
                Nobs = toy_data;
                

                // Вычисляем тестовую статистику для псевдоэксперимента
                TF3 *f_toy = new TF3("f_toy", "doubleNLL(x, y, z, [0], [1])", 0, 5, -5, 5, -5, 5);
                f_toy->SetParameter(0, thetas_0);
                f_toy->SetParameter(1, thetab_0);
                double gmin_toy[3] = {0, 0, 0};
                double global_min_toy = f_toy->GetMinimum(gmin_toy);
                
                TF2 *f_fixed_mu_toy = new TF2("f_fixed_mu_toy", "doubleNLL([0], x, y, [1], [2])", -5, 5, -5, 5);
                f_fixed_mu_toy->SetParameter(0, mu);
                f_fixed_mu_toy->SetParameter(1, thetas_0);
                f_fixed_mu_toy->SetParameter(2, thetab_0);
                double min_val_toy = f_fixed_mu_toy->GetMinimum(min_params);
                
                double qmu_toy = min_val_toy - global_min_toy;
                
                // Восстанавливаем данные
                Nobs = Nobs_orig_temp;

                if(qmu < qmu_toy)
                    probability++;
            }
        }
        
        probability /= n_toys;

        if (probability > 0.05 && !low95_found_toy){
            low95_found_toy = true;
            mu_95_low_toy = mu;
        }
        if (probability < 0.05 && !high95_found_toy && low95_found_toy){
            high95_found_toy = true;
            mu_95_high_toy = mu;
        }

        if (probability > 0.32 && !low68_found_toy){
            low68_found_toy = true;
            mu_68_low_toy = mu;
        }
        if (probability < 0.32 && !high68_found_toy && low68_found_toy){
            high68_found_toy = true;
            mu_68_high_toy = mu;
        }
        graph_pvalue->SetPoint(graph_pvalue->GetN(), mu, pvalue);
        graph_teststat->SetPoint(graph_teststat->GetN(), mu, qmu);
        graph_pvalue_bin_toy->SetPoint(graph_pvalue_bin_toy->GetN(), mu, probability);


    }
    
    // Находим доверительные интервалы
    
    
    std::cout << "68% CL: [" << mu_68_low << ", " << mu_68_high << "]" << std::endl;
    std::cout << "95% CL: [" << mu_95_low << ", " << mu_95_high << "]" << std::endl;
    
    // Восстанавливаем оригинальные данные
    Ns = Ns_orig;
    sigma_s = sigma_s_orig;
    Nb = Nb_orig;
    sigma_b = sigma_b_orig;
    Nobs = Nobs_orig;
}

int main() {
    // Анализ для всех бинов вместе
    std::cout << "=== Analysis using all bins ===" << std::endl;
    TF3 *f_all = new TF3("f_all", "doubleNLL(x, y, z, [0], [1])", 0, 5, -5, 5, -5, 5);
    f_all->SetParameter(0, 0);
    f_all->SetParameter(1, 0);
    double gmin_all[3] = {0, 0, 0};
    double global_min_all = f_all->GetMinimum(gmin_all);
    
    double mu_hat_all = gmin_all[0];
    std::cout << "Global mu_hat = " << mu_hat_all << std::endl;
    
    // Графики для всех бинов
    TGraph *graph_pvalue_all = new TGraph();
    TGraph *graph_teststat_all = new TGraph();
    TH1D *hist_teststat_mc_all = new TH1D("hist_teststat_mc_all", 
                                        "Test statistic distribution (MC) for all bins, #mu = 1;t^{obs}_{#mu}", 
                                        100, 0, 10);

    TGraph *graph_pvalue_all_toy = new TGraph();
    TGraph *graph_teststat_all_toy = new TGraph();
    
    double mu_step = 0.05;
    double mu_test = 1.0;

    double mu_low95, mu_high95, mu_low68, mu_high68;
    double mu_low95_toy, mu_high95_toy, mu_low68_toy, mu_high68_toy;
    bool    low68_found = false,
            high68_found = false,
            low95_found = false,
            high95_found = false;
    bool    low68_found_toy = false,
            high68_found_toy = false,
            low95_found_toy = false,
            high95_found_toy = false;
    for (double mu = 0.; mu < 1.5; mu += mu_step) {
        TF2 *f_fixed_mu = new TF2("f_fixed_mu", "doubleNLL([0], x, y, [1], [2])", -5, 5, -5, 5);
        f_fixed_mu->SetParameter(0, mu);
        f_fixed_mu->SetParameter(1, 0);
        f_fixed_mu->SetParameter(2, 0);
        
        double min_params[2] = {0, 0};
        double min_val = f_fixed_mu->GetMinimum(min_params);

        TRandom3* r_data = new TRandom3(1);
        TRandom3* r_thetas0 = new TRandom3(2);
        TRandom3* r_thetab0 = new TRandom3(3);

        bool draw_flag = true;
        int n_toys = 1e4;
        double probability = 0;
        double qmu = (min_val - global_min_all);

        double pvalue = (qmu <= 0) ? 1.0 : TMath::Prob(qmu, 1);
        if (pvalue > 0.05 && !low95_found){
            low95_found = true;
            mu_low95 = mu;
        }
        if (pvalue < 0.05 && !high95_found && low95_found){
            high95_found = true;
            mu_high95 = mu;
        }

        if (pvalue > 0.32 && !low68_found){
            low68_found = true;
            mu_low68 = mu;
        }
        if (pvalue < 0.32 && !high68_found && low68_found){
            high68_found = true;
            mu_high68 = mu;
        }

        if (abs(mu - mu_test) < 1e-6) {
            n_toys = 1e5;
            for (int i = 0; i < n_toys; ++i) {
                // Генерируем псевдоэксперимент
                std::vector<double> toy_data;
                GenerateToyMC(mu, min_params[0], min_params[1], toy_data);
                double thetas_0 = r_thetas0->Gaus(min_params[0], 1);
                double thetab_0 = r_thetab0->Gaus(min_params[1], 1);
                
                // Сохраняем оригинальные данные
                std::vector<double> Nobs_orig_temp = Nobs;
                
                // Подставляем данные псевдоэксперимента
                Nobs = toy_data;
                
                // Вычисляем тестовую статистику для псевдоэксперимента
                TF3 *f_toy = new TF3("f_toy", "doubleNLL(x, y, z, [0], [1])", 0, 5, -5, 5, -5, 5);
                f_toy->SetParameter(0, thetas_0);
                f_toy->SetParameter(1, thetab_0);
                double gmin_toy[3] = {0, 0, 0};
                double global_min_toy = f_toy->GetMinimum(gmin_toy);
                
                TF2 *f_fixed_mu_toy = new TF2("f_fixed_mu_toy", "doubleNLL([0], x, y, [1], [2])", -5, 5, -5, 5);
                f_fixed_mu_toy->SetParameter(0, mu);
                f_fixed_mu_toy->SetParameter(1, thetas_0);
                f_fixed_mu_toy->SetParameter(2, thetab_0);
                double min_val_toy = f_fixed_mu_toy->GetMinimum(min_params);
                
                double qmu_toy = min_val_toy - global_min_toy;
                hist_teststat_mc_all->Fill(qmu_toy);

                if (qmu_toy > qmu)
                    probability++;
                
                // Восстанавливаем данные
                Nobs = Nobs_orig_temp;
            }
            draw_flag = false;
        }

        if (draw_flag){
            for (int i = 0; i < n_toys; ++i) {
                // Генерируем псевдоэксперимент
                std::vector<double> toy_data;
                GenerateToyMC(mu, min_params[0], min_params[1], toy_data);
                double thetas_0 = r_thetas0->Gaus(min_params[0], 1);
                double thetab_0 = r_thetab0->Gaus(min_params[1], 1);
                
                // Сохраняем оригинальные данные
                std::vector<double> Nobs_orig_temp = Nobs;
                
                // Подставляем данные псевдоэксперимента
                Nobs = toy_data;
                
                // Вычисляем тестовую статистику для псевдоэксперимента
                TF3 *f_toy = new TF3("f_toy", "doubleNLL(x, y, z, [0], [1])", 0, 5, -5, 5, -5, 5);
                f_toy->SetParameter(0, thetas_0);
                f_toy->SetParameter(1, thetab_0);
                double gmin_toy[3] = {0, 0, 0};
                double global_min_toy = f_toy->GetMinimum(gmin_toy);
                
                TF2 *f_fixed_mu_toy = new TF2("f_fixed_mu_toy", "doubleNLL([0], x, y, [1], [2])", -5, 5, -5, 5);
                f_fixed_mu_toy->SetParameter(0, mu);
                f_fixed_mu_toy->SetParameter(1, thetas_0);
                f_fixed_mu_toy->SetParameter(2, thetab_0);
                double min_val_toy = f_fixed_mu_toy->GetMinimum(min_params);
                
                double qmu_toy = min_val_toy - global_min_toy;
                
                // Восстанавливаем данные
                if (qmu_toy > qmu)
                    probability++;
                Nobs = Nobs_orig_temp;
            }
        }

        probability /= n_toys;
        if (probability > 0.05 && !low95_found_toy){
            low95_found_toy = true;
            mu_low95_toy = mu;
        }
        if (probability < 0.05 && !high95_found_toy && low95_found_toy){
            high95_found_toy = true;
            mu_high95_toy = mu;
        }

        if (probability > 0.32 && !low68_found_toy){
            low68_found_toy = true;
            mu_low68_toy = mu;
        }
        if (probability < 0.32 && !high68_found_toy && low68_found_toy){
            high68_found_toy = true;
            mu_high68_toy = mu;
        }

        /*std::cout << qmu << std::endl;
        qmu /= global_min_all;
        std::cout << "div:" << qmu << std::endl;*/
        
        
        graph_pvalue_all->SetPoint(graph_pvalue_all->GetN(), mu, pvalue);
        graph_teststat_all->SetPoint(graph_teststat_all->GetN(), mu, abs(qmu));
        graph_pvalue_all_toy->SetPoint(graph_pvalue_all_toy->GetN(), mu, probability);
    }

    std::cout << "68% CL: [" << mu_low68 << ", " << mu_high68 << "]" << std::endl;
    std::cout << "95% CL: [" << mu_low95 << ", " << mu_high95 << "]" << std::endl;
    
    // Анализ для каждого бина в отдельности
    double mu_hat_bin[3], mu_68_low_bin[3], mu_68_low_bin_toy[3], 
    mu_68_high_bin[3], mu_68_high_bin_toy[3], mu_95_low_bin[3], mu_95_low_bin_toy[3], 
    mu_95_high_bin[3], mu_95_high_bin_toy[3];
    TGraph *graph_pvalue_bin[3];
    TGraph *graph_teststat_bin[3];
    TGraph *graph_pvalue_bin_toy[3];
    TGraph *graph_teststat_bin_toy[3];
    TH1D *hist_teststat_mc_bin[3];
    
    for (int bin = 0; bin < 3; ++bin) {
        AnalyzeSingleBin(bin, mu_hat_bin[bin], mu_68_low_bin[bin], mu_68_low_bin_toy[bin], 
                        mu_68_high_bin[bin], mu_68_high_bin_toy[bin], mu_95_low_bin[bin], 
                        mu_95_low_bin_toy[bin], mu_95_high_bin[bin], mu_95_high_bin_toy[bin], 
                        graph_pvalue_bin[bin], graph_teststat_bin[bin], hist_teststat_mc_bin[bin], 
                        graph_pvalue_bin_toy[bin], graph_teststat_bin_toy[bin]);
    }

    std::cout << "============RESULTS============" << std::endl;
    for (int i = -1; i < 3; ++i){
        if (i < 0){
            std::cout << "AllBins:\n mu_hat = " << mu_hat_all << std::endl;
            std::cout << "68% CL: [" << mu_low68 << ", " << mu_high68 << "]" << "toy: [" << mu_low68_toy << ", " << mu_high68_toy << "]" << std::endl;
            std::cout << "95% CL: [" << mu_low95 << ", " << mu_high95 << "]" << "toy: [" << mu_low95_toy << ", " << mu_high95_toy << "]" << std::endl;
        }
        else{
            std::cout << "Bin " << i << ":\n mu_hat = " << mu_hat_bin[i] << std::endl;
            std::cout << "68% CL: [" << mu_68_low_bin[i] << ", " << mu_68_high_bin[i] << "]" 
            << "toy: [" << mu_68_low_bin_toy[i] << ", " << mu_68_high_bin_toy[i] << "]" << std::endl;
            std::cout << "95% CL: [" << mu_95_low_bin[i] << ", " << mu_95_high_bin[i] << "]" 
            << "toy: [" << mu_95_low_bin_toy[i] << ", " << mu_95_high_bin_toy[i] << "]" << std::endl;
        }
    }
    
    // Сохранение результатов в файл
    TFile *outfile = new TFile("results.root", "RECREATE");
    
    // Записываем результаты для всех бинов
    graph_pvalue_all->Write("pvalue_all");
    graph_teststat_all->Write("teststat_all");
    
    // Записываем результаты для отдельных бинов
    for (int bin = 0; bin < 3; ++bin) {
        graph_pvalue_bin[bin]->Write(Form("pvalue_bin%d", bin));
        graph_teststat_bin[bin]->Write(Form("teststat_bin%d", bin));
    }
    
    outfile->Close();
    
    // Визуализация результатов
    TCanvas *c1 = new TCanvas("c1", "p-values", 1200, 800);
    c1->Divide(2, 2);
    
    // Все бины
    c1->cd(1)->SetLogy();
    TMultiGraph* pValAll = new TMultiGraph();
    graph_pvalue_all->SetTitle("All bins;#mu;p-value");
    graph_pvalue_all->SetLineWidth(2);
    graph_pvalue_all->SetLineColor(kBlack);

    graph_pvalue_all_toy->SetLineWidth(4);
    graph_pvalue_all_toy->SetLineColor(kRed);

    pValAll->SetTitle("All bins;#mu;p-value");
    pValAll->Add(graph_pvalue_all);
    pValAll->Add(graph_pvalue_all_toy);
    pValAll->Draw("AL");
    
    // Отдельные бины
    c1->cd(2)->SetLogy();
    TMultiGraph *mg_pvalue = new TMultiGraph();
    mg_pvalue->SetTitle("Single bins;#mu;p-value");
    
    for (int bin = 0; bin < 3; ++bin) {
        graph_pvalue_bin[bin]->SetLineWidth(2);
        graph_pvalue_bin_toy[bin]->SetLineWidth(4);
        graph_pvalue_bin[bin]->SetLineColor(bin+1);
        graph_pvalue_bin_toy[bin]->SetLineColor(bin+6);
        mg_pvalue->Add(graph_pvalue_bin[bin]);
        mg_pvalue->Add(graph_pvalue_bin_toy[bin]);
    }
    
    mg_pvalue->Draw("AL");
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(graph_pvalue_bin[0], "Bin 0", "l");
    leg->AddEntry(graph_pvalue_bin[1], "Bin 1", "l");
    leg->AddEntry(graph_pvalue_bin[2], "Bin 2", "l");
    leg->AddEntry(graph_pvalue_bin_toy[0], "Bin 0 toy", "l");
    leg->AddEntry(graph_pvalue_bin_toy[1], "Bin 1 toy", "l");
    leg->AddEntry(graph_pvalue_bin_toy[2], "Bin 2 toy", "l");
    leg->Draw();
    
    // Тестовая статистика для всех бинов
    c1->cd(3);
    graph_teststat_all->SetTitle("All bins;#mu;t^{obs}_{#mu}");
    graph_teststat_all->SetLineWidth(2);
    graph_teststat_all->SetLineColor(kBlack);
    graph_teststat_all->Draw("AL");
    
    // Тестовая статистика для отдельных бинов
    c1->cd(4);
    TMultiGraph *mg_teststat = new TMultiGraph();
    mg_teststat->SetTitle("Single bins;#mu;t^{obs}_{#mu}");
    
    for (int bin = 0; bin < 3; ++bin) {
        graph_teststat_bin[bin]->SetLineWidth(2);
        graph_teststat_bin[bin]->SetLineColor(bin+1);
        mg_teststat->Add(graph_teststat_bin[bin]);
    }
    
    mg_teststat->Draw("AL");
    leg->Draw();
    
    c1->SaveAs("results.pdf");

    auto createChi2Graph = [](double xmin, double xmax, double ymax) {
        int n_points = 1000;
        TGraph* chi2_graph = new TGraph(n_points);
        for (int i = 0; i < n_points; ++i) {
            double x = xmin + (xmax - xmin) * i / (n_points - 1);
            double y = TMath::Prob(x, 1); // 1 степень свободы
            chi2_graph->SetPoint(i, x, y * ymax); // Масштабируем под гистограмму
        }
        chi2_graph->SetLineColor(kRed);
        chi2_graph->SetLineWidth(2);
        chi2_graph->SetTitle(";Test statistic;Normalized counts / probability");
        return chi2_graph;
    };

    // Создаем холст для сравнения асимптотического и MC распределений
    TCanvas *c2 = new TCanvas("c2", "Test statistic comparison", 1200, 800);
    c2->Divide(2, 2);
    
    // Для всех бинов
    c2->cd(1)->SetLogy();
    hist_teststat_mc_all->SetLineColor(kBlue);
    hist_teststat_mc_all->Draw("PE");

    TGraph* chi2_graph_all = createChi2Graph(
        hist_teststat_mc_all->GetXaxis()->GetXmin(),
        hist_teststat_mc_all->GetXaxis()->GetXmax(),
        1e4
    );
    chi2_graph_all->Draw("L SAME");
    
    TLegend *leg2 = new TLegend(0.5, 0.5, 0.7, 0.7);
    leg2->AddEntry(hist_teststat_mc_all, "Monte Carlo", "l");
    leg2->AddEntry(chi2_graph_all, "#chi^{2} (1 dof)", "l");
    leg2->Draw();
    
    // Для отдельных бинов
    for (int bin = 0; bin < 3; ++bin) {
        c2->cd(bin + 2)->SetLogy();
        hist_teststat_mc_bin[bin]->SetLineColor(kBlue);
        hist_teststat_mc_bin[bin]->Draw("PE");

        TGraph* chi2_graph_bin = createChi2Graph(
            hist_teststat_mc_bin[bin]->GetXaxis()->GetXmin(),
            hist_teststat_mc_bin[bin]->GetXaxis()->GetXmax(),
            1e4
        );
    chi2_graph_bin->Draw("L SAME");
        
        TLegend *leg_bin = new TLegend(0.5, 0.5, 0.7, 0.7);
        leg_bin->AddEntry(hist_teststat_mc_bin[bin], "Monte Carlo", "l");
        leg_bin->AddEntry(chi2_graph_bin, "#chi^{2} (1 dof)", "l");
        leg_bin->Draw();
    }
    
    c2->SaveAs("teststat_comparison.pdf");
    
    return 0;
}