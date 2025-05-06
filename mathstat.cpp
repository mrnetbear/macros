#include <iostream>
#include <vector>
#include <cmath>
#include "TMath.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"

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

// Логарифм правдоподобия (Poisson + гауссовы ограничения на theta)
double LogLikelihood(double mu, double theta_s, double theta_b) {
    double logL = 0.0;
    for (size_t i = 0; i < Nobs.size(); ++i) {
        double Ni = Npred(i, mu, theta_s, theta_b);
        logL += Nobs[i] * log(Ni) - Ni; // Poisson part (константа опущена)
    }
    logL += -0.5 * theta_s * theta_s; // Гауссово ограничение на theta_s
    logL += -0.5 * theta_b * theta_b; // Гауссово ограничение на theta_b
    return logL;
}

// Обёртка для Minuit (нужна для минимизации)
void FCN(int &npar, double *gin, double &f, double *par, int iflag) {
    double mu = par[0];
    double theta_s = par[1];
    double theta_b = par[2];
    f = -LogLikelihood(mu, theta_s, theta_b); // Minuit минимизирует, поэтому минус
}

// Профилированное правдоподобие (фиксируем mu, минимизируем по theta)
double ProfileLogLikelihood(double mu) {
    TMinuit minuit(3); // 3 параметра: mu (фиксирован), theta_s, theta_b
    minuit.SetFCN(FCN);

    // Устанавливаем параметры
    minuit.DefineParameter(0, "mu", mu, 0.1, 0, 0);
    minuit.FixParameter(0); // Фиксируем mu
    
    minuit.DefineParameter(1, "theta_s", 0, 0.1, -5, 5);
    minuit.DefineParameter(2, "theta_b", 0, 0.1, -5, 5);

    // Настраиваем минимизацию
    double arglist[2];
    arglist[0] = 5000; // Максимальное число вызовов
    arglist[1] = 0.01; // Точность
    int ierr = 0;
    minuit.mnexcm("MIGRAD", arglist, 2, ierr);

    // Получаем статус минимизации (исправленный вызов mnstat)
    double fmin, fedm, errdef;
    int npari, nparx, istat;
    minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);

    if (ierr != 0 || istat < 1) {
        std::cerr << "Warning: Minimization problem for mu=" << mu 
                 << " ierr=" << ierr << " istat=" << istat << std::endl;
    }

    return -fmin; // Возвращаем значение правдоподобия
}

// Тестовая статистика t_mu = -2 (logL(mu) - logL(hat_mu))
double TestStatistic(double mu, double logL_hat) {
    double logL_mu = ProfileLogLikelihood(mu);
    return -2 * (logL_mu - logL_hat);
}

// Основная функция
void Variant4() {
    // 1. Находим глобальный максимум правдоподобия (оценку mu)
    TMinuit minuit(3);
    minuit.SetFCN(FCN);
    minuit.DefineParameter(0, "mu", 1.0, 0.1, 0, 10);
    minuit.DefineParameter(1, "theta_s", 0, 0.1, -5, 5);
    minuit.DefineParameter(2, "theta_b", 0, 0.1, -5, 5);
    minuit.Migrad();

    double mu_hat, mu_err;
    minuit.GetParameter(0, mu_hat, mu_err);
    double logL_hat = LogLikelihood(mu_hat, 0, 0); // Лучшее значение правдоподобия

    std::cout << "Оценка mu_hat = " << mu_hat << " ± " << mu_err << std::endl;

    // 2. Строим графики для разных случаев (все бины, по одному бину)
    std::vector<std::vector<int>> cases = {
        {0, 1, 2},  // Все три бина
        {0},         // Только первый бин
        {1},         // Только второй бин
        {2}          // Только третий бин
    };

    TCanvas *canvas = new TCanvas("canvas", "Variant 4", 800, 600);
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    std::vector<TGraph*> graphs;

    for (size_t icase = 0; icase < cases.size(); ++icase) {
        std::vector<int> bins = cases[icase];
        std::vector<double> mu_values;
        std::vector<double> tmu_values;

        // Сканируем mu от 0 до 2
        for (double mu = 0.0; mu <= 2.0; mu += 0.05) {
            // Временный массив для хранения данных выбранных бинов
            std::vector<double> Ns_case, sigma_s_case, Nb_case, sigma_b_case, Nobs_case;
            for (int bin : bins) {
                Ns_case.push_back(Ns[bin]);
                sigma_s_case.push_back(sigma_s[bin]);
                Nb_case.push_back(Nb[bin]);
                sigma_b_case.push_back(sigma_b[bin]);
                Nobs_case.push_back(Nobs[bin]);
            }

            // Обновляем глобальные переменные для данного случая
            Ns = Ns_case;
            sigma_s = sigma_s_case;
            Nb = Nb_case;
            sigma_b = sigma_b_case;
            Nobs = Nobs_case;

            // Вычисляем тестовую статистику
            double tmu = TestStatistic(mu, logL_hat);
            mu_values.push_back(mu);
            tmu_values.push_back(tmu);
        }

        // Создаём график
        TGraph *graph = new TGraph(mu_values.size(), &mu_values[0], &tmu_values[0]);
        graph->SetLineColor(icase + 1);
        graph->SetLineWidth(2);
        graph->SetTitle(Form("Case \%zu", icase + 1));
        graphs.push_back(graph);

        if (icase == 0) {
            graph->Draw("AL");
            graph->GetXaxis()->SetTitle("#mu");
            graph->GetYaxis()->SetTitle("t_{#mu}");
            graph->GetYaxis()->SetRangeUser(0, 10);
        } else {
            graph->Draw("L SAME");
        }

        legend->AddEntry(graph, Form("Case \%zu", icase + 1), "l");
    }

    legend->Draw();

    // 3. Находим пределы для 68% и 95% CL
    for (size_t icase = 0; icase < cases.size(); ++icase) {
        std::vector<int> bins = cases[icase];
        std::vector<double> Ns_case, sigma_s_case, Nb_case, sigma_b_case, Nobs_case;
        for (int bin : bins) {
            Ns_case.push_back(Ns[bin]);
            sigma_s_case.push_back(sigma_s[bin]);
            Nb_case.push_back(Nb[bin]);
            sigma_b_case.push_back(sigma_b[bin]);
            Nobs_case.push_back(Nobs[bin]);
        }

        Ns = Ns_case;
        sigma_s = sigma_s_case;
        Nb = Nb_case;
        sigma_b = sigma_b_case;
        Nobs = Nobs_case;

        // Находим пределы (перебором)
        double mu_68_low = mu_hat, mu_68_high = mu_hat;
        double mu_95_low = mu_hat, mu_95_high = mu_hat;

        for (double mu = mu_hat; mu >= 0.0; mu -= 0.01) {
            double tmu = TestStatistic(mu, logL_hat);
            if (tmu <= 1.0) mu_68_low = mu;
            if (tmu <= 3.84) mu_95_low = mu;
        }

        for (double mu = mu_hat; mu <= 2.0; mu += 0.01) {
            double tmu = TestStatistic(mu, logL_hat);
            if (tmu <= 1.0) mu_68_high = mu;
            if (tmu <= 3.84) mu_95_high = mu;
        }

        std::cout << "Case " << icase + 1 << ":\n";
        std::cout << "  68\% CL: [" << mu_68_low << ", " << mu_68_high << "]\n";
        std::cout << "  95\% CL: [" << mu_95_low << ", " << mu_95_high << "]\n";
    }

    canvas->SaveAs("Variant4.png");
}