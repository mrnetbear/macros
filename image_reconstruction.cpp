#define Pi 3.1415926535
#define mPi 135.0
#define Epi 6.3e3
#define zpos 312.0

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <time.h>
#include <random>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"

using namespace std;

typedef unsigned int uInt;

std::vector<std::pair<std::pair<double, double>, double>> ReadData(const std::string& filename) {
    std::vector<std::pair<std::pair<double, double>, double>> data;
    std::ifstream file(filename);
    double r, theta;
    while (file >> r >> theta) {
        bool uniqueFlag = true; //флаг на уникальность
        int lorNumber = 0;
        for (uInt i = 0; i < data.size(); i++){
            if (r == data[i].first.first && theta == data[i].first.second){ // поиск уникальности
                uniqueFlag = false;
                lorNumber = i;
                break;
            }
        }
        if (uniqueFlag)
            data.push_back(std::make_pair(std::make_pair(r, theta), 1)); //если уникально, создаём объект с числом появлений 1
        else
            data[lorNumber].second++; //если не уникально, добавляем появление нужному объекту
    }
    return data;
}

// Функция для применения фильтра к данным преобразования Радона
void ApplyFilter(std::vector<std::pair<std::pair<double, double>, double>>& data) {
    //pair<double, double> filteredData;

    //поиск средней частоты
    double avrgFrequency = 0;
    for (const auto& datum : data) {
        avrgFrequency += datum.second;
    }
    cout << avrgFrequency << endl;
    avrgFrequency /= data.size();
    cout << data.size() << "; " << avrgFrequency << endl;
    //for (size_t i = 0; i < data.size(); i++) {
        //data[i].second = data[i].second * data[i].second / avFreq;
        //data[i].second = (data[i].second > avrgFrequency) ?  data[i].second * 10 : data[i].second / 100;
    //}

    // Пример простого фильтра: фильтр Рамахандрана-Лакшминарайана (Ram-Lak)
    const double filterLength = data.size();
    for (size_t i = 0; i < data.size(); i++) {
        double frequency = 2 * M_PI * i / filterLength;
        double filterValue = (i == 0) ? 1.0 : std::abs(frequency);
        data[i].second *= filterValue;
    }
}

// Функция обратного преобразования Радона
void InverseRadonTransform(std::vector<std::pair<std::pair<double, double>, double>>& data, TH2D* outputImage) {
    // Предполагаем, что outputImage уже инициализирован и имеет нужные размеры
    const int nbinsX = outputImage->GetNbinsX();
    const int nbinsY = outputImage->GetNbinsY();
    const double xMin = outputImage->GetXaxis()->GetXmin();
    const double xMax = outputImage->GetXaxis()->GetXmax();
    const double yMin = outputImage->GetYaxis()->GetXmin();
    const double yMax = outputImage->GetYaxis()->GetXmax();

    // Проходим по всем пикселям изображения
    for (int i = 0; i < nbinsX; i++) {
        for (int j = 0; j < nbinsY; j++) {
            double x = xMin + (xMax - xMin) * i / nbinsX;
            double y = yMin + (yMax - yMin) * j / nbinsY;
            // Для каждого пикселя вычисляем значение, основываясь на данных преобразования Радона
            double sum = 0;
            for (const auto& datum : data) {
                double r = datum.first.first;
                double theta = datum.first.second;
                double phi = atan2(y, x);
                // Проверяем, лежит ли точка (x, y) на линии, заданной r и theta
                if (std::abs(r - (x * cos(theta) + y * sin(theta))) < 0.5) {
                    double weight = cos(theta - phi);
                    sum += datum.second * weight;
                    //sum += datum.second; // Увеличиваем сумму, если точка лежит на линии
                }
            }
            if (i % 5 == 0 && j % 5 == 0)
                cout <<"(" << x << "; " << y << ") " << sum << endl;
            outputImage->SetBinContent(i + 1, j + 1, sum);
        }
    }
}

int magic() {
    const std::string inputFilename = "datGAGG.dat"; // Имя файла с входными данными
    const std::string outputFilename = "outGAGG.dat"; // Имя файла для выходных данных

    // Чтение данных из файла
    std::vector<std::pair<std::pair<double,double>, double>> data = ReadData(inputFilename);

    // Создание изображения для обратного преобразования Радона
    TH2D* image = new TH2D("InverseRadon", "Inverse Radon Transform", 100, -50, 50, 100, -50, 50);


    // Применение фильтра к данным преобразования Радона
    //ApplyFilter(data);
    //ApplyFilter(data);

    // Выполнение обратного преобразования Радона
    InverseRadonTransform(data, image);

    // Создание проекций изображения
    TH1D* imageX = image -> ProjectionX();
    imageX->SetTitle("Inverse Radon Transform. X projection");
    TH1D* imageY = image -> ProjectionY();
    imageY->SetTitle("Inverse Radon Transform. Y projection");
    // Сохранение результата в текстовый файл
    std::ofstream outputFile(outputFilename);
    for (int i = 0; i < image->GetNbinsX(); ++i) {
        for (int j = 0; j < image->GetNbinsY(); ++j) {
            outputFile << image->GetXaxis()->GetBinCenter(i) << " "
                       << image->GetYaxis()->GetBinCenter(j) << " "
                       << image->GetBinContent(i, j) << std::endl;
        }
    }

    // Отображение результата в виде изображения
    TCanvas* canvas = new TCanvas("Canvas", "Canvas");
    canvas -> Divide(3,1);
    canvas -> cd(1);
    image->Draw("COLZ");
    canvas -> cd(2);
    imageX->Draw("HIST");
    canvas -> cd(3);
    imageY->Draw("HIST");

    // Сохранение канвы в файл (например, в формате .pdf)
    canvas->SaveAs("output.pdf");

    return 0;
}