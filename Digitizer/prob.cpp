#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>




void prob()
{
    // Открываем ROOT-файл с деревом
    TFile* file = TFile::Open("psdData_3.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file psdData.root" << std::endl;
        return;
    }

    // Получаем дерево из файла
    TTree* tree = (TTree*)file->Get("psdTree");
    if (!tree) {
        std::cerr << "Error: Cannot find tree 'psdTree' in file" << std::endl;
        file->Close();
        return;
    }

    // Создаем гистограмму
    TH1F* hQLong = new TH1F("hQLong", "Distribution of qLong;qLong [a.u.];Counts", 10000, 1000, 3000); // Автоматический диапазон

    // Заполняем гистограмму из дерева
    tree->Draw("height>>hQLong", "", "goff"); // "goff" - без отрисовки

    // Настраиваем внешний вид гистограммы
    hQLong->SetLineColor(kBlue);
    hQLong->SetLineWidth(2);
    hQLong->SetFillStyle(3004);
    hQLong->SetFillColor(kBlue);

    // Создаем canvas для отрисовки
    TCanvas* c1 = new TCanvas("c1", "qLong Distribution", 800, 600);
    c1->SetGrid();

    // Рисуем гистограмму
    hQLong->Draw();

    // Добавляем статистику
    //gStyle->SetOptStat(1111);

    // Обновляем canvas
    c1->Update();
}