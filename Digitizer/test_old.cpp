#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>

#pragma pack(push, 1)
struct PsdDataPacket {
    char header[4];          // AB 57 46 41
    uint32_t deviceId;       // 0F 00 00 35
    uint16_t channelId;      // 00 00
    uint64_t timestamp;      // 28 80 58 3F 84 00 00 00
    int16_t cfd_y1;          // FF 01
    int16_t cfd_y2;          // 50 FF
    int16_t height;          // 32 09
    int16_t baseline;        // Следующие 2 байта
    int32_t qLong;           // Следующие 4 байта
    int32_t qShort;          // Следующие 4 байта
    int16_t psdValue;        // Следующие 2 байта
    uint32_t eventCounter;   // Следующие 4 байта
    uint32_t eventCounterPsd;// Следующие 4 байта
    uint16_t decimationFactor; // Следующие 2 байта
    char footer[4];          // 57 46 50 BB
};
#pragma pack(pop)

void readPsdData(const char* inputFile, const char* outputFile = "psdData.root") {
    std::ifstream inFile(inputFile, std::ios::binary);
    if (!inFile) {
        std::cerr << "Error: Cannot open input file " << inputFile << std::endl;
        return;
    }

    // Проверяем модифицированную сигнатуру файла
    char fileSig[8];
    inFile.read(fileSig, 8);
    
    // Ваша реальная сигнатура: 25 44 47 53 00 01 03 DB
    const char expectedFileSig[] = {
        static_cast<char>(0x25), static_cast<char>(0x44), 
        static_cast<char>(0x47), static_cast<char>(0x53),
        static_cast<char>(0x00), static_cast<char>(0x01),
        static_cast<char>(0x03), static_cast<char>(0xDB)
    };

    if (memcmp(fileSig, expectedFileSig, 8) != 0) {
        std::cerr << "Error: Invalid file signature!" << std::endl;
        inFile.close();
        return;
    }

    // Создаем ROOT файл и дерево
    TFile outFile(outputFile, "RECREATE");
    TTree tree("psdTree", "PSD Data");

    PsdDataPacket packet;
    
    // Настраиваем ветви дерева
    tree.Branch("deviceId", &packet.deviceId, "deviceId/i");
    tree.Branch("channelId", &packet.channelId, "channelId/s");
    tree.Branch("timestamp", &packet.timestamp, "timestamp/l");
    tree.Branch("cfd_y1", &packet.cfd_y1, "cfd_y1/S");
    tree.Branch("cfd_y2", &packet.cfd_y2, "cfd_y2/S");
    tree.Branch("height", &packet.height, "height/S");
    tree.Branch("baseline", &packet.baseline, "baseline/S");
    tree.Branch("qLong", &packet.qLong, "qLong/I");
    tree.Branch("qShort", &packet.qShort, "qShort/I");
    tree.Branch("psdValue", &packet.psdValue, "psdValue/S");
    tree.Branch("eventCounter", &packet.eventCounter, "eventCounter/i");
    tree.Branch("eventCounterPsd", &packet.eventCounterPsd, "eventCounterPsd/i");
    tree.Branch("decimationFactor", &packet.decimationFactor, "decimationFactor/s");

    // Сигнатуры пакетов
    const char expectedHeader[] = {
        static_cast<char>(0xAB), static_cast<char>(0x57),
        static_cast<char>(0x46), static_cast<char>(0x41)
    };
    const char expectedFooter[] = {
        static_cast<char>(0x57), static_cast<char>(0x46),
        static_cast<char>(0x50), static_cast<char>(0xBB)
    };

    // Читаем пакеты данных
    int packetCount = 0;
    while (inFile.read(reinterpret_cast<char*>(&packet), sizeof(PsdDataPacket))) {
        if (memcmp(packet.header, expectedHeader, 4) != 0 || 
            memcmp(packet.footer, expectedFooter, 4) != 0) {
            std::cerr << "Warning: Invalid packet at position " << packetCount << std::endl;
            continue;
        }
        
        tree.Fill();
        packetCount++;
        
        // Вывод прогресса
        if (packetCount % 1000 == 0) {
            std::cout << "Processed " << packetCount << " packets\r" << std::flush;
        }
    }

    std::cout << "\nSuccessfully read " << packetCount << " data packets." << std::endl;
    
    tree.Write();
    outFile.Close();
    inFile.close();
}

void test() {
    const char* inputFilename = "/usr/local/Cellar/root/6.32.08/share/root/macros/Digitizer/data_psd_2025_07_18__16_47_59.bin";
    const char* outputFilename = "/usr/local/Cellar/root/6.32.08/share/root/macros/Digitizer/psdData_3.root";
    
    readPsdData(inputFilename, outputFilename);
    
    std::cout << "\nAnalysis commands:" << std::endl;
    std::cout << "  TFile f(\"" << outputFilename << "\");" << std::endl;
    std::cout << "  psdTree->Print();" << std::endl;
    std::cout << "  psdTree->Draw(\"psdValue\");" << std::endl;
    std::cout << "  psdTree->Draw(\"qLong:qShort\", \"\", \"colz\");" << std::endl;
    std::cout << "  psdTree->Draw(\"timestamp>>hTS(100,0,1e10)\", \"\", \"\");" << std::endl;
}