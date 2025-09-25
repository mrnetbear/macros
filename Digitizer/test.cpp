#include <iostream>
#include <fstream>
#include <vector>
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

#pragma pack(push, 1)
struct PsdWaveformPacket{
    char header[4];              // AB 57 46 41
    uint32_t numberOfValues;     // ammount of waveforms in a packet 4 bytes
    int16_t baseline;            // baseline value 2 bytes
    uint16_t decimationFactor;   // decimation factor 2 bytes
    uint64_t timestamp;          // timestamp since power on, ns 8 bytes
    uint16_t channelId;          // channel ID 2 bytes
};
#pragma pack(pop)

#pragma pack(push, 1)
struct DevicePsdSettings{
    uint16_t filterType; // Channel agreement by triggers mode 0 == LED, 1 == CFD 2 bytes
    uint16_t channelId; // channel Id 2 bytes
    uint32_t psdWaveLength; // ammount of waveform values 4 bytes
    uint32_t psdPreTriggerLength; // number of waveform value for Trigger Hold Off gate beginning
    uint32_t triggerHoldOff; // ammount of gate-hitted waveforms (Trigger Hold Off gate)
    uint32_t psdPreGateLength; // Number of waveform value fot Long & Short Gates beginning
    uint32_t psdShortGateLength; // ammount of Short Gate hitted waveform values 
    uint32_t psdLongGateLength; // ammount of Long Gate hitted waveform values
    uint32_t padding; 
};
#pragma pack(pop)

#pragma pack(push, 1)
struct DevicePsdLEDSettings{
    uint32_t ledThresholdUp; // LED activation Threshold
    uint32_t ledThresholdDown; // LED disactivation Threshold
    double padding;
};
#pragma pack(pop)

#pragma pack(push, 1)
struct DevicePsdCFDSettings{
    uint32_t cfdDelay; // CFD signal Delay
    uint32_t cfdThreshold; // CFD Threshold activalion
    double cfdFraction; // Signal weakness for CFD
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

void readWaveformPsd(const char* inputFile, const char* outputFile = "psdWaveform.root"){
    std::ifstream inFile(inputFile, std::ios::binary);
    if (!inFile) {
        std::cerr << "Error: Cannot open input file " << inputFile << std::endl;
        return;
    };

    //Checking for file signature
    char fileSig[8];
    inFile.read(fileSig, 8);
    
    // Expected signature 25 44 47 53 00 01 03 DB
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

    // Читаем количество каналов
    uint16_t numberOfChannels;
    inFile.read(reinterpret_cast<char*>(&numberOfChannels), sizeof(uint16_t));

    // Читаем настройки каналов
    std::vector<DevicePsdSettings> channelSettings(numberOfChannels);
    for (int i = 0; i < numberOfChannels; i++) {
        inFile.read(reinterpret_cast<char*>(&channelSettings[i]), sizeof(DevicePsdSettings));
        
        // Читаем дополнительные настройки в зависимости от filterType
        if (channelSettings[i].filterType == 0) { // LED
            DevicePsdLEDSettings ledSettings;
            inFile.read(reinterpret_cast<char*>(&ledSettings), sizeof(DevicePsdLEDSettings));
            // Здесь можно сохранить LED настройки если нужно
        } else if (channelSettings[i].filterType == 1) { // CFD
            DevicePsdCFDSettings cfdSettings;
            inFile.read(reinterpret_cast<char*>(&cfdSettings), sizeof(DevicePsdCFDSettings));
            // Здесь можно сохранить CFD настройки если нужно
        }
    }

    // Пропускаем MD5 хеш (16 байт)
    inFile.seekg(16, std::ios::cur);

    // Создаем ROOT файл и дерево
    TFile outFile(outputFile, "RECREATE");
    TTree tree("waveformTree", "PSD Waveform Data");

    // Переменные для ветвей дерева
    uint32_t numberOfValues;
    int16_t baseline;
    uint16_t decimationFactor;
    uint64_t timestamp;
    uint16_t channelId;
    std::vector<int16_t> waveformValues;
    
    // Настраиваем ветви дерева
    tree.Branch("numberOfValues", &numberOfValues, "numberOfValues/i");
    tree.Branch("baseline", &baseline, "baseline/S");
    tree.Branch("decimationFactor", &decimationFactor, "decimationFactor/s");
    tree.Branch("timestamp", &timestamp, "timestamp/l");
    tree.Branch("channelId", &channelId, "channelId/s");
    tree.Branch("waveformValues", &waveformValues);

    // Сигнатуры пакетов
    const char expectedHeader[] = {
        static_cast<char>(0xAB), static_cast<char>(0x57),
        static_cast<char>(0x46), static_cast<char>(0x41)
    };
    const char expectedFooter[] = {
        static_cast<char>(0x57), static_cast<char>(0x46),
        static_cast<char>(0x50), static_cast<char>(0xBB)
    };

    // Читаем пакеты waveform данных
    int packetCount = 0;
    PsdWaveformPacket packetHeader;
    
    while (inFile.read(reinterpret_cast<char*>(&packetHeader), sizeof(PsdWaveformPacket))) {
        // Проверяем заголовок пакета
        if (memcmp(packetHeader.header, expectedHeader, 4) != 0) {
            std::cerr << "Warning: Invalid packet header at position " << inFile.tellg() << std::endl;
            continue;
        }

        // Читаем значения waveform
        waveformValues.resize(packetHeader.numberOfValues);
        inFile.read(reinterpret_cast<char*>(waveformValues.data()), 
                   packetHeader.numberOfValues * sizeof(int16_t));

        // Читаем и проверяем footer
        char actualFooter[4];
        inFile.read(actualFooter, 4);
        
        if (memcmp(actualFooter, expectedFooter, 4) != 0) {
            std::cerr << "Warning: Invalid packet footer at position " << inFile.tellg() << std::endl;
            continue;
        }

        // Заполняем переменные для дерева
        numberOfValues = packetHeader.numberOfValues;
        baseline = packetHeader.baseline;
        decimationFactor = packetHeader.decimationFactor;
        timestamp = packetHeader.timestamp;
        channelId = packetHeader.channelId;

        tree.Fill();
        packetCount++;
        
        // Вывод прогресса
        if (packetCount % 100 == 0) {
            std::cout << "Processed " << packetCount << " waveform packets\r" << std::flush;
        }
    }

    std::cout << "\nSuccessfully read " << packetCount << " waveform packets." << std::endl;
    
    tree.Write();
    outFile.Close();
    inFile.close();
}

void test() {
    const char* inputDataFilename = "/usr/local/root/build/macros/Digitizer/data_psd_2025_07_18__16_47_59.bin";
    const char* outputDataFilename = "/usr/local/root/build/macros/Digitizer/psdData_3.root";
    
    readPsdData(inputDataFilename, outputDataFilename);

    const char* inputWaveformFilename = "/usr/local/root/build/macros/Digitizer/waveform_psd_2025_07_18__16_47_59.bin";
    const char* outputWaveformFilename = "/usr/local/root/build/macros/Digitizer/psdWaveform_3.root";
    
    readWaveformPsd(inputWaveformFilename, outputWaveformFilename);

    std::cout << "\nAnalysis commands:" << std::endl;
    std::cout << "  TFile f(\"" << outputDataFilename << "\");" << std::endl;
    std::cout << "  psdTree->Print();" << std::endl;
    std::cout << "  psdTree->Draw(\"psdValue\");" << std::endl;
    std::cout << "  psdTree->Draw(\"qLong:qShort\", \"\", \"colz\");" << std::endl;
    std::cout << "  psdTree->Draw(\"timestamp>>hTS(100,0,1e10)\", \"\", \"\");" << std::endl;

    std::cout << "\nAnalysis commands:" << std::endl;
    std::cout << "  TFile f(\"" << outputWaveformFilename << "\");" << std::endl;
    std::cout << "  psdTree->Print();" << std::endl;
    std::cout << "  psdTree->Draw(\"psdValue\");" << std::endl;
    std::cout << "  psdTree->Draw(\"qLong:qShort\", \"\", \"colz\");" << std::endl;
    std::cout << "  psdTree->Draw(\"timestamp>>hTS(100,0,1e10)\", \"\", \"\");" << std::endl;
}