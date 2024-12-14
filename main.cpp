#include <iostream>
#include <csignal>
#include <cstdio>

#include "FluidEngine.h"
#include "inputParser.h"

bool isSave = false;

void signalHandler(int nsig) {
    isSave = true;
}


int main(int argc, char* argv[]) {
    (void) signal(SIGQUIT, signalHandler);

    OptionsParser optsParser(argc, argv);
    auto inputFile = optsParser.getOptVal("--input-file");
    auto saveFileName = optsParser.getOptVal("--save-file");
    int pTypeCode = GetTypeCode(optsParser.getOptVal("--p-type"));
    int vTypeCode = GetTypeCode(optsParser.getOptVal("--v-type"));
    int vFlowTypeCode = GetTypeCode(optsParser.getOptVal("--v-flow-type"));

    std::ifstream input(inputFile);
    if (!input.is_open()) {
        throw std::invalid_argument("can't open file");
    }
    int n, m;
    input >> n >> m;
    input.seekg(0, ios::beg);

    auto engine = ProduceEngine(pTypeCode, vTypeCode, vFlowTypeCode, n, m);

    std::ofstream saveFile(saveFileName);
    if (!input.is_open()) {
        throw std::invalid_argument("can't open file");
    }


    engine->load(input);
    int T = 1'000'000;
    for (int i = 0; i < T; ++i) {
        if (isSave) {
            std::cout << std::endl << "Saving..." << std::endl;
            engine->save(saveFile);
            isSave = false;
            std::cout << "Saved to " + saveFileName << std::endl;
            std::cout << "Press any button to continue" << std::endl;
            getchar();
        }
        engine->next(std::cout);
    }

    input.close();
    saveFile.close();
    return 0;
}