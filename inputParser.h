
#ifndef FLUID_PROJECT_INPUTPARSER_H
#define FLUID_PROJECT_INPUTPARSER_H

#include <string>
#include <unordered_map>
#include <algorithm>
#include <cstdio>
#include "engineFactory.h"

enum class EType {
    FIXED,
    F_FIXED,
    FL,
    DB,
};

class OptionsParser{
public:
    OptionsParser(int argc, char* argv[]) {
        for (int i = 1; i < argc; ++i) {
            std::string line = argv[i];
            auto eqPlace = std::find(line.begin(), line.end(), '=');
            if (eqPlace == line.end()) {
                throw std::invalid_argument("Invalid option '" + line + "'");
            }
            auto key = line.substr(0, eqPlace - line.begin());
            auto val = line.substr(eqPlace - line.begin() + 1, line.end() - eqPlace - 1);
            opts_.insert({key, val});
        }
    }

    std::string getOptVal(const std::string& opt) {
        if (!opts_.contains(opt)) {
            throw std::invalid_argument("No such option '" + opt + "'");
        }
        return opts_[opt];
    }
private:
    std::unordered_map<std::string, std::string> opts_;
};


int GetTypeCode(const std::string& typeName) {
    int n, k;
    if (sscanf(typeName.c_str(), "FIXED(%d,%d)", &n, &k) == 2) {
        return FIXED(n, k);
    }
    if (sscanf(typeName.c_str(), "FAST_FIXED(%d,%d)", &n, &k) == 2) {
        return FAST_FIXED(n, k);
    }
    if (typeName == "DOUBLE") {
        return DOUBLE;
    }
    if (typeName == "FLOAT") {
        return FLOAT;
    }
    throw std::invalid_argument("Unknown type '" + typeName + "'");
}

#endif //FLUID_PROJECT_INPUTPARSER_H
