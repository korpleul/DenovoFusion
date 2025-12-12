//
// Created by xinwei on 6/6/24.
//

#ifndef FUSION_DETECTION_2_LOG_H
#define FUSION_DETECTION_2_LOG_H

#include <iostream>
#include <string>
#include "options.h"
enum class LogLevel {
    INFO,
    DEBUG,
    WARNING,
    ERROR
};

class Logger {
public:
    static void init(const options_t& options) ;
    static void close();

    static void Log(const std::string& message, LogLevel level = LogLevel::INFO);


    static void Debug(const std::string& message);
    static void Info(const std::string& message);
    static void Warning(const std::string& message);
    static void Error(const std::string& message);


    static std::ofstream logFile;
};
#endif //FUSION_DETECTION_2_LOG_H
