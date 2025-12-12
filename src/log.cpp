//
// Created by xinwei on 6/6/24.
//

#include <fstream>
#include "log.h"

std::ofstream Logger::logFile;

void Logger::init(const options_t& options) {
    logFile.open(options.output + "/" + options.prefix + ".log", std::ios::out); // Clear existing content
    logFile.close(); // Close the file
    logFile.open(options.output + "/" +options.prefix + ".log", std::ios::out | std::ios::app); // Reopen in append mode
}

void Logger::close() {
    if (logFile.is_open()) {
        logFile.close();
    }
}

void Logger::Log(const std::string& message, LogLevel level) {
    // Prepare the log level prefix as a string
    std::string logLevelPrefix;
    switch (level) {
        case LogLevel::DEBUG:
            logLevelPrefix = "[DEBUG]: ";
            break;
        case LogLevel::INFO:
            logLevelPrefix = "[INFO]: ";
            break;
        case LogLevel::WARNING:
            logLevelPrefix = "[WARNING]: ";
            break;
        case LogLevel::ERROR:
            logLevelPrefix = "[ERROR]: ";
            break;
    }

    // Write to the log file with level prefix and message
    if (logFile.is_open()) {
        logFile << logLevelPrefix << message << std::endl;
    }
}


void Logger::Debug(const std::string& message) {
    Log(message, LogLevel::DEBUG);
}

void Logger::Info(const std::string& message) {
    Log(message, LogLevel::INFO);
}

void Logger::Warning(const std::string& message) {
    Log(message, LogLevel::WARNING);
}

void Logger::Error(const std::string& message) {
    Log(message, LogLevel::ERROR);
}