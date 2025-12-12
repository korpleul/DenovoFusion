//
// Created by xinwei on 12/4/24.
//

#include "utils.h"

#include <ctime>
#include <iomanip>
#include <sstream>

std::string get_time_string() {
    time_t now = time(0);
    char buffer[100];
    strftime(buffer, sizeof(buffer), "[%Y-%m-%d %X]", localtime(&now));
    return buffer;
}

std::string get_hhmmss_string(unsigned long long seconds) {
    std::ostringstream oss;
    oss << std::setfill('0');
    oss << std::setw(2) << (seconds / 3600) << ":";
    seconds %= 3600;
    oss << std::setw(2) << (seconds / 60) << ":";
    seconds %= 60;
    oss << std::setw(2) << seconds;
    return oss.str();
}