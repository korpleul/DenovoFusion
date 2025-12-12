//
// Created by xinwei on 5/19/24.
//

#ifndef ERROR_H
#define ERROR_H

#include <exception>
#include <string>

// Base custom exception class
class MyError : public std::exception {
private:
    std::string message;

public:
    explicit MyError(const std::string& msg) : message(msg) {}
    virtual ~MyError() noexcept {}

    virtual const char* what() const noexcept override {
        return message.c_str();
    }
};

// Specific parsing error
class ParsingError : public MyError {
public:
    using MyError::MyError;  // Inherit constructor
};

// Specific parsing error
class PslParserError : public MyError {
public:
    using MyError::MyError;  // Inherit constructor
};

// Specific parsing error
class AlignmentError : public MyError {
public:
    using MyError::MyError;  // Inherit constructor
};

// Specific parsing error
class NoAlignmentError : public MyError {
public:
    using MyError::MyError;  // Inherit constructor
};
// Specific parsing error
class CoordPairError : public MyError {
public:
    using MyError::MyError;  // Inherit constructor
};

#endif // ERROR_H
