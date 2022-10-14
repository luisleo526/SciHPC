//
// Created by 溫晧良 on 2022/10/14.
//

#ifndef SCIHPC_SPARSEMATRIXEXCEPTIONS_H
#define SCIHPC_SPARSEMATRIXEXCEPTIONS_H


#include <exception>


class Exception : public std::exception {

public:

    explicit Exception(const std::string &message) : exception(), message(message) {}


    virtual ~Exception(void) throw() {}


    inline std::string getMessage(void) const {
        return this->message;
    }


protected:

    std::string message;

};


class InvalidDimensionsException : public Exception {

public:

    InvalidDimensionsException(const std::string &message) : Exception(message) {}

};


class InvalidCoordinatesException : public Exception {

public:

    InvalidCoordinatesException(const std::string &message) : Exception(message) {}

};


#endif //SCIHPC_SPARSEMATRIXEXCEPTIONS_H
