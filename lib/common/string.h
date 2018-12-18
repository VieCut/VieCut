/*******************************************************************************
 * string.h
 *
 * Some string helper functions
 *
 * From Project Thrill - http://project-thrill.org
 *
 * Copyright (C) 2015 Timo Bingmann <tb@panthema.net>
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 ******************************************************************************/

#pragma once

//! Logging helper to print vectors as [a1,a2,a3,...]
template <typename T>
static std::string VecToStr(const std::vector<T>& data) {
    std::ostringstream oss;
    oss << '[';
    for (typename std::vector<T>::const_iterator it = data.begin();
         it != data.end(); ++it)
    {
        if (it != data.begin()) oss << ',';
        oss << *it;
    }
    oss << ']';
    return oss.str();
}

template <typename T>
static typename std::enable_if<std::is_integral<T>::value, std::string>::type
NonZeroVecToStr(const std::vector<T>& data) {
    std::ostringstream oss;
    oss << '[';

    bool first = true;
    for (typename std::vector<T>::const_iterator it = data.begin();
         it != data.end(); ++it)
    {
        size_t i = it - data.begin();
        if (*it != 0) {
            if (first) {
                first = false;
            }
            else {
                oss << ',';
            }
            oss << i << ":";
            oss << *it;
        }
    }
    oss << ']';
    return oss.str();
}

template <typename T>
static typename std::enable_if<!std::is_integral<T>::value, std::string>::type
NonZeroVecToStr(const std::vector<T>& data) {
    return "This method only works for integral types. Use VecToStr()";
}
