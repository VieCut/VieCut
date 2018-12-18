//
// Created by noea45cs on 03.10.18.
//

#pragma once

#include <algorithm>
#include <vector>

class vector
{
public:
    template <typename T>
    static bool contains(const std::vector<T>& vec, const T& elem) {
        return (std::find(vec.begin(), vec.end(), elem) != vec.end());
    }
};
