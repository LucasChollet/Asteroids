#pragma once

#include "Constants.h"

struct Body {
    vec3 position{};
    vec3 speed{};

    constexpr Body operator+(Body const& other) const {
        return {position + other.position, speed + other.speed};
    }

    constexpr Body operator*(double p) const {
        return {position * p, speed * p};
    }
};
