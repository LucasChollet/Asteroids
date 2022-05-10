#pragma once

#include "Constants.h"

/**
 * Struct to represent a body in the system. It has a vector for its
 * position and one for its velocity
 */
struct Body {
    vec3 position{};
    vec3 speed{};

    /**
     * Define the addition of two bodies (used in Runge-Kutta)
     * @param other
     * @return
     */
    constexpr Body operator+(Body const& other) const {
        return {position + other.position, speed + other.speed};
    }

    /**
     * Define the multiplication by a constant (used in Runge-Kutta)
     * @param p
     * @return
     */
    constexpr Body operator*(double p) const {
        return {position * p, speed * p};
    }
};
