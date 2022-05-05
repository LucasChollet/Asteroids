#pragma once

#include "Body.h"
#include "SpaceMechanics.h"

constexpr
vec3 compute_gravity_3body(Body b, double time) {
    // TODO: Enable n-body simulation
    auto const position_jupiter = kep2cart(aJ, eJ, iJ, wJ, WJ, time, M0J, T0J).position;
    auto const &position = b.position;

    auto const distance_to_sun = norm(b.position);
    auto const solar_perturbation = -G * (Ms + mA) / cube(distance_to_sun) * position;

    auto const jupiter_to_pos = position - position_jupiter;
    auto const jupiter_perturbation = -G * mJ * (
            jupiter_to_pos / cube(norm((jupiter_to_pos)))
            + position_jupiter / cube(norm(position_jupiter)));
    return solar_perturbation + jupiter_perturbation;
}

constexpr
Body derivate(Body const &b, double time) {
    return {b.speed, compute_gravity_3body(b, time)};
}
