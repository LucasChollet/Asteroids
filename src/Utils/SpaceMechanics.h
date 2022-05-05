#pragma once

#include "Body.h"

constexpr
double compute_semi_major_axis(Body b) {
    auto const distance = norm(b.position);
    auto const squared_speed = norm2(b.speed);

    return 1. / (2. / distance - squared_speed / mu_s);
}

constexpr
double compute_eccentricity(Body b) {
    auto const tmp1 = glm::cross(b.speed, glm::cross(b.position, b.speed)) / mu_s;
    auto const tmp2 = b.position / norm(b.position);

    return norm(tmp1 - tmp2);
}

constexpr
double compute_inclination(Body b) {
    auto const tmp = glm::cross(b.position, b.speed) / (norm(b.position) * norm(b.speed));
    return gcem::acos(tmp.z);
}

constexpr
Body kep2cart(double semi_major, double eccentricity, double inclination, double omega,
              double Omega, double time, double M0, double initial_time) {

    auto const n = k / gcem::pow(semi_major, 1.5);
    auto const M = M0 + n * (time - initial_time);

    // Newton
    auto const E = [=]() {
        auto tmp = M;
        while (gcem::abs(tmp - eccentricity * gcem::sin(tmp) - M) > 1e-8)
            tmp -= (tmp - eccentricity * gcem::sin(tmp) - M) / (1 - eccentricity * gcem::cos(tmp));
        return tmp;
    }();

    auto const X = semi_major * (gcem::cos(E) - eccentricity);
    auto const Y = semi_major * gcem::sqrt(1 - eccentricity * eccentricity) * gcem::sin(E);

    auto const n_a2_on_r = n * semi_major / (1 - eccentricity * gcem::cos(E));
    auto const X_dot = -n_a2_on_r * gcem::sin(E);
    auto const Y_dot = n_a2_on_r * gcem::sqrt(1 - eccentricity * eccentricity) * gcem::cos(E);

    // Rotation Matrix helper
    auto const create_rot1 = [](double param) -> mat3 {
        auto const tmp1 = gcem::cos(param);
        auto const tmp2 = gcem::sin(param);
        return {1, 0, 0, 0, tmp1, -tmp2, 0, tmp2, tmp1};
    };

    auto const create_rot3 = [](double param) -> mat3 {
        auto const tmp1 = gcem::cos(param);
        auto const tmp2 = gcem::sin(param);
        return {tmp1, -tmp2, 0, tmp2, tmp1, 0, 0, 0, 1};
    };

    // Rotation matrix
    auto const R3W = create_rot3(-Omega);
    auto const R1i = create_rot1(-inclination);
    auto const R3w = create_rot3(-omega);

    auto const transformation_matrix = R3W * R1i * R3w;
    return {transformation_matrix * vec3{X, Y, 0}, transformation_matrix * vec3{X_dot, Y_dot, 0}};
}
