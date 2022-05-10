#pragma once

#include "Body.h"

/**
 * Runge Kutta method to solve the differential equation
 * @tparam Derivator
 * @param b Body
 * @param time
 * @param step
 * @param f function used for the derivation (compute gravity on the body)
 * @return body at time t+1
 */
template <typename Derivator>
constexpr
Body runge_kutta4(Body const& b, double time, double step, Derivator f) {
    auto const half_step = step / 2;
    auto const k1 = f(b, time);
    auto const k2 = f(b + k1 * half_step, time + half_step);
    auto const k3 = f(b + k2 * half_step, time + half_step);
    auto const k4 = f(b + k3 * step, time + step);
    return b + (k1 + k2 * 2 + k3 * 2 + k4) * (step / 6);
}
