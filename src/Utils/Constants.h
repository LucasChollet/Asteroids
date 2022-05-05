#pragma once

#include <cmath>
#include <numbers>

#include <gcem.hpp>
#include <glm/vec3.hpp>
#include <glm/gtx/norm.hpp>

// Universal
double constexpr G = 2.95824e-4;
double constexpr k = gcem::sqrt(G);
double constexpr UA = 149.5978707;

double constexpr year = 365.25;

// Sun
double constexpr Ms = 1;
double constexpr mu_s = G * Ms;

// Jupiter
double constexpr mJ = 1 / 1047.348625;
double constexpr aJ = 5.202575;
double constexpr eJ = 0.048908;
double constexpr iJ = 1.3038 * std::numbers::pi / 180;
double constexpr WJ = 100.5145 * std::numbers::pi / 180;
double constexpr wJ = 273.8752 * std::numbers::pi / 180;
double constexpr M0J = 80.0392 * std::numbers::pi / 180;
double constexpr T0J = 2456600.5;

//  Asteroid 2007 VW266
double constexpr mA = 0;
double constexpr aA = 5.454;
double constexpr eA = 0.3896;
double constexpr iA = 108.358 * std::numbers::pi / 180;
double constexpr wA = 226.107 * std::numbers::pi / 180;
double constexpr WA = 276.509 * std::numbers::pi / 180;
double constexpr M0A = 146.88 * std::numbers::pi / 180;
double constexpr T0A = T0J;

template<typename T>
constexpr
T to_degrees(T t) {
    return t * 180 / std::numbers::pi_v<T>;
}

template<typename T>
constexpr
T cube(T t) {
    return t * t * t;
}

// Define GLM constexpr equivalent
using vec3 = glm::vec<3, double>;
using mat3 = glm::mat<3, 3, double>;

template<typename T>
constexpr
typename T::value_type norm2(T const &t) {
    return glm::dot(t, t);
}

template<typename T>
constexpr
typename T::value_type norm(T const &t) {
    return gcem::sqrt(norm2(t));
}

constexpr
mat3 operator*(mat3 const &mat1, mat3 const &mat2) {
    auto const SrcA00 = mat1[0][0];
    auto const SrcA01 = mat1[0][1];
    auto const SrcA02 = mat1[0][2];
    auto const SrcA10 = mat1[1][0];
    auto const SrcA11 = mat1[1][1];
    auto const SrcA12 = mat1[1][2];
    auto const SrcA20 = mat1[2][0];
    auto const SrcA21 = mat1[2][1];
    auto const SrcA22 = mat1[2][2];

    auto const SrcB00 = mat2[0][0];
    auto const SrcB01 = mat2[0][1];
    auto const SrcB02 = mat2[0][2];
    auto const SrcB10 = mat2[1][0];
    auto const SrcB11 = mat2[1][1];
    auto const SrcB12 = mat2[1][2];
    auto const SrcB20 = mat2[2][0];
    auto const SrcB21 = mat2[2][1];
    auto const SrcB22 = mat2[2][2];

    return {SrcA00 * SrcB00 + SrcA01 * SrcB10 + SrcA02 * SrcB20,
            SrcA00 * SrcB01 + SrcA01 * SrcB11 + SrcA02 * SrcB21,
            SrcA00 * SrcB02 + SrcA01 * SrcB12 + SrcA02 * SrcB22,
            SrcA10 * SrcB00 + SrcA11 * SrcB10 + SrcA12 * SrcB20,
            SrcA10 * SrcB01 + SrcA11 * SrcB11 + SrcA12 * SrcB21,
            SrcA10 * SrcB02 + SrcA11 * SrcB12 + SrcA12 * SrcB22,
            SrcA20 * SrcB00 + SrcA21 * SrcB10 + SrcA22 * SrcB20,
            SrcA20 * SrcB01 + SrcA21 * SrcB11 + SrcA22 * SrcB21,
            SrcA20 * SrcB02 + SrcA21 * SrcB12 + SrcA22 * SrcB22,};
}

constexpr
vec3 operator*(mat3 const &m, vec3 const &v) {
    return {
            m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z,
            m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z,
            m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z
    };
}


