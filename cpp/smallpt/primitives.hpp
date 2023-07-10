/**
 * @file sphere.hpp
 */
#pragma once

#include <smallpt/span.hpp>
#include <smallpt/vec.hpp>

#include <cstddef>

namespace smallpt {

/// @brief Material types.
enum class ReflectionType { Diffuse, Specular, Reflective };

struct Sphere {
    double radius{0.0};

    Vec position;
    Vec emission;
    Vec colour;

    ReflectionType reflection;

    constexpr Sphere(
        double radius,
        Vec position,
        Vec emission,
        Vec colour,
        ReflectionType reflection) noexcept
        : radius{radius}, position{position}, emission{emission},
          colour{colour}, reflection{reflection} {}
};

struct Ray {
    Vec position;
    Vec direction;

    constexpr Ray(Vec position, Vec direction) noexcept
        : position{position}, direction{direction} {}
};

/// @brief Calculates the distance to the intersection.
/// @returns The distance to the ray, 0 if no hit
double intersection(Sphere const& sphere, Ray const& ray) noexcept;

bool intersection(
    span<Sphere const> spheres,
    Ray const& ray,
    double& distance,
    std::size_t& id) noexcept;

}  // namespace smallpt
