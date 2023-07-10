/**
 * @file vec.hpp
 * @brief Defines type for vector maths
 */
#pragma once

namespace smallpt {

struct Vec {
    double x{};
    double y{};
    double z{};

    constexpr Vec() noexcept = default;
    constexpr Vec(double x, double y, double z) noexcept : x{x}, y{y}, z{z} {}
    constexpr Vec(double x, double y) noexcept : x{x}, y{y} {}
    constexpr Vec(double x) noexcept : x{x} {}

    double length() const noexcept;
    Vec normalize() const noexcept;

    constexpr double dot(Vec const& right) const noexcept {
        return x * right.x + y * right.y + z * right.z;
    }

    constexpr Vec cross(Vec const& right) const noexcept {
        return Vec(
            y * right.z - z * right.y,
            z * right.x - x * right.z,
            x * right.y - y * right.x);
    }
};

constexpr Vec operator+(Vec const& left, Vec const& right) noexcept {
    return Vec(left.x + right.x, left.y + right.y, left.z + right.z);
}

constexpr Vec operator-(Vec const& left, Vec const& right) noexcept {
    return Vec(left.x - right.x, left.y - right.y, left.z - right.z);
}

constexpr Vec operator*(Vec const& left, Vec const& right) noexcept {
    return Vec(left.x * right.x, left.y * right.y, left.z * right.z);
}

constexpr Vec operator/(Vec const& left, Vec const& right) noexcept {
    return Vec(left.x / right.x, left.y / right.y, left.z / right.z);
}

constexpr Vec operator+(Vec const& left, const double right) noexcept {
    return Vec(left.x + right, left.y + right, left.z + right);
}

constexpr Vec operator-(Vec const& left, const double right) noexcept {
    return Vec(left.x - right, left.y - right, left.z - right);
}

constexpr Vec operator*(Vec const& left, const double right) noexcept {
    return Vec(left.x * right, left.y * right, left.z * right);
}

constexpr Vec operator/(Vec const& left, const double right) noexcept {
    return Vec(left.x / right, left.y / right, left.z / right);
}

constexpr Vec& operator+=(Vec& left, Vec const& right) noexcept {
    left.x += right.x;
    left.y += right.y;
    left.z += right.z;
    return left;
}

constexpr Vec& operator-=(Vec& left, Vec const& right) noexcept {
    left.x -= right.x;
    left.y -= right.y;
    left.z -= right.z;
    return left;
}

constexpr Vec& operator*=(Vec& left, Vec const& right) noexcept {
    left.x *= right.x;
    left.y *= right.y;
    left.z *= right.z;
    return left;
}

constexpr Vec& operator/=(Vec& left, Vec const& right) noexcept {
    left.x /= right.x;
    left.y /= right.y;
    left.z /= right.z;
    return left;
}

constexpr Vec& operator+=(Vec& left, const double right) noexcept {
    left.x += right;
    left.y += right;
    left.z += right;
    return left;
}

constexpr Vec& operator-=(Vec& left, const double right) noexcept {
    left.x -= right;
    left.y -= right;
    left.z -= right;
    return left;
}

constexpr Vec& operator*=(Vec& left, const double right) noexcept {
    left.x *= right;
    left.y *= right;
    left.z *= right;
    return left;
}

constexpr Vec& operator/=(Vec& left, const double right) noexcept {
    left.x /= right;
    left.y /= right;
    left.z /= right;
    return left;
}

}  // namespace smallpt
