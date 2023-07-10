#include <smallpt/vec.hpp>

#include <cmath>

using namespace smallpt;

double Vec::length() const noexcept {
    return std::sqrt(x * x + y * y + z * z);
}

Vec Vec::normalize() const noexcept {
    return *this / length();
}
