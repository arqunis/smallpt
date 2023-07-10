/**
 * @file span.hpp
 * @brief Defines a contiguous view to a list of items.
 */
#pragma once

#include <cstddef>
#include <exception>

namespace smallpt {

class IndexError : public std::exception {
public:
    IndexError() noexcept = default;

    char const* what() const noexcept override {
        return "out of bounds";
    }
};

template <typename T>
struct span {
    using value_type = T;
    using iterator = T*;
    using const_iterator = T*;

    constexpr span() noexcept = default;
    constexpr span(T* data, std::size_t size) noexcept
        : data_{data}, size_{size} {}

    template <std::size_t N>
    constexpr span(T (&data)[N]) noexcept : data_{data}, size_{N} {}

    constexpr T* data() const noexcept {
        return data_;
    }

    constexpr std::size_t size() const noexcept {
        return size_;
    }

    constexpr T& get_unchecked(std::size_t index) const noexcept {
        return data_[index];
    }

    constexpr T& get(std::size_t index) const {
        if (index >= size_) {
            throw IndexError();
        }

        return get_unchecked(index);
    }

    constexpr T& operator[](std::size_t index) const {
        return get(index);
    }

    constexpr T* begin() const noexcept {
        return data_;
    }

    constexpr T* end() const noexcept {
        return data_ + size_;
    }

    constexpr T const* cbegin() const noexcept {
        return data_;
    }

    constexpr T const* cend() const noexcept {
        return data_ + size_;
    }

private:
    T* data_{nullptr};
    std::size_t size_{0};
};

}  // namespace smallpt
