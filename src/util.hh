#pragma once

//==================================================
// Utility Templates
//==================================================

template<int iters, class func>
struct Unroll;

template<int Value>
constexpr int ceil_log2();

template<class T, bool Left, int32_t exp>
constexpr T static_shift(const T init);

template<class T, int32_t exp>
constexpr T static_left_shift(const T init);
template<class T, int32_t exp>
constexpr T static_right_shift(const T init);

template<class T, bool Left>
T shift(const T init, const int32_t exp);

template<class T>
T left_shift(const T init, const int32_t exp);
template<class T>
T right_shift(const T init, const int32_t exp);

template<class T>
constexpr T countl_zero(const T x) noexcept;

template<class T>
constexpr T countl_one(const T x) noexcept;

template<class T>
constexpr T countr_zero(const T x) noexcept;

template<class T>
constexpr T countr_one(const T x) noexcept;

#include "util.inl"