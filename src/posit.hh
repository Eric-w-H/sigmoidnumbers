#pragma once

//==================================================
// Standard library includes
//==================================================

#include <bitset>
#include <bit>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <type_traits>
#include <array>
#include <iostream>

//==================================================
// Internal includes
//==================================================

#include "util.hh"

//==================================================
// Library namespace
//==================================================
namespace posit {


//==================================================
// Forward declarations
//==================================================
template<int32_t, int32_t, int32_t = 30> struct quire;
template<int32_t, int32_t> struct posit;


//==================================================
// Standard aliases
//==================================================

using posit8_t = posit<8, 2>;
using posit16_t = posit<16, 2>;
using posit32_t = posit<32, 2>;
using posit64_t = posit<64, 2>;
using posit128_t = posit<128, 2>;

using quire8_t = quire<8, 2>;
using quire16_t = quire<16, 2>;
using quire32_t = quire<32, 2>;
using quire64_t = quire<64, 2>;
using quire128_t = quire<128, 2>;


//==================================================
// Stubbing
//==================================================

/** 
 * Quire
 */
template<int32_t BITS, int32_t ES, int32_t C> struct quire
{
//--------------------------------------------------
// Friends
//--------------------------------------------------
template<int32_t _BITS, int32_t _ES, int32_t _C>
friend struct quire;
friend struct posit<BITS, ES>;

/**
 * Equal-precision posit type
 */
using posit_t = posit<BITS, ES>;


//--------------------------------------------------
// Internal data
//--------------------------------------------------
private:

// Calculate the size of the whole quire, including carry guard
constexpr static uint32_t ideal_quiresize = std::__bit_ceil( (4*BITS - 8)*(1 << ES) + 1 + C );

// divide by the bitwidth of our type, rounded up
constexpr static uint32_t array_size = (ideal_quiresize + 7) / 8; 
std::array<int8_t, array_size> m_data;

constexpr static uint32_t quiresize = array_size * 8;
constexpr static uint32_t impl_overshoot = quiresize - ideal_quiresize;

public:

//--------------------------------------------------
// Constructors
//--------------------------------------------------
// Quires may only be constructed from posits, or from other quires
constexpr quire(const posit_t& p) noexcept;
quire(const quire& q) = default;
~quire() = default;

// Quires may not be default constructed
quire() = delete;


//--------------------------------------------------
// Convenience functions
//--------------------------------------------------
// Constants
static constexpr quire NaR() noexcept;
static constexpr int32_t precision() noexcept;
static constexpr int32_t exponent() noexcept;
static constexpr int32_t sumLimit() noexcept;

// methods
constexpr bool isNaR() const noexcept;
constexpr posit_t sign() const noexcept;

// friend methods
constexpr friend bool isNaR(const quire& obj) noexcept;
constexpr friend posit_t& sign(const quire& obj) noexcept;


//--------------------------------------------------
// Conversion
//--------------------------------------------------
constexpr explicit operator posit_t() const noexcept;
template<int32_t _BITS, int32_t _ES, int32_t _C>
constexpr explicit operator quire<_BITS, _ES, _C>() const noexcept;
template<int32_t _BITS>
constexpr explicit operator quire<_BITS, ES, C>() const noexcept;
template<int32_t _ES>
constexpr explicit operator quire<BITS, _ES, C>() const noexcept;
template<int32_t _C>
constexpr explicit operator quire<BITS, ES, _C>() const noexcept;

//--------------------------------------------------
// Functions involving quires as member functions
//--------------------------------------------------
constexpr quire qNegate() const noexcept;
constexpr quire qAbs() const noexcept;

constexpr quire qAddQ(const quire& toAdd) const noexcept;
constexpr quire qSubQ(const quire& toSub) const noexcept;

constexpr quire qAddP(const posit_t& toAdd) const noexcept;
constexpr quire qSubP(const posit_t& toSub) const noexcept;

constexpr quire qMulAdd(const posit_t& toMulA, const posit_t& toMulB) const noexcept;
constexpr quire qMulSub(const posit_t& toMulA, const posit_t& toMulB) const noexcept;

constexpr posit_t qToP() const noexcept;


//--------------------------------------------------
// Functions involving quires as friend functions
//--------------------------------------------------
constexpr friend quire qNegate(const quire& obj) noexcept;
constexpr friend quire qAbs(const quire& obj) noexcept;

constexpr friend quire qAddQ(const quire& obj, const quire& toAdd) noexcept;
constexpr friend quire qSubQ(const quire& obj, const quire& toSub) noexcept;

constexpr friend quire qAddP(const quire& obj, const posit_t& toAdd) noexcept;
constexpr friend quire qSubP(const quire& obj, const posit_t& toSub) noexcept;

constexpr friend quire qMulAdd(const quire& obj, const posit_t& toMulA, const posit_t& toMulB) noexcept;
constexpr friend quire qMulSub(const quire& obj, const posit_t& toMulA, const posit_t& toMulB) noexcept;

constexpr friend posit_t qToP(const quire& obj) noexcept;


//--------------------------------------------------
// operators
//--------------------------------------------------
constexpr quire operator=(const quire& other) = default;
constexpr quire operator=(quire&& other) = default;

constexpr friend quire operator+(const quire& lhs, const quire& rhs) noexcept;
constexpr friend quire operator-(const quire& lhs, const quire& rhs) noexcept;

constexpr friend quire operator+(const quire& lhs, const posit_t& rhs) noexcept;
constexpr friend quire operator-(const quire& lhs, const posit_t& rhs) noexcept;

// all other operators are `= delete' but since a quire is non-trivial, don't need to be listed
};


/**
 * Struct to get internal representation type
 */
template<int32_t _BITS, class = void>
struct posit_internal_repr;

template<int32_t _BITS>
struct posit_internal_repr<_BITS, std::enable_if_t<
  (_BITS <= 8)
>> { using type = int8_t; };

template<int32_t _BITS>
struct posit_internal_repr<_BITS, std::enable_if_t<
  (8 < _BITS && _BITS <= 16)
>> { using type = int16_t; };

template<int32_t _BITS>
struct posit_internal_repr<_BITS, std::enable_if_t<
  (16 < _BITS && _BITS <= 32)
>> { using type = int32_t; };

template<int32_t _BITS>
struct posit_internal_repr<_BITS, std::enable_if_t<
  (32 < _BITS && _BITS <= 64)
>> { using type = int64_t; };

template<int32_t _BITS>
struct posit_internal_repr<_BITS, std::enable_if_t<
  (64 < _BITS && _BITS <= 128)
>> { using type = __int128_t; };

template<int32_t _BITS>
using posit_internal_repr_t = typename posit_internal_repr<_BITS>::type;

/** 
 * Posit
 */
template<int32_t BITS, int32_t ES>
struct posit
{
//--------------------------------------------------
// Friends
//--------------------------------------------------
template<int32_t _BITS, int32_t _ES>
friend struct posit;
friend struct quire<BITS, ES>;


/**
 * Internal integer type used to store posit
 */
using posit_data_t = posit_internal_repr_t<BITS>;

/**
 * Equivalent-precision quire type
 */
using quire_t = quire<BITS, ES>;


//--------------------------------------------------
// Member data
//--------------------------------------------------
private:
posit_data_t m_data;

static constexpr posit_data_t power_to_bit_encoding(posit_data_t pow) noexcept;
static constexpr posit data_to_posit(posit_data_t dat) noexcept;
public:


//--------------------------------------------------
// Constructors
//--------------------------------------------------
constexpr posit() noexcept;
constexpr posit(const posit& other) noexcept;
constexpr posit(const std::bitset<BITS>& val) noexcept;

template<typename Integer, std::enable_if_t<std::is_integral<Integer>::value, bool> = true>
constexpr posit(Integer val) noexcept;
template<typename Float, std::enable_if_t<std::is_floating_point<Float>::value, bool> = true>
constexpr posit(Float val) noexcept;

// Destructor
~posit() = default;


//--------------------------------------------------
// Conversion and Casting
//--------------------------------------------------
// posits can convert to each other explicitly
template<int32_t _BITS, int32_t _ES>
constexpr explicit operator posit<_BITS, _ES>() const noexcept;
template<int32_t _BITS>
constexpr explicit operator posit<_BITS, ES>() const noexcept;
template<int32_t _ES>
constexpr explicit operator posit<BITS, _ES>() const noexcept;

// and to integers
template<typename Integer, std::enable_if_t<std::is_integral<Integer>::value && std::is_signed<Integer>::value, bool> = true>
constexpr explicit operator Integer() const noexcept;
// and to floats
template<typename Float, std::enable_if_t<std::is_floating_point<Float>::value, bool> = true>
constexpr explicit operator Float() const noexcept;
// and to bitsets
constexpr explicit operator std::bitset<BITS>() const noexcept;
// and to quires
constexpr explicit operator quire_t() const noexcept;

// Standard quire conversions
constexpr quire_t ptoQ() const noexcept;
friend constexpr quire_t ptoQ(posit& obj) noexcept;


//--------------------------------------------------
// Convenience Functions
//--------------------------------------------------
// Constants
static constexpr posit minPos() noexcept;
static constexpr posit maxPos() noexcept;
static constexpr posit NaR() noexcept;
static constexpr posit_data_t pIntMax() noexcept;
static constexpr posit_data_t pIntMin() noexcept;
static constexpr int32_t precision() noexcept;

// Methods
constexpr bool isNaR() const noexcept;
constexpr posit_data_t raw() const noexcept;
constexpr posit_data_t fixed_significand() const noexcept; // returns the fixed-point repr of the significand
constexpr posit_data_t fixed_fraction() const noexcept; // returns the fixed-point repr of the fraction
constexpr posit_data_t fixed_sign() const noexcept; // returns an integegral 1, 0, -1, NaR like sign
constexpr posit_data_t regime() const noexcept; // (-BITS < regime < BITS)*(2<<ES)
constexpr posit_data_t regime_bitwidth() const noexcept; // 2 <= regime_bitwidth < BITS
constexpr posit_data_t exponent() const noexcept; // 0 <= exponent < (1<<ES)
constexpr posit_data_t power() const noexcept; // -(2<<ES)+(1<<ES) <= power <= (2<<ES)-(1<<ES)
constexpr posit rtz() const noexcept; // Round a posit towards 0.
constexpr posit fraction() const noexcept;    // 0 <= fraction < 1. Returns 0.B for some posit A.B, B>=0.
constexpr posit significand() const noexcept; // -2 <= significand < -1 OR 0 OR 1 <= significand < 2
constexpr posit sign() const noexcept; // -1 if this < 0; 1 if this > 0; 0 if this == 0; NaR if this is NaR
constexpr posit reciprocal() const noexcept; // 1/this

// Methods as friends
constexpr friend bool isNaR(const posit& obj) noexcept;
constexpr friend posit_data_t raw(const posit& obj) noexcept;
constexpr friend posit_data_t fixed_significand(const posit& obj) noexcept; // returns the fixed-point repr of the significand
constexpr friend posit_data_t regime(const posit& obj) noexcept; // (-BITS < regime < BITS)*(2<<ES)
constexpr friend posit_data_t regime_bitwidth(const posit& obj) noexcept; // 2 <= regime_bitwidth < BITS
constexpr friend posit_data_t exponent(const posit& obj) noexcept; // 0 <= exponent < (1<<ES)
constexpr friend posit fraction(const posit& obj) noexcept;    // 0 <= fraction < 1
constexpr friend posit significand(const posit& obj) noexcept; // -2 <= significand < -1 OR 0 OR 1 <= significand < 2
constexpr friend posit sign(const posit& obj) noexcept; // -1 if obj < 0; 1 if obj > 0; 0 if obj == 0; NaR if obj is NaR
constexpr friend posit reciprocal(const posit& obj) noexcept; // 1/obj


//--------------------------------------------------
// Standard Math Functions
//--------------------------------------------------
// Basic functions of one argument as member methods
constexpr posit negate() const noexcept; 
constexpr posit abs() const noexcept; 

constexpr posit nearestInt() const noexcept; 
constexpr posit ceil() const noexcept; 
constexpr posit floor() const noexcept; 

constexpr posit next() const noexcept; 
constexpr posit prior() const noexcept; 

// Elementary functions of one argument as member methods
constexpr posit sqrt() const noexcept;
constexpr posit rSqrt() const noexcept;

constexpr posit exp() const noexcept;
constexpr posit expMinus1() const noexcept;
constexpr posit exp2() const noexcept;
constexpr posit exp2Minus1() const noexcept;
constexpr posit exp10() const noexcept;
constexpr posit exp10Minus1() const noexcept;

constexpr posit log() const noexcept;
constexpr posit logPlus1() const noexcept;
constexpr posit log2() const noexcept;
constexpr posit log2Plus1() const noexcept;
constexpr posit log10() const noexcept;
constexpr posit log10Plus1() const noexcept;

constexpr posit sin() const noexcept;
constexpr posit sinPi() const noexcept;
constexpr posit cos() const noexcept;
constexpr posit cosPi() const noexcept;
constexpr posit tan() const noexcept;
constexpr posit tanPi() const noexcept;

constexpr posit arcSin() const noexcept;
constexpr posit arcSinPi() const noexcept;
constexpr posit arcCos() const noexcept;
constexpr posit arcCosPi() const noexcept;
constexpr posit arcTan() const noexcept;
constexpr posit arcTanPi() const noexcept;

constexpr posit sinH() const noexcept;
constexpr posit cosH() const noexcept;
constexpr posit tanH() const noexcept;

constexpr posit arcSinH() const noexcept;
constexpr posit arcCosH() const noexcept;
constexpr posit arcTanH() const noexcept;

// Functions of one posit and one integer argument, as member methods
template<typename Integer, std::enable_if_t<std::is_integral<Integer>::value, bool> = true>
constexpr posit compound(Integer i) const noexcept;

template<typename Integer, std::enable_if_t<std::is_integral<Integer>::value, bool> = true>
constexpr posit rootN(Integer i) const noexcept;


//--------------------------------------------------
// Free (friend) Standard Math Functions
//--------------------------------------------------
// Basic functions of one argument as free methods
friend constexpr posit negate(const posit& obj) noexcept; 
friend constexpr posit abs(const posit& obj) noexcept; 

friend constexpr posit nearestInt(const posit& obj) noexcept; 
friend constexpr posit ceil(const posit& obj) noexcept; 
friend constexpr posit floor(const posit& obj) noexcept; 

friend constexpr posit next(const posit& obj) noexcept; 
friend constexpr posit prior(const posit& obj) noexcept; 

// Elementary functions of one argument as free methods
friend constexpr posit sqrt(const posit& obj) noexcept;
friend constexpr posit rSqrt(const posit& obj) noexcept;

friend constexpr posit exp(const posit& obj) noexcept;
friend constexpr posit expMinus1(const posit& obj) noexcept;
friend constexpr posit exp2(const posit& obj) noexcept;
friend constexpr posit exp2Minus1(const posit& obj) noexcept;
friend constexpr posit exp10(const posit& obj) noexcept;
friend constexpr posit exp10Minus1(const posit& obj) noexcept;

friend constexpr posit log(const posit& obj) noexcept;
friend constexpr posit logPlus1(const posit& obj) noexcept;
friend constexpr posit log2(const posit& obj) noexcept;
friend constexpr posit log2Plus1(const posit& obj) noexcept;
friend constexpr posit log10(const posit& obj) noexcept;
friend constexpr posit log10Plus1(const posit& obj) noexcept;

friend constexpr posit sin(const posit& obj) noexcept;
friend constexpr posit sinPi(const posit& obj) noexcept;
friend constexpr posit cos(const posit& obj) noexcept;
friend constexpr posit cosPi(const posit& obj) noexcept;
friend constexpr posit tan(const posit& obj) noexcept;
friend constexpr posit tanPi(const posit& obj) noexcept;

friend constexpr posit arcSin(const posit& obj) noexcept;
friend constexpr posit arcSinPi(const posit& obj) noexcept;
friend constexpr posit arcCos(const posit& obj) noexcept;
friend constexpr posit arcCosPi(const posit& obj) noexcept;
friend constexpr posit arcTan(const posit& obj) noexcept;
friend constexpr posit arcTanPi(const posit& obj) noexcept;

friend constexpr posit sinH(const posit& obj) noexcept;
friend constexpr posit cosH(const posit& obj) noexcept;
friend constexpr posit tanH(const posit& obj) noexcept;

friend constexpr posit arcSinH(const posit& obj) noexcept;
friend constexpr posit arcCosH(const posit& obj) noexcept;
friend constexpr posit arcTanH(const posit& obj) noexcept;

// Functions of one posit and one integer argument, as friend methods
template<typename Integer, std::enable_if_t<std::is_integral<Integer>::value, bool>>
friend constexpr posit compound(const posit& obj, Integer i) noexcept;

template<typename Integer, std::enable_if_t<std::is_integral<Integer>::value, bool>>
friend constexpr posit rootN(const posit& obj, Integer i) noexcept;

// Comparison functions of two posit arguments, as friend methods
friend constexpr posit compareEqual(const posit& lhs, const posit& rhs) noexcept;
friend constexpr posit compareNotEqual(const posit& lhs, const posit& rhs) noexcept;
friend constexpr posit compareGreater(const posit& lhs, const posit& rhs) noexcept;
friend constexpr posit compareGreaterEqual(const posit& lhs, const posit& rhs) noexcept;
friend constexpr posit compareLess(const posit& lhs, const posit& rhs) noexcept;
friend constexpr posit compareLessEqual(const posit& lhs, const posit& rhs) noexcept;

// Arithmetic functions of two posit values
friend constexpr posit addition(const posit& lhs, const posit& rhs) noexcept;
friend constexpr posit subtraction(const posit& lhs, const posit& rhs) noexcept;
friend constexpr posit multiplication(const posit& lhs, const posit& rhs) noexcept;
friend constexpr posit division(const posit& lhs, const posit& rhs) noexcept;

// Elementary functions of two posit values
friend constexpr posit hypot(const posit& first, const posit& second) noexcept;
friend constexpr posit pow(const posit& first, const posit& second) noexcept;
friend constexpr posit arcTan2(const posit& first, const posit& second) noexcept;
friend constexpr posit arcTan2Pi(const posit& first, const posit& second) noexcept;

// Functions of three posit values
friend constexpr posit fMM(const posit& first, const posit& second, const posit& third) noexcept;


//--------------------------------------------------
// operators
//--------------------------------------------------
// assignment operators
constexpr posit& operator=(const posit& other) = default;
constexpr posit& operator=(posit&& other) = default;

// negation
constexpr posit operator-() const noexcept;

// basic arithmetic
constexpr posit operator+(const posit& rhs) const noexcept;
constexpr posit operator-(const posit& rhs) const noexcept;
constexpr posit operator*(const posit& rhs) const noexcept;
constexpr posit operator/(const posit& rhs) const noexcept;

// compound assignment
constexpr posit operator+=(const posit& rhs) noexcept;
constexpr posit operator-=(const posit& rhs) noexcept;
constexpr posit operator*=(const posit& rhs) noexcept;
constexpr posit operator/=(const posit& rhs) noexcept;

// comparisons are trivial
constexpr auto operator<=>(const posit& rhs) const noexcept = default;
constexpr bool operator==(const posit& rhs) const noexcept = default;

//--------------------------------------------------
// Undefined operators
//--------------------------------------------------
// increment and decrement are probably confusing. Disable.
posit& operator++() = delete;    //prefix
posit& operator--() = delete;    //prefix
posit& operator++(int) = delete; //postfix
posit& operator--(int) = delete; //postfix

// logical, remainder, and bitwise operators are undefined. Disable.
posit& operator%(const posit&) = delete;
posit& operator&(const posit&) = delete;
posit& operator|(const posit&) = delete;
posit& operator^(const posit&) = delete;

posit& operator&&(const posit&) = delete;
posit& operator||(const posit&) = delete;

posit& operator~() = delete;
posit& operator!() = delete;

posit& operator<<(const posit&) = delete;
posit& operator>>(const posit&) = delete;

// assignment operators of the above
posit& operator%=(const posit&) = delete;
posit& operator&=(const posit&) = delete;
posit& operator|=(const posit&) = delete;
posit& operator^=(const posit&) = delete;

posit& operator<<=(const posit&) = delete;
posit& operator>>=(const posit&) = delete;

};

//===============================================
// Misc Structs
//===============================================

template<int BITS, int ES, class = void>
struct power_of_pIntMax;

// Provides A and B, the smallest and largest, respectively,
// ES for which pIntMax is maximized for a given number of BITS
template<int32_t BITS>
struct posit_max_consec_int_range_t;

template<int32_t BITS, int32_t ES_TEST, int32_t DIR, class = void>
struct count_until_decr;



template<int32_t BITS>
constexpr int32_t max_exponent()
{
  int32_t ceil = 8*sizeof(posit_internal_repr_t<BITS>);
  int32_t adj = std::__bit_floor(BITS-1);
  return ceil-adj;
}

template<int32_t BITS>
using max_posit_t = posit<BITS,max_exponent<BITS>()>;
template<int32_t BITS>
using max_quire_t = quire<BITS,max_exponent<BITS>()>;


} // end namespace posit

// ostreams
template<int32_t BITS, int32_t ES>
std::ostream& operator<<(std::ostream& os, const posit::template posit<BITS,ES>& obj);

#include "posit.inl"