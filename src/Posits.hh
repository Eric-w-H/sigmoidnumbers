#pragma once

// fixed width integer types
#include <cstdint>
#include <type_traits>
#include <bit>
#include <bitset>
#include <algorithm>
#include <iostream>

template<uint32_t BITS, uint32_t ES>
struct posit_traits;

template<uint32_t BITS, uint32_t ES>
class Posit;

template<uint32_t BITS>
struct quire_traits;

template<uint32_t BITS>
class Quire;

#include <Posits.inl>