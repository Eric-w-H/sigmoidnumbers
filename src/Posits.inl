#include "Posits.hh"


template<uint32_t BITS>
using get_fast_uint_type_t = typename
std::conditional<BITS <= 8, uint_fast8_t, typename
  std::conditional<BITS <= 16, uint_fast16_t, typename
    std::conditional< BITS <= 32, uint_fast32_t, uint_fast64_t >::type 
  >::type 
>::type;

// {
//   if constexpr (BITS <= 8 ) return uint_fast8_t{};
//   if constexpr (BITS <= 16) return uint_fast16_t{};
//   if constexpr (BITS <= 32) return uint_fast32_t{};
//   if constexpr (BITS <= 64) return uint_fast64_t{};
//   // return std::bitset<BITS>;
//}


template<uint32_t BITS, uint32_t ES>
struct posit_traits
{
  using underlying_type = get_fast_uint_type_t<BITS>;

  struct unpacked_posit
  {
    // Unpacked data
    bool is_negative;
    typename std::make_signed<underlying_type>::type regime;
    underlying_type exponent;
    underlying_type fraction;

    // data widths
    underlying_type regime_bits;
    underlying_type exponent_bits;
    underlying_type fraction_bits;

    bool body_is_zero()
    {
      return ( regime | exponent | fraction | regime_bits | exponent_bits | fraction_bits ) == 0; 
    }
  };

  static constexpr underlying_type nbits = BITS;
  static constexpr underlying_type max_exponent_width = ES;
  static constexpr underlying_type underlying_nbits = sizeof(underlying_type) * 8;

  static constexpr underlying_type useed = 1 << (2 * ES);

  // Useful constants
  static constexpr underlying_type infinity = 1 << (BITS - 1);
  static constexpr underlying_type zero = 0;

  // Useful masks
  static constexpr underlying_type sign_bit   = infinity;
  static constexpr underlying_type regime_bit = infinity >> 1;
  static constexpr underlying_type posit_bits = ~static_cast<underlying_type>(static_cast<typename std::make_signed<underlying_type>::type>(-1) << BITS);
  
  // Useful functions
  static constexpr unpacked_posit  get_unpacked (const underlying_type& in_val) {
    if (in_val == 0) return unpacked_posit{false, 0, 0, 0, 0, 0, 0};
    underlying_type val = in_val;
    // p = ((1-3s) + f) x 2^((1-2s)*((2^es)*r+e+s))
    
    // S
    bool is_negative = (val & sign_bit) > 0;
    
    // R
    underlying_type regime_intermediate = (val << 1) & posit_bits; // remove sign bit
    bool regime_bit_set = (regime_intermediate & sign_bit) > 0;    // check whether the regime is 1-leading or 0-leading
    underlying_type regime_bits = regime_intermediate << (underlying_nbits - nbits); // pack the posit to the MSB of the underlying_type to use countl_*
    
    typename std::make_signed<underlying_type>::type regime = 0;
    if (regime_bit_set)
    {
      regime_bits = std::countl_one(regime_bits);
      regime_bits = std::min(regime_bits, static_cast<underlying_type>(nbits - 1) );
      regime = regime_bits - 1;
    } else {
      regime_bits = std::countl_zero(regime_bits);
      regime_bits = std::min(regime_bits, static_cast<underlying_type>(nbits - 1) );
      regime = -regime_bits;
    }


    // Checking for not-a-real
    bool NaR = (!regime_bit_set) && (regime_bits >= nbits - 1) && is_negative;
    if (NaR) return unpacked_posit{true, 0, 0, 0, 0, 0, 0};

    // exponent
    underlying_type nbits_sub_regime_and_sign_bits = nbits - 1 - (regime_bits + ((regime_bits == nbits - 1) ? 0 : 1));
    underlying_type actual_exponent_width = std::min(max_exponent_width, nbits_sub_regime_and_sign_bits);
    underlying_type exponent_bit_mask = ((1 << actual_exponent_width) - 1) << (nbits_sub_regime_and_sign_bits - actual_exponent_width);
    underlying_type exponent = (exponent_bit_mask & val) >> (nbits_sub_regime_and_sign_bits - actual_exponent_width);
    
    underlying_type non_fraction_bits = (((regime_bits == nbits - 1) ? 1 : 2) + regime_bits + actual_exponent_width);
    underlying_type fraction_mask = posit_bits >> non_fraction_bits;
    underlying_type fraction = val & fraction_mask;

    return unpacked_posit{is_negative, regime, exponent, fraction, regime_bits, actual_exponent_width, static_cast<underlying_type>(nbits - non_fraction_bits)};
  }
};


template<uint32_t BITS, uint32_t ES>
class Posit
{
private:
  using traits = posit_traits<BITS, ES>;
  using dtype = typename traits::underlying_type;

  dtype data_bits;

  Posit(dtype data) : data_bits(data) { };

public:
  std::ostream& print_posit(std::ostream& stream, std::string_view sep = " | ")
  {
    auto unpacked = traits::get_unpacked(data_bits);
    stream << (unpacked.is_negative ? '1' : '0');
    stream << sep;
    stream << +unpacked.regime;
    stream << sep;
    stream << +unpacked.exponent;
    stream << sep;
    stream << +unpacked.fraction;

    return stream;
  }

  std::ostream& print_posit_bitwidths(std::ostream& stream, std::string_view sep = " | ")
  {
    auto unpacked = traits::get_unpacked(data_bits);
    stream << 1;
    stream << sep;
    stream << +unpacked.regime_bits;
    stream << sep;
    stream << +unpacked.exponent_bits;
    stream << sep;
    stream << +unpacked.fraction_bits;

    return stream;
  }

  std::ostream& print_fraction_derivation(std::ostream& stream)
  {
    auto pu = traits::get_unpacked(data_bits);
    if (pu.is_negative && pu.body_is_zero())
    {
      stream << "p is Not a Real (NaR)\n";
      return stream;
    }
    if(pu.body_is_zero())
    {
      stream << "p = 0\n";
      return stream;
    }
    long long power = (pu.is_negative ? -1 : 1) * ((1 << ES) * pu.regime + pu.exponent + (pu.is_negative ? 1 : 0));

    long long frac_denom = 1 << pu.fraction_bits;
    long long frac_num = (pu.is_negative ? -2 : 1) * frac_denom + pu.fraction;

    stream  << "p = ((1-3s) + f) x 2^((1-2s)*((2^es)*r+e+s))\n"
            << "  = (" << (pu.is_negative ? "-2" : "1") << " + " << +pu.fraction << "/" << +frac_denom << ") x 2^(" << +power << ')' << '\n'
            << "  = (" << +frac_num << "/" << +frac_denom << ")" << "*2**" << +power << '\n';
    return stream;
  }

  friend std::ostream& operator<<(std::ostream& lhs, Posit<BITS, ES>& posit)
  {
    auto pu = traits::get_unpacked(posit.data_bits);
    if (pu.is_negative && pu.body_is_zero())
    {
      lhs << "NaR";
      return lhs;
    }
    if(pu.body_is_zero())
    {
      lhs << '0';
      return lhs;
    }

    long long power = (pu.is_negative ? -1 : 1) * ((1 << ES) * pu.regime + pu.exponent + (pu.is_negative ? 1 : 0));

    long long frac_denom = 1 << pu.fraction_bits;
    long long frac_num = (pu.is_negative ? -2 : 1) * frac_denom + pu.fraction;

    lhs << frac_num << '/' << frac_denom << "*2**" << +power;
    return lhs;
  }

  // static Posit from_bits(const dtype& bits) { return Posit{bits}; }
  // static Posit from_bits(const std::bitset<BITS>& bits) { return Posit{static_cast<dtype>(bits.to_ullong())}; }
  // template<typename I, std::enable_if_t<std::is_integral<I>::value, bool> = true>
  // static Posit from_int(I num)
  // {
  //   return Posit{};
  // }

  // template<typename F, std::enable_if_t<std::is_floating_point<F>::value, bool> = true>
  // static Posit from_float(F num)
  // {
  //   return Posit{};
  // }

  //--------------------------------------------------
  // basic arithmetic
  //--------------------------------------------------
  // Unary Prefix
  Posit operator-() const { return Posit{static_cast<dtype>(-static_cast<typename std::make_signed_t<dtype>>(data_bits))}; }
  Posit operator+() const { return Posit{static_cast<dtype>(+static_cast<typename std::make_signed_t<dtype>>(data_bits))}; } 

};


template<uint32_t BITS>
struct quire_traits
{
  static constexpr long nbits = 16 * BITS;
  static constexpr long integal_bits = 8*BITS - 16;
  static constexpr long fractional_bits = 8*BITS - 16;
  static constexpr long carry_bits = 31;

  static constexpr long exponent = 16 - 8*BITS;
};


template<uint32_t BITS>
class Quire
{
private:
  using traits = quire_traits<BITS>;
  int carry;
  std::array<int8_t, BITS-2> integal;
  std::array<int8_t, BITS-2> fraction;

public:

};


