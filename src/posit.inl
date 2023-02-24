#include "posit.hh"

//==================================================
// library namespace
//==================================================
namespace posit {

//==================================================
// implementations
//==================================================

/** 
 * Construct a quire from a posit.
 */
template<int32_t BITS, int32_t ES, int32_t C>
constexpr quire<BITS,ES,C>::quire(const posit_t& p) noexcept
  : m_data{0}
{
  using dtype = typename posit_t::posit_data_t;

  // handle exception values
  if constexpr (p == 0) return;
  if constexpr (p.isNaR()) { m_data[array_size - 1] = 0x80; return; }

  constexpr dtype fraction_bitwidth = std::max(BITS - 1 - p.regime_bitwidth() - ES, 0);
  constexpr dtype fraction_mask = (1 << fraction_bitwidth) - 1;
  constexpr dtype shift_into_quire = p.power() - quire::exponent();
  
  dtype index_into_quire = shift_into_quire / 8;
  dtype fraction_bits = p.m_data & fraction_mask;

  for(; fraction_bits > 0; fraction_bits >>= 8)
  {
    m_data[index_into_quire] = fraction_bits & 0xff;
    index_into_quire += 1;
  }
}


//--------------------------------------------------
// Quire Convenience Functions
//--------------------------------------------------

/**
 * NaR
 * returns the Quire{NaR}
 */
template<int32_t BITS, int32_t ES, int32_t C>
constexpr quire<BITS, ES, C> quire<BITS,ES,C>::NaR() noexcept
{ return quire<BITS,ES>{posit<BITS,ES>::NaR()}; }


/**
 * precision
 * returns the bitwidth of the quire
 */
template<int32_t BITS, int32_t ES, int32_t C>
constexpr int32_t quire<BITS,ES,C>::precision() noexcept { return quiresize; }


/**
 * exponent
 * returns the power of two of the LSB of the quire
 */
template<int32_t BITS, int32_t ES, int32_t C>
constexpr int32_t quire<BITS,ES,C>::exponent() noexcept { return posit_t::minPos().power()*2; }


/**
 * sumLimit
 * returns the power of two of the largest representable number
 */
template<int32_t BITS, int32_t ES, int32_t C>
constexpr int32_t quire<BITS, ES, C>::sumLimit() noexcept { return posit_t::maxPos().power()*2 + 1 + C; }


/**
 * isNaR
 * returns whether the quire is NaR
 */
template<int32_t BITS, int32_t ES, int32_t C>
constexpr bool quire<BITS, ES, C>::isNaR() const noexcept {
  bool isNaR = false;
  isNaR = (m_data[array_size - 1] == 0x80); // 8'b1000_0000
  for(int i = 0; i < array_size - 1; ++i)
  {
    int8_t val = m_data[i];
    isNaR |= (val == 0);
  }
  return isNaR;
}


/**
 * sign
 * Returns a posit that is:
 *  1 if this > 0; 
 *  0 if this == 0;
 * -1 if this < 0;
 *  and NaR if this is NaR
 */
template<int32_t BITS, int32_t ES, int32_t C>
constexpr posit<BITS,ES> quire<BITS,ES,C>::sign() const noexcept
{
  uint8_t msB = bit_cast<uint8_t>(m_data[array_size - 1]);
  bool isNaR      = (msB == 0x80);
  bool isZero     = (msB == 0);
  bool isNegative = (msB < 0);

  // finish out isNaR or isZero
  bool isBodyZero = isNaR || isZero;
  for(int i = array_size - 2; i >= 0 && isBodyZero; --i)
  {
    isBodyZero &= (m_data[i] == 0);
  }

  if(isNaR && isBodyZero) return posit_t::NaR();
  if(isZero && isBodyZero) return posit_t{0};
  if(isNegative) return posit_t{-1};
  return posit_t{1};
}


/**
 * isNaR
 * returns whether the quire is NaR
 */
template<int32_t BITS, int32_t ES, int32_t C>
constexpr bool isNaR(const quire<BITS,ES,C>& obj) noexcept { return obj.isNaR(); }


/**
 * sign
 * Returns a posit that is 1 if this > 0; 0 if this == 0; -1 if this < 0; and NaR if this is NaR
 */
template<int32_t BITS, int32_t ES, int32_t C>
constexpr posit<BITS,ES>& sign(const quire<BITS,ES,C>& obj) noexcept { return obj.sign(); }


/**
 * static_cast<posit_t>
 * Rounds the quire into the corresponding posit type.
 */
template<int32_t BITS, int32_t ES, int32_t C>
constexpr quire<BITS,ES,C>::operator posit<BITS,ES>() const noexcept
{
  posit_t sign = sign();
  if(sign.isNaR() || sign == posit_t{0}) return sign;

  // we are now guaranteed that the number is nonzero
  return sign;
}


//=========================================================
// Explicit casting functions
//=========================================================

template<int BITS, int ES, int C>
template<int C2>
constexpr quire<BITS, ES, C>::operator quire<BITS, ES, C2>() const noexcept
{
  using type_from = quire<BITS, ES, C >;
  using type_to   = quire<BITS, ES, C2>;
  
  // exception handling
  constexpr posit_t sign = this->sign();
  if constexpr (sign.isNaR()) return type_to{typename type_to::posit_t::NaR()};
  if constexpr (sign == 0) return type_to{typename type_to::posit_t{0}};

  constexpr uint32_t quiresize_from = type_from::precision();
  constexpr uint32_t quiresize_to = type_to::precision();

  if constexpr (quiresize_from > quiresize_to)
  {
    // check whether we're truncating data
    bool wasDataTruncated = false;
    constexpr int8_t msB_extend = (sign < 0) ? 0xFF : 0x00; 
    for(int i = type_to::array_size; i < type_from::array_size; ++i) wasDataTruncated |= (msB_extend == m_data[i]);

    // round to maximum (or minimum) quire
    if /* constexpr */ (wasDataTruncated)
    {
      type_to result{posit_t{0}};
      result.m_data[result.array_size - 1] = (sign < 0) ? 0x80 : 0x7F;
      for(int i = 0; i < result.array_size - 1; ++i) result.m_data[i] = (sign < 0) ? 0x00 : 0xFF;
      if constexpr (sign < 0) result.m_data[0] += 1;
      return result;
    }

    // nothing is truncated, just copy
    type_to result{posit_t{0}};
    for(int i = 0; i < result.array_size; ++i) result.m_data[i] = m_data[i];
    return result;
  }

  // need to add carry guard bits
  // constexpr uint32_t expand = quiresize_to - quiresize_from;
  
  type_to result{posit_t{0}};

  // copy
  for(int i = 0; i < array_size; ++i) result.m_data[i] = m_data[i];

  // sign extend
  int8_t sext = (sign < 0) ? 0xFF : 0x00;
  for(int i = array_size; i < result.array_size; ++i) result.m_data[i] = sext;
  return result;
}

template<int BITS, int ES, int C>
template<int BITS2>
constexpr quire<BITS, ES, C>::operator quire<BITS2, ES, C>() const noexcept
{
  constexpr quire<ES   ,BITS ,C> reinterpret = reinterpret_cast<quire<ES,BITS,C>>(*this);
  constexpr quire<ES   ,BITS2,C> cast        = static_cast<quire<ES,BITS2,C>>(reinterpret);
  constexpr quire<BITS2,ES   ,C> result      = reinterpret_cast<quire<BITS2,ES,C>>(cast);
  return result;
}

template<int BITS, int ES, int C>
template<int ES2>
constexpr quire<BITS, ES, C>::operator quire<BITS, ES2, C>() const noexcept
{
  using type_from = quire<BITS, ES , C>;
  using type_to   = quire<BITS, ES2, C>;

  // exception handling
  constexpr posit_t sign = sign();
  if constexpr (sign.isNaR()) return type_to{typename type_to::posit_t::NaR()};
  if constexpr (sign == 0) return type_to{typename type_to::posit_t{0}};

  // maxpos = useed ^ (nbits-2)
  // useed = 2^2^ES
  // maxpos = 2^2^ES^(nbits-2)
  // maxpos = 2^(2*ES*(nbits-2))
  // maxpos_exp = 2*ES*(nbits-2)
  constexpr int32_t exp_stepsize = 2*(BITS - 2);
  
  if constexpr (ES > ES2) {
    // scaling down
    constexpr int32_t power_diff = exp_stepsize * (ES - ES2);
    constexpr int32_t array_diff_bits = type_from::quiresize - type_to::quiresize;
    constexpr int32_t off_the_bottom_bits = 1 << power_diff;
    // total difference minus the bit from the bottom
    constexpr int32_t off_the_top_bits = array_diff_bits - (1 << power_diff);

    constexpr int32_t bottom_start_idx = off_the_bottom_bits/8;
    constexpr int32_t top_start_idx = type_to::array_size + off_the_bottom_bits/8;

    // check whether we're truncating data--top first
    bool wasDataTruncated_top = false;
    constexpr int8_t msB_extend = (sign < 0) ? 0xFF : 0x00;
    for(int i = 0; i < off_the_top_bits/8; ++i) 
      wasDataTruncated_top |= (msB_extend == m_data[i + top_start_idx]);

    // check whether there's data in the body
    bool fromIsSmallerThanToCanRepresent = true;
    for(int i = bottom_start_idx; i < top_start_idx; ++i)
      fromIsSmallerThanToCanRepresent &= (msB_extend == m_data[i]);

    // iff truncating data on the top, set to largest abs value * sign
    // iff not truncating and fromIsSmaller, set to smallest abs value * sign
    // else just copy over data body.
    type_to result{typename type_to::posit_t{0}};

    if (wasDataTruncated_top)
    {
      type_to result{typename type_to::posit_t{0}};
      result.m_data[result.array_size - 1] = (sign < 0) ? 0x80 : 0x7F;
      for(int i = 0; i < result.array_size - 1; ++i) result.m_data[i] = (sign < 0) ? 0x00 : 0xFF;
      if constexpr (sign < 0) result.m_data[0] += 1;
      return result;
    } else if (fromIsSmallerThanToCanRepresent) {
      if (sign < 0)
        for(auto& elem : result.m_data) elem = 0xFF;
      else
        result.m_data[0] = 0x01;
      return result;
    }
    for(int i = 0; i < type_to::array_size; ++i)
      result.m_data[i] = m_data[i + bottom_start_idx];
    return result;
  }
  // scaling up
  constexpr int32_t power_diff = exp_stepsize * (ES2 - ES);
  // constexpr int32_t array_diff_bits = type_to::quiresize - type_from::quiresize;
  constexpr int32_t off_the_bottom_bits = 1 << power_diff;
  // total difference minus the bit from the bottom
  // constexpr int32_t off_the_top_bits = array_diff_bits - (1 << power_diff);

  constexpr int32_t bottom_start_idx = off_the_bottom_bits/8;
  constexpr int32_t top_start_idx = type_to::array_size + off_the_bottom_bits/8;
  
  // copy over data
  type_to result{typename type_to::posit_t{0}};
  for(int i = 0; i < type_from::array_size; ++i)
    result.m_data[i + bottom_start_idx] = m_data[i];

  // sign extend negative
  if(sign < 0)
    for(int i = top_start_idx; i < type_to::array_size; ++i)
      result.m_data[i] = 0xFF;

  return result;
}

template<int BITS, int ES, int C>
template<int BITS2, int ES2, int C2>
constexpr quire<BITS, ES, C>::operator quire<BITS2, ES2, C2>() const noexcept
{
  quire<BITS, ES, C2> trim = static_cast<quire<BITS, ES, C2>>(*this);
  quire<BITS2, ES, C2> scale = static_cast<quire<BITS2, ES, C2>>(trim);
  quire<BITS2, ES2, C2> shift = static_cast<quire<BITS2, ES2, C2>>(scale);
  return shift;
}

//=========================================================
// Other Functions
//=========================================================
template<int BITS, int ES, int C>
constexpr quire<BITS,ES,C> quire<BITS, ES, C>::qNegate() const noexcept
{
  quire<BITS,ES,C> result{posit_t{0}};

  bool found_first_set = false;
  for(int i = 0; i < m_data.size(); ++i)
  {
    if(!found_first_set)
    {
      found_first_set = (m_data[i] != 0);
      if(found_first_set) result.m_data[i] = -m_data[i];
      continue;
    } else {
      result.m_data[i] = ~m_data[i];
    }
  }
  
  return result;
}

template<int BITS, int ES, int C>
constexpr quire<BITS, ES, C> quire<BITS, ES, C>::qAbs() const noexcept
{
  return sign() < 0 ? qNegate() : *this;
}

template<int BITS, int ES, int C>
constexpr quire<BITS, ES, C> quire<BITS, ES, C>::qAddQ(const quire<BITS, ES, C>& toAdd) const noexcept
{
  quire<BITS, ES, C> result{posit_t{0}};
  int16_t carry = 0;
  
  for(int i = 0; i < m_data.size(); ++i)
  {
    carry = m_data[i] + toAdd.m_data[i] + carry;
    result.m_data[i] = (carry & 0xFF);
    carry = (carry >> 8);
  }

  return result;
}

template<int BITS, int ES, int C>
constexpr quire<BITS, ES, C> quire<BITS, ES, C>::qSubQ(const quire<BITS, ES, C>& toSub) const noexcept
{
  return qAddQ(qNegate(toSub));
}

template<int BITS, int ES, int C>
constexpr quire<BITS, ES, C> quire<BITS, ES, C>::qAddP(const posit_t& toAdd) const noexcept
{
  return qAddQ(quire{toAdd});
}

template<int BITS, int ES, int C>
constexpr quire<BITS, ES, C> quire<BITS, ES, C>::qSubP(const posit_t& toSub) const noexcept
{
  return qAddQ(quire{-toSub});
}

template<int BITS, int ES, int C>
constexpr quire<BITS, ES, C> quire<BITS, ES, C>::qMulAdd(const posit_t& toMulA, const posit_t& toMulB) const noexcept
{
  using dtype = typename posit_t::posit_data_t;

  constexpr posit_t signA = toMulA.sign();
  constexpr posit_t signB = toMulB.sign();

  // exception cases
  if(signA.isNaR() || signA == 0) return quire<BITS,ES,C>{signA};
  if(signB.isNaR() || signB == 0) return quire<BITS,ES,C>{signB};
  if(isNaR()) return *this;

  constexpr bool negate_result = (signA < 0) ^ (signB < 0);

  constexpr int32_t FRAC_BITWIDTH = BITS - 3 - ES;
  constexpr dtype powerA = toMulA.power();
  constexpr dtype powerB = toMulB.power();
  constexpr dtype fracA = toMulA.fraction().m_data & (static_shift<dtype,true,FRAC_BITWIDTH>(1) - 1);
  constexpr dtype fracB = toMulB.fraction().m_data & (static_shift<dtype,true,FRAC_BITWIDTH>(1) - 1);

  dtype finalPower = powerA + powerB;

  constexpr int32_t quire_array_bit_offset = finalPower - exponent();
  constexpr int32_t quire_array_index = quire_array_bit_offset / 8;
  constexpr int32_t byte_offset = quire_array_bit_offset - quire_array_index * 8;
  constexpr int32_t byte_span = sizeof(dtype);
  constexpr int32_t nonzero_byte_offset = (byte_offset == 0) ? 0 : 1;

  // if multiplication doesn't fully fit into the internal type
  if constexpr (2 * FRAC_BITWIDTH > sizeof(dtype) * 8) {
    // anonymous functions to help with long multiplication
    auto hi = [](dtype x) { return       static_shift<dtype,false,FRAC_BITWIDTH/2>(x);       };
    auto lo = [](dtype x) { return (x & (static_shift<dtype,true ,FRAC_BITWIDTH/2>(1) - 1)); };

    constexpr dtype x = lo(fracA) * lo(fracB);
    constexpr dtype s0 = lo(x);

    constexpr dtype x2 = hi(fracA) * lo(fracB) + hi(x);
    constexpr dtype s1 = lo(x2);
    constexpr dtype s2 = hi(x2);

    constexpr dtype x3 = s1 + lo(fracA) * hi(fracB);
    constexpr dtype s1_2 = lo(x3);

    constexpr dtype x4 = s2 + hi(fracA) * hi(fracB) + hi(x3);
    constexpr dtype s2_2 = lo(x4);
    constexpr dtype s3 = hi(x4);

    // result is 2xFRAC_BITWIDTH long;
    dtype result = static_shift<dtype,true,FRAC_BITWIDTH/2>(s1_2) | s0;
    dtype carry  = static_shift<dtype,true,FRAC_BITWIDTH/2>(s3  ) | s2_2;

    if constexpr (byte_offset != 0)
      result.m_data[quire_array_index] = static_shift<dtype,true,byte_offset>(result) & 0xFF;
    for(int byte = 0; byte < byte_span; ++byte)
      result.m_data[quire_array_index + byte + nonzero_byte_offset] = static_shift<dtype,false,(8 * byte + 8 - byte_offset)>(result) & 0xFF;

    if constexpr (byte_offset != 0)
      result.m_data[quire_array_index + byte_span] |= static_shift<dtype,true,byte_offset>(carry) & 0xFF;
    for(int byte = 0; byte < byte_span; ++byte)
      result.m_data[quire_array_index + byte_span + byte + nonzero_byte_offset] = static_shift<dtype,false,(8 * byte + 8 - byte_offset)>(carry) & 0xFF;

    if constexpr (negate_result) return qSubQ(result);
    return qAddQ(result);
  }

  dtype result = fracA * fracB;
  if constexpr(byte_offset != 0)
    result.m_data[quire_array_index] = static_shift<dtype,true,byte_offset>(result) & 0xFF;
  for(int byte = 0; byte < byte_span; ++byte)
    result.m_data[quire_array_index + byte + nonzero_byte_offset] = (static_shift<dtype,false,(8 * byte + 8 - byte_offset)>(result) & 0xFF);

  if constexpr (negate_result) return qSubQ(result);
  return qAdd(result);
}

template<int BITS, int ES, int C>
constexpr quire<BITS,ES,C> quire<BITS,ES,C>::qMulSub(const posit_t& toMulA, const posit_t& toMulB) const noexcept
{ return qMulAdd(toMulA, -toMulB); }

template<int BITS, int ES, int C>
constexpr typename quire<BITS,ES,C>::posit_t quire<BITS,ES,C>::qToP() const noexcept
{ return static_cast<posit_t>(*this); }

//-----------------------------------------------
// Friends
//-----------------------------------------------

template<int BITS, int ES, int C>
constexpr quire<BITS,ES,C> qNegate(const quire<BITS,ES,C>& obj) noexcept
{ return obj.qNegate(); }

template<int BITS, int ES, int C>
constexpr quire<BITS,ES,C> qAbs(const quire<BITS,ES,C>& obj) noexcept
{ return obj.qAbs(); }

template<int BITS, int ES, int C>
constexpr quire<BITS,ES,C> qAddQ(const quire<BITS,ES,C>& obj, const quire<BITS,ES,C>& toAdd) noexcept
{ return obj.qAddQ(toAdd); }

template<int BITS, int ES, int C>
constexpr quire<BITS,ES,C> qSubQ(const quire<BITS,ES,C>& obj, const quire<BITS,ES,C>& toSub) noexcept
{ return obj.qSubQ(toSub); }

template<int BITS, int ES, int C>
constexpr quire<BITS,ES,C> qMulAdd(const quire<BITS,ES,C>& obj, const posit<BITS,ES>& toMulA, const posit<BITS,ES>& toMulB) noexcept
{ return obj.qMulAdd(toMulA, toMulB); }

template<int BITS, int ES, int C>
constexpr quire<BITS,ES,C> qMulSub(const quire<BITS,ES,C>& obj, const posit<BITS,ES>& toMulA, const posit<BITS,ES>& toMulB) noexcept
{ return obj.qMulSub(toMulA, toMulB); }

template<int BITS, int ES, int C>
constexpr posit<BITS,ES> qToP(quire<BITS,ES,C>& obj) noexcept
{ return static_cast<posit<BITS,ES>>(obj); }

//-----------------------------------------------
// Operators
//-----------------------------------------------
template<int32_t BITS, int32_t ES, int32_t C>
constexpr quire<BITS,ES,C> operator+(const quire<BITS,ES,C>& lhs, const quire<BITS,ES,C>& rhs) noexcept
{ return qAddQ(lhs, rhs); }

template<int32_t BITS, int32_t ES, int32_t C>
constexpr quire<BITS,ES,C> operator-(const quire<BITS,ES,C>& lhs, const quire<BITS,ES,C>& rhs) noexcept
{ return qSubQ(lhs, rhs); }

template<int32_t BITS, int32_t ES, int32_t C>
constexpr quire<BITS,ES,C> operator+(const quire<BITS,ES,C>& lhs, const posit<BITS,ES>& rhs) noexcept
{ return qAddP(lhs,rhs); }

template<int32_t BITS, int32_t ES, int32_t C>
constexpr quire<BITS,ES,C> operator-(const quire<BITS,ES,C>& lhs, const posit<BITS,ES>& rhs) noexcept
{ return qSubP(lhs,rhs); }

/////////////////////////////////////////////////
// Posits
/////////////////////////////////////////////////

template<int BITS, int ES>
constexpr typename posit<BITS,ES>::posit_data_t 
power_to_bit_encoding(typename posit<BITS,ES>::posit_data_t pow) noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  constexpr posit_data_t regime_coef = static_left_shift<posit_data_t, ES>(1);
  
  // p = (2^ES)*r + e
  // r is signed, e is unsigned
  // (p - e) / (2^ES) = r
  // 2^ES > e, therefore
  // p / (2^ES) = r

  // r = -k if k < 0 else r = k - 1
  // therefore
  // k = -r if r < 0 else k = r + 1
  // and 
  // k0 = 0 if r < 0 else k0 = 1

  // since 0 <= e < 2^ES, and r * (2^ES) <= p
  // 

  constexpr posit_data_t runlen_test = pow / regime_coef;
  constexpr posit_data_t runlen_signed = pow < 0 ? runlen_test - regime_coef : runlen_test;
  constexpr posit_data_t exp = pow - regime_coef * runlen_signed;
  constexpr posit_data_t runlen = pow < 0 ? -runlen_signed : runlen_signed;

  constexpr posit_data_t run_left_shift = BITS - 2 - runlen;
  
  constexpr bool exp_shifts_left = run_left_shift > ES;

  constexpr posit_data_t run_right_encoding = pow < 0 ? 1 : left_shift<posit_data_t>(1, runlen + 1) - 1;
  constexpr posit_data_t run_encoding = left_shift<posit_data_t>(run_right_encoding, run_left_shift);
  constexpr posit_data_t exp_shifted = exp_shifts_left ? left_shift<posit_data_t>(exp, run_left_shift - ES)
                                                      : right_shift<posit_data_t>(exp, ES - run_left_shift);

  // implement posit rounding for the case of a right shift on ES
  constexpr posit_data_t posit_of_power = run_encoding | exp_shifted;
  // rules for promotion: if a bit was shifted out of ES and posit_of_power isn't max_pos
  // AND either:
  //  (shifted out bits of exp = 1) and lsb of posit_of_pwer is 0
  //  OR
  //  shifted out bits of exp > 1
  constexpr posit_data_t shifted_out_bits_of_exp = exp_shifts_left ? 0 : exp & (left_shift<posit_data_t>(2, ES - run_left_shift) - 1);
  constexpr bool promote_to_next_posit = (posit_of_power == posit<BITS,ES>::maxPos().m_data) ? false
    : (shifted_out_bits_of_exp > 1) || ((shifted_out_bits_of_exp == 1) && ((posit_of_power & 1) == 1));

  return posit_of_power + (promote_to_next_posit ? 1 : 0);
}

template<int BITS, int ES>
inline constexpr posit<BITS,ES> data_to_posit(typename posit<BITS,ES>::posit_data_t dat) noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  posit<BITS,ES> res{};
  res.m_data = dat & (static_left_shift<posit_data_t, BITS>(1) - 1);
  return res;
}

template<int BITS, int ES>
constexpr posit<BITS,ES>::posit() noexcept : m_data{0} { }

template<int BITS, int ES>
constexpr posit<BITS,ES>::posit(const posit& other) noexcept : m_data{other.m_data} { }

template<int BITS, int ES>
constexpr posit<BITS,ES>::posit(const std::bitset<BITS>& val) noexcept : m_data{static_cast<posit_data_t>(val)} { }

template<int BITS, int ES>
template<typename Integer, std::enable_if_t<std::is_integral<Integer>::value, bool>>
constexpr posit<BITS,ES>::posit(const Integer val) noexcept
  : m_data{0}
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  // check exception cases
  if(val == 0) return;
  if(val == static_shift<Integer,true,sizeof(Integer)>(1)) return;

  // compute ES and regime using abs(val)
  bool sign = val < 0;
  Integer abs_val = sign ? -val : val;
  posit_data_t val_bitwidth = 8 * sizeof(Integer);
  posit_data_t leading_zeroes = std::__countl_zero(abs_val);
  posit_data_t power = val_bitwidth - leading_zeroes;

  posit_data_t regime_coef = static_shift<posit_data_t,true,ES>(1);
  posit_data_t regime = power / regime_coef;
  posit_data_t exponent = power - regime * regime_coef;

  posit_data_t regime_width = regime + 1;
  posit_data_t regime_mask = (static_shift<posit_data_t,true,regime_width>(1) - 1);

  posit_data_t fraction_bits = BITS - 2 - regime_width - ES;
  posit_data_t fraction_int_mask = static_shift<posit_data_t,true,fraction_bits>(1) - 1;
  posit_data_t fraction_int_mask_offset = static_shift<posit_data_t,true,(val_bitwidth - leading_zeroes)>(fraction_int_mask);

  posit_data_t u_bits_of_int = static_shift<posit_data_t,false,(val_bitwidth - leading_zeroes)>(fraction_int_mask_offset & val);
  posit_data_t w_bits_of_int = u_bits_of_int + 1;
  posit_data_t v_bits_of_int = (u_bits_of_int << 1) | 1;
  posit_data_t v_bits_of_int_shifted = static_shift<posit_data_t,true,(val_bitwidth - leading_zeroes - 1)>(v_bits_of_int);

  bool u_is_even = ((u_bits_of_int & 0x1) == 0);
  bool x_is_lt_v = v_bits_of_int_shifted < val;
  bool x_is_eq_v = v_bits_of_int_shifted == val;

  posit_data_t fraction = (u_is_even && x_is_eq_v) || x_is_lt_v ? u_bits_of_int : v_bits_of_int;
  
  bool include_exponent = (BITS - 2 - regime_width) >  0;
  bool include_full_exp = (BITS - 2 - regime_width) >= ES;
  bool include_fraction = (BITS - 2 - regime_width - ES) > 0;

  constexpr posit_data_t regime_shifted     = static_shift<posit_data_t,true ,      (BITS - 1 - regime_width      )>(regime_mask);
  constexpr posit_data_t exponent_shifted   = static_shift<posit_data_t,true ,      (BITS - 2 - regime_width - ES )>(exponent   );
  constexpr posit_data_t exponent_shifted_r = static_shift<posit_data_t,false,(ES - (BITS - 2 - regime_width     ))>(exponent   );

  if(include_fraction)
  {
    m_data = regime_shifted
            | exponent_shifted
            | fraction;
    return;
  } 
  if(include_full_exp)
  { 
    m_data = regime_shifted
            | exponent;
    return;
  }
  if(include_exponent)
  {
    m_data = regime_shifted
            | exponent_shifted_r;
    return;
  } 

  m_data = regime_mask;
}

template<int BITS, int ES>
template<typename Float, std::enable_if_t<std::is_floating_point<Float>::value, bool>>
constexpr posit<BITS,ES>::posit(Float val) noexcept
  : m_data{0}
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  // check exception cases
  if constexpr (std::isnan(val)) { m_data = NaR().m_data; return; }
  if constexpr (std::abs(val) == 0) return;

  using v_int_t = posit_internal_repr_t<8 * sizeof(Float)>;

  // decompose the float
  int32_t power;
  constexpr Float significand = std::frexp(val, &power);
  constexpr bool sign = (significand < 0);

  // we're just going to assume that FLT_RADIX is 2
  constexpr int32_t significand_bits = std::numeric_limits<Float>::digits;
  constexpr v_int_t significand_mask = static_left_shift<posit_data_t,significand_bits>(1) - 1;
  constexpr v_int_t significand_int  = significand_mask & bit_cast<v_int_t>(significand);

  constexpr int32_t significand_msb = std::__countl_zero(significand_int);
  constexpr int32_t significand_pwr = significand_bits - significand_msb - 1;
  constexpr int32_t total_power = -significand_pwr + power;

  constexpr int32_t regime_coef = static_left_shift<posit_data_t,ES>(1);
  constexpr int32_t regime_test = total_power / regime_coef;
  constexpr int32_t exponent_test = total_power - regime_test * regime_coef;
  constexpr int32_t regime = exponent_test >= 0 ? regime_test : regime_test - 1;
  constexpr posit_data_t exponent = total_power - regime * regime_coef;

  constexpr int32_t regime_width = regime < 0 ? -regime : regime + 1;
  constexpr posit_data_t regime_mask = regime < 0 ? 1 : left_shift<posit_data_t>(1, regime_width)-1;

  constexpr posit_data_t fraction_bits = BITS - 2 - regime_width - ES;
  constexpr posit_data_t fraction_int_mask = left_shift<posit_data_t>(1, fraction_bits) - 1;
  constexpr posit_data_t fraction_int_mask_offset = left_shift<posit_data_t>(fraction_int_mask, (significand_bits - significand_msb));

  constexpr posit_data_t u_bits_of_int = right_shift<posit_data_t>(fraction_int_mask_offset & significand_int, (significand_bits - significand_msb));
  constexpr posit_data_t w_bits_of_int = u_bits_of_int + 1;
  constexpr posit_data_t v_bits_of_int = (u_bits_of_int << 1) | 1;
  constexpr posit_data_t v_bits_of_int_shifted = left_shift<posit_data_t>(v_bits_of_int,(significand_bits - significand_msb - 1));

  constexpr bool u_is_even = ((u_bits_of_int & 0x1) == 0);
  constexpr bool x_is_lt_v = v_bits_of_int_shifted < significand_int;
  constexpr bool x_is_eq_v = v_bits_of_int_shifted == significand_int;

  constexpr posit_data_t fraction = (u_is_even && x_is_eq_v) || x_is_lt_v ? u_bits_of_int : v_bits_of_int;
  
  constexpr bool include_exponent = (BITS - 2 - regime_width) >  0;
  constexpr bool include_full_exp = (BITS - 2 - regime_width) >= ES;
  constexpr bool include_fraction = (BITS - 2 - regime_width - ES) > 0;

  constexpr posit_data_t regime_shifted     = left_shift<posit_data_t>(regime_mask,      (BITS - 1 - regime_width      ));
  constexpr posit_data_t exponent_shifted   = left_shift<posit_data_t>(exponent   ,      (BITS - 2 - regime_width - ES ));
  constexpr posit_data_t exponent_shifted_r = right_shift<posit_data_t>(exponent  ,(ES - (BITS - 2 - regime_width     )));

  if constexpr (include_fraction)
  {
    m_data = regime_shifted
            | exponent_shifted
            | fraction;
  } 
  else if constexpr (include_full_exp)
  { 
    m_data = regime_shifted
            | exponent;
  }
  else if constexpr (include_exponent)
  {
    m_data = regime_shifted
            | exponent_shifted_r;
  } else {
    m_data = regime_mask;
  }

  if constexpr (sign) m_data = -m_data;
}

//=========================================================
// Explicit casting functions
//=========================================================

template<int BITS, int ES>
template<int _BITS>
constexpr posit<BITS,ES>::operator posit<_BITS,ES>() const noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  constexpr posit_data_t data = m_data;
  if constexpr (_BITS > BITS)
    return data_to_posit(static_left_shift< posit_data_t,(_BITS -  BITS)>(data));
  else
    return data_to_posit(static_right_shift<posit_data_t,( BITS - _BITS)>(data));
}

template<int BITS, int ES>
template<int _ES>
constexpr posit<BITS,ES>::operator posit<BITS,_ES>() const noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  constexpr posit_data_t data_bit_length = sizeof(posit_data_t) * 8;

  // calculate new regime run & exponent value
  constexpr posit_data_t old_power = power();
  constexpr posit_data_t new_regime_coef = static_shift<posit_data_t,true,_ES>(1);
  constexpr posit_data_t new_regime_value = old_power / new_regime_coef;

  // regime run lengths are -value if vlaue < 0 else 1 + value
  // regime runs are 0's if new_regime_value < 0 else 1's
  constexpr posit_data_t new_regime_run = new_regime_value < 0 ? -new_regime_value : new_regime_value + 1;
  constexpr posit_data_t new_regime_mask = new_regime_value < 0
                                          ? 0 
                                          : static_shift<posit_data_t,false, (data_bit_length - BITS + 1)>(
                                            static_shift<posit_data_t,true , (data_bit_length - new_regime_run)>(~0));

  constexpr posit_data_t new_exp_val = old_power - new_regime_run * new_regime_coef;
  constexpr posit_data_t new_regime_and_exp = new_regime_run | new_exp_val;
  constexpr posit_data_t new_exp_regime_len = 1 + new_regime_run + _ES;
  constexpr posit_data_t new_frac_bitwidth = BITS - new_exp_regime_len;

  if constexpr (new_frac_bitwidth <= 0) {
    if constexpr (m_data < 0) return data_to_posit(-new_regime_run);
    else return data_to_posit(new_regime_run);
  }

  // calculate the new fraction
  constexpr posit<BITS,ES> old_fraction = fraction();
  constexpr posit_data_t old_frac_bitwidth = BITS - old_fraction.regime_bitwidth() - ES - 1; 
  constexpr posit_data_t old_frac_mask = static_shift<posit_data_t,false,(data_bit_length - old_frac_bitwidth)>(~0);
  constexpr posit_data_t frac_width_diff = new_frac_bitwidth - old_frac_bitwidth;
  constexpr posit_data_t new_fraction = frac_width_diff < 0
                                      ? static_shift<posit_data_t,false,-frac_width_diff>(old_frac_mask & old_fraction.m_data)
                                      : static_shift<posit_data_t,true , frac_width_diff>(old_frac_mask & old_fraction.m_data);

  constexpr posit_data_t constructed_data = static_shift<posit_data_t,true,new_frac_bitwidth>(new_regime_and_exp) | new_fraction;

  if constexpr (m_data < 0) return posit<BITS,_ES>{std::bitset<BITS>{-constructed_data}};
  return posit<BITS,_ES>{std::bitset<BITS>{constructed_data}};
}

template<int BITS, int ES>
template<int _BITS, int _ES>
constexpr posit<BITS,ES>::operator posit<_BITS,_ES>() const noexcept
{
  if constexpr (_BITS > BITS) 
  {
    // promote bitwidth before shifting exponent to preserve precision when it matters
    constexpr posit<_BITS,ES> promote_bits = static_cast<posit<_BITS,ES>>(*this);
    return static_cast<posit<_BITS,_ES>>(promote_bits);
  }
  // promote exponent size first in all other cases
  constexpr posit<BITS,_ES> shift_ES = static_cast<posit<BITS,_ES>>(*this);
  return static_cast<posit<_BITS,_ES>>(shift_ES);
}

template<int BITS, int ES>
template<typename Integer, std::enable_if_t<std::is_integral<Integer>::value && std::is_signed<Integer>::value, bool>>
constexpr posit<BITS,ES>::operator Integer() const noexcept
{
  if constexpr (isNaR() || (power() >= 8 * sizeof(Integer)))
    return static_shift<posit_data_t,true,8 * sizeof(Integer) - 1>(1);

  constexpr posit<BITS,ES> significand = significand();
  constexpr posit<BITS,ES> significand_frac_mask = static_shift<posit_data_t,true,(BITS - ES - 2)>(1) - 1;
  constexpr Integer intermediate = power() >= 0 ? (m_data < 0 ? -static_shift<posit_data_t,false,power()>(static_cast<Integer>(-significand.m_data & significand_frac_mask))
                                                              :  static_shift<posit_data_t,false,power()>(static_cast<Integer>( significand.m_data & significand_frac_mask)))
                                                : (m_data < 0 ? -static_shift<posit_data_t,true ,power()>(static_cast<Integer>(-significand.m_data & significand_frac_mask))
                                                              :  static_shift<posit_data_t,true ,power()>(static_cast<Integer>( significand.m_data & significand_frac_mask)));
  if constexpr (intermediate == 0 && m_data != 0)
    return power() >= 0 ? static_shift<posit_data_t,true,8 * sizeof(Integer) - 1>(1)
                        : ((m_data < 0) ? -1 : 1);
  return intermediate;
}

template<int BITS, int ES>
template<typename Float, std::enable_if_t<std::is_floating_point<Float>::value, bool>>
constexpr posit<BITS,ES>::operator Float() const noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  
  constexpr posit_data_t frac_denom = BITS - ES - fraction().regime_bitwidth();
  constexpr posit_data_t frac_mask = left_shift<posit_data_t>(1,frac_denom) - 1;
  constexpr posit_data_t int_significand = (fraction() & frac_mask) + (frac_mask + 1);
  constexpr Float significand = static_cast<Float>(int_significand);
  constexpr Float power = std::exp2(power() + frac_denom);
  return significand * power;
}

template<int BITS, int ES>
constexpr posit<BITS,ES>::operator std::bitset<BITS>() const noexcept
{
  return std::bitset<BITS>(m_data);
}

template<int BITS, int ES>
constexpr posit<BITS,ES>::operator quire_t() const noexcept
{
  return quire_t{*this};
}

template<int BITS, int ES>
constexpr quire<BITS,ES> posit<BITS,ES>::ptoQ() const noexcept
{
  return quire_t{*this};
}

template<int BITS, int ES>
constexpr quire<BITS,ES> ptoQ(const posit<BITS,ES>& obj) noexcept
{
  return quire<BITS,ES>{obj};
}

//--------------------------------------------------
// Convenience Functions
//--------------------------------------------------
// Constants
template<int BITS, int ES>
constexpr posit<BITS,ES> posit<BITS,ES>::minPos() noexcept
{
  return posit<BITS,ES>{0}.next();
}

template<int BITS, int ES>
constexpr posit<BITS,ES> posit<BITS,ES>::maxPos() noexcept
{
  return posit<BITS,ES>{0}.prev().negate();
}

template<int BITS, int ES>
constexpr posit<BITS,ES> posit<BITS,ES>::NaR() noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  return data_to_posit(static_left_shift<posit_data_t,BITS-1>(1));
}

template<int BITS, int ES>
constexpr typename posit<BITS,ES>::posit_data_t posit<BITS,ES>::pIntMax() noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  constexpr posit_data_t raw_power = power_of_pIntMax<BITS,ES>::value;
  constexpr posit_data_t power = raw_power < 1 ? 1 : raw_power;
  return static_left_shift<posit_data_t,power>(1);
}

template<int BITS, int ES>
constexpr typename posit<BITS,ES>::posit_data_t posit<BITS,ES>::pIntMin() noexcept
{
  return -pIntMax();
}

template<int BITS, int ES>
constexpr int32_t posit<BITS,ES>::precision() noexcept
{
  return BITS;
}

// Methods
template<int BITS, int ES>
constexpr bool posit<BITS,ES>::isNaR() const noexcept
{
  return m_data == NaR().m_data;
}

template<int BITS, int ES>
constexpr typename posit<BITS,ES>::posit_data_t 
posit<BITS,ES>::regime() const noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  if(isNaR()) return m_data;
  constexpr posit_data_t regime_coef = static_left_shift<posit_data_t,ES>(2);
  constexpr posit_data_t regime_on_left = static_left_shift<posit_data_t,(8*sizeof(posit_data_t)-BITS)>(m_data);

  constexpr bool regime_is_ones = regime_on_left < 0;
  constexpr posit_data_t regime_bits = regime_is_ones ? std::__countl_one(regime_on_left)
                                                      : std::__countl_zero(regime_on_left);

  constexpr posit_data_t regime_value = regime_is_ones  ?  regime_value - 1 
                                                        : -regime_value;

  return regime_value * regime_coef;
}

template<int BITS, int ES>
constexpr typename posit<BITS,ES>::posit_data_t 
posit<BITS,ES>::regime_bitwidth() const noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  if(isNaR()) return m_data;
  constexpr posit_data_t regime_on_left = static_left_shift<posit_data_t,(8*sizeof(posit_data_t)-BITS)>(m_data);

  constexpr bool regime_is_ones = regime_on_left < 0;
  constexpr posit_data_t regime_bits = regime_is_ones ? std::__countl_one(regime_on_left)
                                                      : std::__countl_zero(regime_on_left);

  return regime_bits;
}

template<int BITS, int ES>
constexpr typename posit<BITS,ES>::posit_data_t 
posit<BITS,ES>::exponent() const noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  if(isNaR()) return m_data;
  constexpr posit_data_t regime_on_left = static_shift<posit_data_t,true,(8*sizeof(posit_data_t)-BITS)>(m_data);

  constexpr bool regime_is_ones = regime_on_left < 0;
  constexpr posit_data_t regime_bits = regime_is_ones ? std::__countl_one(regime_on_left)
                                                      : std::__countl_zero(regime_on_left);

  constexpr posit_data_t exponent_on_left = left_shift<posit_data_t>(regime_on_left,(regime_bits+1));
  constexpr posit_data_t exponent_on_right = static_right_shift<posit_data_t,(8*sizeof(posit_data_t) - ES)>(exponent_on_left);

  return exponent_on_right;
}

template<int BITS, int ES>
constexpr typename posit<BITS,ES>::posit_data_t 
posit<BITS,ES>::power() const noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  if(isNaR()) return m_data;
  constexpr posit_data_t regime_coef = static_shift<posit_data_t,true,ES>(2);
  constexpr posit_data_t regime_on_left = static_shift<posit_data_t,true,(8*sizeof(posit_data_t)-BITS)>(m_data);

  constexpr bool regime_is_ones = regime_on_left < 0;
  constexpr posit_data_t regime_bits = regime_is_ones ? std::__countl_one(regime_on_left)
                                                      : std::__countl_zero(regime_on_left);

  constexpr posit_data_t regime_value = regime_is_ones  ?  regime_bits - 1 
                                                        : -regime_bits;

  constexpr posit_data_t exponent_on_left = left_shift<posit_data_t>(regime_on_left, (regime_bits+1));
  constexpr posit_data_t exponent_on_right = static_right_shift<posit_data_t,(8*sizeof(posit_data_t) - ES)>(exponent_on_left);

  return m_data < 0 ? -(regime_value * regime_coef + exponent_on_right + 1)
                    :  (regime_value * regime_coef + exponent_on_right);
}

template<int BITS, int ES>
constexpr posit<BITS,ES>
posit<BITS,ES>::rtz() const noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  if(isNaR() || m_data == 0) return *this;
  constexpr posit_data_t lower_frac_bitwidth = BITS - regime_bitwidth() - ES + 1 - power();
  // if the power is greater than the fraction denominator
  if constexpr (lower_frac_bitwidth < 0) return *this;
  // if the power is very small
  if constexpr (lower_frac_bitwidth >= BITS - 3) return posit<BITS,ES>{0};
  constexpr posit_data_t lower_bit_mask = left_shift<posit_data_t>(2, lower_frac_bitwidth) - 1;
  constexpr posit_data_t modified_posit = m_data & ~lower_bit_mask;
  return data_to_posit(modified_posit);
}

template<int BITS, int ES>
constexpr posit<BITS,ES> 
posit<BITS,ES>::fraction() const noexcept
{
  if(isNaR() || m_data == 0) return *this;
  return (*this - this->rtz()).abs();
}

template<int BITS, int ES>
constexpr posit<BITS,ES>
posit<BITS,ES>::significand() const noexcept
{
  if(isNaR() || m_data == 0) return *this;
  return sign() * fraction() + sign();
}

template<int BITS, int ES>
constexpr posit<BITS,ES>
posit<BITS,ES>::sign() const noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  if(isNaR() || m_data == 0) return *this;
  constexpr posit_data_t one = static_left_shift<posit_data_t, BITS - 2>(1);
  constexpr posit_data_t n_1 = static_left_shift<posit_data_t, BITS - 2>(3);
  constexpr posit_data_t result = m_data < 0 ? n_1 : one;
  return data_to_posit(result);
}

template<int BITS, int ES>
constexpr posit<BITS,ES>
posit<BITS,ES>::reciprocal() const noexcept
{
  using posit_data_t = typename posit<BITS,ES>::posit_data_t;
  if(isNaR() || m_data == 0) return NaR();

  // a posit is (1-3s + f) * 2^p
  // where p = (1-2s) * ((2^ES)*r+e+s)
  // and f = (2^-m) * (unsigned int F of fraction bits)

  // the reciprocal of a posit is the rounded result of:
  // 1 / ((1-3s + f) * 2^p)
  // = (2^-p) / ((1-3s) + f)
  // = (2^-p) / ((1-3s) + (2^-m)*(F))
  // = (2^(m-p)) / ((1-3s)*2^m + F)
  // = 2^P / D

  constexpr posit_data_t s = m_data < 0 ? -1 : 1;
  constexpr posit_data_t m_test = BITS - ES - regime_bitwidth() - 1;
  constexpr posit_data_t m = m_test < 0 ? 0 : m_test;
  constexpr posit_data_t p = power();

  constexpr posit_data_t f_mask = left_shift<posit_data_t>(1,m) - 1;
  constexpr posit_data_t F = m_data & f_mask;

  // numerator N = 2^P
  constexpr posit_data_t P = m - power();
  constexpr posit_data_t D = left_shift<posit_data_t>(1 - 3*s, m) + F;

  constexpr posit_data_t abs_D = D < 0 ? -D : D;
  // smallest power of 2 greater or equal to D
  constexpr posit_data_t ceil_D = std::__bit_ceil(abs_D);
  constexpr posit_data_t ceil_powD = std::__countr_zero(ceil_D);
  // largest power of 2 not greater than D
  constexpr posit_data_t floor_D = std::__bit_floor(abs_D);
  constexpr posit_data_t floor_powD = std::__countr_zero(floor_D);

  if constexpr (ceil_D == floor_D) return data_to_posit<BITS,ES>(power_to_regime<BITS,ES>(P - floor_D));

  // need to pick the posits above & below (2^P)/D
  // so p_0 <= 2^P/D <= p_1
  // and p_01 = p_0 append 1

  constexpr posit_data_t upper_est_pow = P - floor_powD;
  constexpr posit_data_t lower_est_pow = P - ceil_powD;

  // obvious conditions:
  // if the num - (underestimate denom) is less than the smallest representable power, the true result is as well.
  // if the num - (overestimate denom) is greater than the smallest representable power, the true result is as well.
  if constexpr (upper_est_pow <= posit<BITS,ES>::minPos().power()) 
    return m_data < 0 ? posit<BITS,ES>::minPos().negate() : posit<BITS,ES>::minPos();
  if constexpr (lower_est_pow >= posit<BITS,ES>::maxPos().power()) 
    return m_data < 0 ? posit<BITS,ES>::maxPos().negate() : posit<BITS,ES>::maxPos();

  using u_posit_data_t = std::make_unsigned_t<posit_data_t>;

  // set the MSB to 1, then do division of the surrogate by the fraction
  constexpr u_posit_data_t pow2_surrogate = ~(~0 >> 1);
  constexpr posit_data_t surrogate_div_res = pow2_surrogate / abs_D;

  constexpr posit_data_t surrogate_pow_complement = std::__countl_zero(surrogate_div_res);
  constexpr posit_data_t surrogate_res_bitwidth = BITS - surrogate_pow_complement;
  constexpr posit_data_t surrogate_frac = surrogate_div_res ^ left_shift<posit_data_t>(1, surrogate_res_bitwidth);

  constexpr posit_data_t final_pow = P - surrogate_pow_complement;
  constexpr posit_data_t result_pow_encoding = power_to_bit_encoding(final_pow);
  constexpr auto result_pow_encoding_posit = data_to_posit(result_pow_encoding);
  if constexpr ((result_pow_encoding_posit == posit<BITS,ES>::minPos()) || (result_pow_encoding_posit == posit<BITS,ES>::maxPos())) return result_pow_encoding_posit;

  constexpr posit_data_t bits_of_frac_in_encoding = BITS - 2 - result_pow_encoding_posit.regime_bitwidth() - ES;
  constexpr posit_data_t bits_of_frac_in_encoding_capped = bits_of_frac_in_encoding < 1 ? 0 : bits_of_frac_in_encoding;
  constexpr posit_data_t rounding_bit_idx = surrogate_res_bitwidth - 1 - bits_of_frac_in_encoding_capped;
  constexpr posit_data_t rounding_bit = (surrogate_frac & left_shift<posit_data_t>(1, rounding_bit_idx)) > 0 ? 1 : 0;

  constexpr posit_data_t frac_encoding_shift_test = surrogate_res_bitwidth - bits_of_frac_in_encoding_capped;
  constexpr posit_data_t frac_encoding_shift = frac_encoding_shift_test < 0 ? 0 : frac_encoding_shift_test;
  constexpr posit_data_t frac_encoding = right_shift<posit_data_t>(surrogate_frac, frac_encoding_shift);
  constexpr posit_data_t result_full_encoding = (result_pow_encoding | frac_encoding) + rounding_bit;

  return data_to_posit(result_full_encoding);
}

//===============================================
// Misc Structs
//===============================================

template<int32_t BITS, int32_t ES>
struct power_of_pIntMax<BITS,ES,std::enable_if_t<
  (ceil_log2<BITS-ES-1>() <= ES)
>>
{
  static constexpr int32_t value = BITS-ES-2;
};

template<int32_t BITS, int32_t ES>
struct power_of_pIntMax<BITS,ES,std::enable_if_t<
  (ceil_log2<BITS-ES-1>() > ES)
>>
{
  static constexpr int32_t value = (BITS-ES-1) - ((BITS-ES-1) + (1 << ES)) / ((1 << ES) + 1);
};


template<int32_t BITS, int32_t ES_TEST, int32_t DIR>
struct count_until_decr<BITS,ES_TEST,DIR,std::enable_if_t<
      (ES_TEST+DIR < BITS) && (ES_TEST+DIR > 0)
  &&  (posit<BITS, ES_TEST+DIR>::pIntMax() >= posit<BITS, ES_TEST>::pIntMax())
>>
{
  static constexpr int32_t value = count_until_decr<BITS, ES_TEST+DIR, DIR>::value;
};

template<int32_t BITS, class V>
struct count_until_decr<BITS,1,-1,V> 
{
  static constexpr int32_t value = posit<BITS, 0>::pIntMax() >= posit<BITS, 1>::pIntMax() ? 0 : 1;
};

template<int32_t BITS, int32_t ES_TEST, int32_t DIR, class V>
struct count_until_decr 
{
  static constexpr int32_t value = ES_TEST;
};

template<int32_t BITS>
struct posit_max_consec_int_range_t
{
  static constexpr int32_t B = count_until_decr<BITS,0, 1>::value;
  static constexpr int32_t A = count_until_decr<BITS,B,-1>::value;

  using lower_bound_type = posit<BITS,A>;
  using upper_bound_type = posit<BITS,B>;
};


} // end namespace posit