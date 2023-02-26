#include "util.hh"

template<int iters, class func>
struct Unroll
{
  constexpr inline static void run(func& body) { body(iters); Unroll<iters-1, func>::run(body); }
};

template<class func>
struct Unroll<0, func>
{
  constexpr inline static void run(func& body) { body(0); }
};

template<int Value>
constexpr inline int ceil_log2()
{
  return 1 + ceil_log2<(Value >> 1) + 1>();
}

template<>
constexpr inline int ceil_log2<2>()
{
  return 1;
}

template<>
constexpr inline int ceil_log2<1>()
{
  return 0;
}

template<class T, bool Left, int32_t exp>
constexpr inline T static_shift(const T init)
{
  if constexpr (exp <= 0) return init;
  if constexpr (Left)
    return static_shift<T,true,(exp>>5)>(init << exp);
  else
    return static_shift<T,false,(exp>>5)>(init >> exp);
}

template<class T, int32_t exp>
constexpr inline T static_left_shift(const T init)
{ return static_shift<T, true, exp>(init); }
template<class T, int32_t exp>
constexpr inline T static_right_shift(const T init)
{ return static_shift<T, false, exp>(init); }

template<class T, bool Left>
inline T shift(T init, int32_t exp)
{
  for(; exp > 0; exp >>= 5)
    init = Left ? init << exp : init >> exp;
  return init;
}

template<class T>
inline T left_shift(const T init, const int32_t exp)
{ return shift<T,true>(init,exp); }
template<class T>
inline T right_shift(const T init, const int32_t exp)
{ return shift<T,false>(init,exp); }
