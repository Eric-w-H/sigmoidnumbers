#include "posit.hh"

#define BITS 8
#define ES_MAX (BITS-1-std::__bit_floor(BITS-1))

template<uint32_t ES>
void print_posits(const std::bitset<BITS>& bits)
{
  using posit_t = posit::posit<BITS, ES>;
  posit_t posit{bits};
  std::cout << posit << ',';
  print_posits<ES+1>(bits);
}

template<>
void print_posits<ES_MAX>(const std::bitset<BITS>& bits)
{
  using posit_t = posit::posit<BITS, BITS-2>;
  posit_t posit{bits};
  std::cout << posit << '\n';
}


int main()
{
  for(int i = 0; i <= ES_MAX; ++i)
    std::cout << ',' << i;
  std::cout << '\n';

  for(int i = 0; i < (1 << BITS); ++i)
  {
    auto bits = std::bitset<BITS>(i);
    std::cout << bits << ',';
    print_posits<0>(bits);
  }
}

