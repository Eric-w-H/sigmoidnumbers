#include <Posits.hh>

#define BITS 16

template<uint32_t ES>
void print_posits(const std::bitset<BITS>& bits)
{
  using posit_t = Posit<BITS, ES>;
  auto posit = posit_t::from_bits(bits);
  std::cout << posit << ',';
  print_posits<ES+1>(bits);
}

template<>
void print_posits<BITS-2>(const std::bitset<BITS>& bits)
{
  using posit_t = Posit<BITS, BITS-2>;
  auto posit = posit_t::from_bits(bits);
  std::cout << posit << '\n';
}


int main()
{
  for(int i = 0; i < BITS-1; ++i)
    std::cout << ',' << i;
  std::cout << '\n';

  for(int i = 0; i < (1 << BITS); ++i)
  {
    auto bits = std::bitset<BITS>(i);
    std::cout << bits << ',';
    print_posits<0>(bits);
  }
}

