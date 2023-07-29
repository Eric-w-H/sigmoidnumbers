#include "posit.hh"

template<uint8_t BITS, uint8_t ES>
void print_all_posits()
{
  std::cout << "BITS: " << +BITS << ", ES: " << +ES << '\n';
  for(int i = 0; i < (1<<BITS); ++i)
  {
    auto bits = std::bitset<BITS>(i);
    auto posit = posit::posit<BITS, ES>{bits};
    std::cout << bits << " => " << +posit.raw() << " @ " << posit << '\n';
    // posit.print_fraction_derivation(std::cout) << '\n';
  }
}

template<uint8_t BITS>
void print_max_posit()
{
  print_all_posits<BITS, posit::max_exponent<BITS>()>();
}

int main()
{
  // print_all_posits<6,4>();
  print_max_posit<8>();
}
