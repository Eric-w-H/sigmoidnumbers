#include "posit.hh"

template<uint8_t BITS, uint8_t ES>
void print_all_posits()
{
  std::cout << "BITS: " << +BITS << ", ES: " << +ES << '\n';
  for(int i = 0; i < (1<<BITS); ++i)
  {
    auto bits = std::bitset<BITS>(i);
    auto posit = posit::posit<BITS, ES>{bits};
    std::cout << bits << " => " << posit << '\n';
    // posit.print_fraction_derivation(std::cout) << '\n';
  }
}

int main()
{
  print_all_posits<8,6>();
}
