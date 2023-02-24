#include "Posits.hh"

template<uint8_t BITS, uint8_t ES>
void print_all_posits()
{
  std::cout << "BITS: " << +BITS << ", ES: " << +ES << '\n';
  for(int i = 0; i < (1<<BITS); ++i)
  {
    auto bits = std::bitset<BITS>(i);
    auto posit = Posit<BITS, ES>::from_bits(bits);
    std::cout << bits << " => " << posit << '\n';
    // posit.print_fraction_derivation(std::cout) << '\n';
  }
}

int main()
{
  using posit8_t = Posit<8, 6>;

  auto bits = std::bitset<8>(0x10);
  auto posit = posit8_t::from_bits(bits);

  std::cout << bits << std::endl;
  posit.print_posit_bitwidths(std::cout) << std::endl;
  posit.print_posit(std::cout) << std::endl;
  posit.print_fraction_derivation(std::cout) << std::endl;

  bits = std::bitset<8>(0x20);
  posit = posit8_t::from_bits(bits);

  std::cout << bits << std::endl;
  posit.print_posit_bitwidths(std::cout) << std::endl;
  posit.print_posit(std::cout) << std::endl;
  posit.print_fraction_derivation(std::cout) << std::endl;

  print_all_posits<8,6>();
}
