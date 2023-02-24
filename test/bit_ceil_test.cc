#include <bit>
#include <iostream>
#include <bitset>

int qsize(int b, int e, int c)
{
  return std::__bit_ceil( (4*b - 8)*(1<<e) + 1 + c);
}

int main() 
{

  for(int b = 4; b < 64; b*=2)
    for(int e = 0; e < b/4; ++e )
      for(int c = 10; c < 31; c += 10)
        std::cout << b << ' ' 
                  << e << ' '  
                  << c << ' ' 
                  << qsize(b, e, c) << '\n';
}