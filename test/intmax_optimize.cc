#include <cmath>
#include <iostream>
#include "posit.hh"

std::ostream&
operator<<( std::ostream& dest, __int128_t value )
{
    std::ostream::sentry s( dest );
    if ( s ) {
        __uint128_t tmp = value < 0 ? -value : value;
        char buffer[ 128 ];
        char* d = std::end( buffer );
        do
        {
            -- d;
            *d = "0123456789"[ tmp % 10 ];
            tmp /= 10;
        } while ( tmp != 0 );
        if ( value < 0 ) {
            -- d;
            *d = '-';
        }
        int len = std::end( buffer ) - d;
        if ( dest.rdbuf()->sputn( d, len ) != len ) {
            dest.setstate( std::ios_base::badbit );
        }
    }
    return dest;
}

template<int i>
struct print;

template<int i>
struct print
{
  static void run()
  {
    posit::posit_max_consec_int_range_t<i> pair;
    posit::posit<i,pair.A> val{};
    std::cout << i << ' ' 
              << pair.A << ' ' 
              << pair.B << ' ' 
              << pair.B - pair.A << ' '
              << posit::power_of_pIntMax<i,pair.A>::value << ' '
              << val.pIntMax() << '\n';
    print<i+1>::run();
  }
};

template<>
struct print<129>
{
  static void run()
  {  }
};

int main()
{
  print<4>::run();
}