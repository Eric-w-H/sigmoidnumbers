#!/bin/bash

clang++ -Wall -Werror -Wextra --std=c++20 -O3 -I ./src/ "test/$1.cc" -o build/$1