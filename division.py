
# bits=8
# es=2

# for k in range(1,bits - 2):
#   for e in range(0, 2**min(es,bits-2-k)):
#     m = bits - es - 2 - k
    
#     p = (2**es)*(k-1) + e
#     for f in range(0, int(2**m)):
#       print(f'F({p},{m},{f}) = 2^{m-p} / (2^{m} + {f})')

# print('---------')
# for k in range(1,bits - 2):
#   for e in range(0, 2**min(es,bits-2-k)):
#     m = bits - es - 2 - k

#     p = (2**es)*(-k) + e
#     for f in range(0, int(2**m)):
#       print(f'F({p},{m},{f}) = 2^{m-p} / (2^{m} + {f})')

def countl_zero(n:int, bits:int):
  zeros = 0
  mask = 1 << (bits-1)
  while(n & mask == 0):
    zeros += 1
    n <<= 1
  return zeros

def floor_log2(n:int):
  zeros = 0
  while(n > 0):
    zeros += 1
    n >>= 1
  return zeros

bits = 8
n = 1 << (bits - 1);
d = 3
print(f'{hex(n)} / {hex(d)} = {hex(n//d)} [{floor_log2(n)} ; {floor_log2(d)} => {countl_zero(n//d, bits)}]')