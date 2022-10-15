import gmpy2
from gmpy2 import mpfr

gmpy2.get_context().precision=200

def to_bf(num):
  (m, e, p) = num.digits()
  e -= 44
  sig = 1
  if m[0] == '-':
    sig = -1
    m = m[1::]
  digits = [int(m[i*4:i*4+4]) for i in range(11)][::-1]
  z = int(m[44:45])
  if z > 5 or (z == 5 and digits[0] & 1 != 0):
    digits[0] += 1
  return 'BigFloatInc { sign: ' + str(sig) + ', e: ' + str(e) + ', n: 44, m: ' + str(digits) + ' },'

# asin polynomial
one = mpfr(1)
two = mpfr(2)
d1 = one
d2 = two
p1 = one
p2 = one

for i in range(100):
  p1 = p1 * d1
  p2 = p2 * d2
  p = p1 / p2 / (d1 + two)
  print(to_bf(p))
  d1 = d1 + two
  d2 = d2 + two

# atan polynomial
one = mpfr(1)
two = mpfr(2)
d1 = one
ng = True
for i in range(100):
  d1 += two
  if ng:
    print(to_bf(mpfr(-1) / d1))
  else:
    print(to_bf(mpfr(1) / d1))
  ng = not ng

# atan consts
with open("atan.txt", "wt") as f:
  for i in range(10000):
    n = mpfr(i+1)/mpfr(10000)
    f.write(to_bf(gmpy2.atan(n)) + '\n')


# sqrt
with open("sqrt.txt", "wt") as f:
  for i in range(11):
    d = pow(mpfr(10000), mpfr(i))
    for j in range(99):
      n = mpfr(j*100+100)*d
      f.write(to_bf(gmpy2.sqrt(n)) + '\n')

# cbrt
with open("cbrt.txt", "wt") as f:
  for i in range(11):
    d = pow(mpfr(10000), mpfr(i))
    for j in range(99):
      n = mpfr(j*100+100)*d
      f.write(to_bf(pow(n, mpfr(1)/mpfr(3))) + '\n')


# ln_const
with open("ln.txt", "wt") as f:
  for i in range(8192):
    a = mpfr(i+1) / mpfr(10000)
    f.write(to_bf(gmpy2.atanh(a)) + "\n")


# exp
with open("exp.txt", "wt") as f:
  for i in range(999):
    a = mpfr(i+1) / mpfr(1000)
    f.write(to_bf(gmpy2.exp(a)) + "\n")


# sin
with open("sin.txt", "wt") as f:
  for i in range(1571):
    a = mpfr(i+1) / mpfr(1000)
    f.write(to_bf(gmpy2.sin(a)) + "\n")
  f.write('\n\n\n\n\n\n\n\n')
  for i in range(1571):
    a = mpfr(i+1) / mpfr(1000) + gmpy2.const_pi()/mpfr(2)
    f.write(to_bf(gmpy2.sin(a)) + "\n")

# inverse factorial
d = mpfr(1)
f = d
for i in range(50):
  d += mpfr(1)
  f *= d
  print(to_bf(mpfr(1)/f))
  

gmpy2.get_context().precision=32000
rs = gmpy2.random_state(hash(gmpy2.random_state()))
n = [gmpy2.mpfr_random(rs) for i in range(10000)]
f1 = n[0]
one = mpfr(1)
start = round(time.time() * 1000)
for a in n:
  f = one / a

round(time.time() * 1000)-start
