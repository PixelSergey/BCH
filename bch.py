from sympy import GF, ZZ, Matrix, poly, Poly, div
from sympy.abc import x
from sympy.polys.galoistools import gf_gcdex

m = 3
n = 2**m - 1

t = 2
generator = poly(x**8 + x**7 + x**6 + x**4 + 1, domain=GF(2))
alpha = poly(x**4 + x**3 + x**2 + x + 1, domain=GF(2))


def polynomify(integer):
    return Poly([int(x) for x in list(bin(integer)[2:])], x, domain=GF(2))

def integerify(polynomial):
    return sum([coeff*2**i for i,coeff in enumerate(polynomial.all_coeffs()[::-1])])


def substitute(polynomial, substitution):
    result = poly(0, x, domain=GF(2))
    for i, coeff in enumerate(polynomial.all_coeffs()[::-1]):
        if coeff == 0:
            continue
        result += substitution**i
    result %= generator
    return result


def find_inverse(polynomial):
    inv, _, gcd = gf_gcdex(polynomial.all_coeffs(), generator.all_coeffs(), 2, ZZ)
    assert(gcd == [1])
    return Poly(inv, x, domain=GF(2))


def expand_division(fraction):
    numerator = Poly(fraction.args[1], x, domain=GF(2))
    denominator = Poly(1/fraction.args[0], x, domain=GF(2))
    inv_denominator = find_inverse(denominator)
    return (numerator * inv_denominator) % alpha


def find_all_roots(polynomial, mod=generator):
    roots = []
    for i in range(16):
        pi = polynomify(i)
        substitution = substitute(polynomial, pi) % mod
        if substitution == 0:
            roots.append(pi)
    
    return roots


def encode_bch(bits):
    plaintext = Poly(bits, x, domain=GF(2))
    encoded = plaintext * generator
    return encoded.all_coeffs()


def find_error_locator(syndromes):
    for i in range(t):
        syndrome_matrix = Matrix(t-i, t-i, lambda x,y: syndromes[x+y])
        detection = syndrome_matrix.det() % generator
        if detection == 0:
            continue

        syndrome_vector = Matrix(t-i, 1, lambda x,_: -syndromes[t+1+x])
        augmented = syndrome_matrix.col_insert(i, syndrome_vector)
        locator = augmented.rref()[0].col(t-i)
        result = []
        for row in range(t-i):
            result.append(expand_division(locator[row]))

        return result


def find_error_pos(locator):
    print(locator)
    locator_poly = Poly(1, x, domain=GF(2))
    for i, lam in enumerate(locator[::-1]):
        locator_poly += Poly(lam**(i+1), x, domain=GF(2))

    print(locator_poly)
    roots = find_all_roots(locator_poly, alpha)
    print(roots)
    
    pos = []
    for i in range(1,16):
        ai = Poly(alpha**i, x, domain=GF(2))
        if substitute(locator_poly, ai) == 0:
            pos.append(i)
    print(pos)
    print([integerify(find_inverse(polynomify(a)%alpha)%alpha) for a in pos])
    return pos


def decode_bch(bits):
    encoded = Poly(bits, x, domain=GF(2))

    syndromes = []
    for root in generator_roots:
        syndrome = substitute(encoded, root) % generator
        if syndrome == 0:
            continue
        syndromes.append(syndrome)

    if not syndromes:
        decoded, _ = div(encoded, generator)
        return decoded.all_coeffs()
    
    locator = find_error_locator(syndromes)
    errors = find_error_pos(locator)


if __name__ == "__main__":
    generator_roots = find_all_roots(generator)

    #print(encode_bch([1,0,1,1,1,0]))
    #print(decode_bch([0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0]))
    print(decode_bch([0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0]))
    #                 ^15th pos      ^--^ errors (9 & 10)          ^ 0th position
