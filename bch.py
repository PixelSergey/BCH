from sympy import GF, ZZ, Matrix, poly, Poly, div, Symbol, Add, Mul, Pow, Integer
from sympy.abc import x, z
from sympy.polys.galoistools import gf_gcdex, gf_lcm
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


def expand_expression(expression, var):
    if isinstance(expression, Symbol):
        return Poly(expression, var, domain=GF(2))

    if isinstance(expression, Integer):
        return Poly(expression, var, domain=GF(2))

    if isinstance(expression, Mul):
        result = Poly(1, var, domain=GF(2))
        for term in expression.args:
            result *= expand_expression(term, var)
        return result
    
    if isinstance(expression, Add):
        result = Poly(0, var, domain=GF(2))
        for term in expression.args:
            result += expand_expression(term, var)
        return result
    
    if isinstance(expression, Pow):
        base, exponent = expression.args
        base = expand_expression(base, var)
        exponent = int(exponent)
        if exponent < 0:
            base = find_inverse(base)
            exponent = abs(exponent)

        return Poly(base**exponent, var, domain=GF(2))


def find_minimal_polynomial(element, reducing):
    seen = set()
    i = 0
    result = Poly(1, x)
    while True:
        current = Poly(element**(2**i), z) % reducing
        if current in seen:
            break
        result *= Poly(x, x) - current.set_domain(ZZ)
        seen.add(current)
        i += 1

    result = Poly(result, x, domain=GF(2)[z])
    reduced = [expand_expression(coefficient, z)%reducing for coefficient in result.all_coeffs()]
    assert(all([coefficient==0 or coefficient==1 for coefficient in reduced]))

    result = Poly(reduced, x, domain=GF(2))
    return result


def find_generator(alpha, reducing, t):
    generator = Poly(1, x, domain=GF(2))
    for i in range(1, 2*t):
        current = find_minimal_polynomial(alpha**i, reducing)
        generator = Poly(gf_lcm(generator.all_coeffs(), current.all_coeffs(), 2, ZZ), x, domain=GF(2))
    return generator


def polynomify(integer):
    return Poly([int(x) for x in list(bin(integer)[2:])], x, domain=GF(2))

def integerify(polynomial):
    return sum([coeff*2**i for i,coeff in enumerate(polynomial.all_coeffs()[::-1])])


def substitute(polynomial, substitution):
    result = poly(0, x, domain=GF(2))
    for i, coeff in enumerate(polynomial.all_coeffs()[::-1]):
        if coeff == 0:
            continue
        result += Poly(substitution**i, x, domain=GF(2))
    return result


def find_inverse(polynomial):
    inv, _, gcd = gf_gcdex(polynomial.all_coeffs(), generator.all_coeffs(), 2, ZZ)
    assert(gcd == [1])
    return Poly(inv, x, domain=GF(2))


def find_all_powers(element, reducing):
    result = dict()
    for i in range(1, 16):
        power = Poly(element**i, z, domain=GF(2)) % reducing
        if tuple(power.all_coeffs()) in result:
            continue
        result[tuple(power.all_coeffs())] = i
    return result


def encode_bch(bits):
    plaintext = Poly(bits, x, domain=GF(2))
    encoded = plaintext * generator
    return encoded.all_coeffs()


def find_error_locator(syndromes):
    for i in range(t):
        nu = t-i  # Number of errors
        syndrome_matrix = Matrix(nu, nu, lambda x,y: syndromes[x+y])
        detection = syndrome_matrix.det() % reducing
        if detection == 0:
            continue

        syndrome_vector = Matrix(nu, 1, lambda x,_: -syndromes[nu+x])
        augmented = syndrome_matrix.col_insert(i, syndrome_vector)
        locator = augmented.rref(pivots=False).col(nu)
        result = []
        for row in range(nu):
            result.append(expand_expression(locator[row]))

        return result


def find_error_pos(locator):
    print(locator)
    locator_poly = Poly(1, x, domain=GF(2))
    for i, lambda_i in enumerate(locator[::-1], start=1):
        locator_poly += lambda_i * Poly(x**i, x, domain=GF(2))
    
    powers = []
    for power in alphapowers:
        root = substitute(locator_poly, Poly(power, x, domain=GF(2))) % reducing
        if root == 0:
            inverse = find_inverse(Poly(power, x, domain=GF(2)))
            powers.append(alphapowers[tuple((inverse%reducing).all_coeffs())])

    print(powers)
    return powers


def decode_bch(bits):
    encoded = Poly(bits, x, domain=GF(2))

    syndromes = []
    for i in range(1,2*t+1):
        syndrome = substitute(encoded, x**i) % reducing
        syndromes.append(syndrome)

    if all((syndrome == 0 for syndrome in syndromes)):
        decoded, _ = div(encoded, generator)
        return decoded.all_coeffs()
    
    locator = find_error_locator(syndromes)
    errors = find_error_pos(locator)


if __name__ == "__main__":
    m = 3
    n = 2**m - 1
    t = 2
    
    reducing = Poly(z**4 + z + 1, domain=GF(2))
    alpha = Poly(z, z, domain=GF(2))

    generator = find_generator(alpha, reducing, t)
    print(generator)

    alphapowers = find_all_powers(alpha, reducing)
    # print(alphapowers)
    # print(alphapowers)
    # for power in alphapowers:
    #     print(alphapowers[power], Poly(power, x, domain=GF(2)), substitute(generator, Poly(power, x, domain=GF(2)))%reducing)

    #print(encode_bch([1,0,1,1,1,0]))
    #print(decode_bch([0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0]))
    #print(decode_bch([0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0]))
    #                 ^15th pos      ^--^ errors (9 & 10)          ^ 0th position
    #print(decode_bch([0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0]))
    #                          ^ 12           ^ 7


    # correct = [0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0]
    # # Test all 1-bit errors
    # for i in range(len(correct)-1):
    #         incorrect = correct[:i]+[1-correct[i]]+correct[i+1:]
    #         decode_bch(incorrect)
    # print("---")
    # # Test all 2-bit errors
    # for i in range(len(correct)-1):
    #     for j in range(i+1,len(correct)):
    #         incorrect = correct[:i]+[1-correct[i]]+correct[i+1:j]+[1-correct[j]]+correct[j+1:]
    #         decode_bch(incorrect)
    
