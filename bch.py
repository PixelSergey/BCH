import random
from sympy import GF, ZZ, Matrix, poly, Poly, div, Symbol, Add, Mul, Pow, Integer
from sympy.abc import x, z
from sympy.polys.galoistools import gf_gcdex, gf_lcm
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


def expand_expression(expression, var, generator, reducing):
    if isinstance(expression, Symbol):
        return Poly(expression, var, domain=GF(2)) % reducing

    if isinstance(expression, Integer):
        return Poly(expression, var, domain=GF(2)) % reducing

    if isinstance(expression, Mul):
        result = Poly(1, var, domain=GF(2))
        for term in expression.args:
            result *= expand_expression(term, var, generator, reducing)
        return result % reducing
    
    if isinstance(expression, Add):
        result = Poly(0, var, domain=GF(2))
        for term in expression.args:
            result += expand_expression(term, var, generator, reducing)
        return result % reducing
    
    if isinstance(expression, Pow):
        base, exponent = expression.args
        base = expand_expression(base, var, generator, reducing)
        exponent = int(exponent)
        if exponent < 0:
            assert(generator is not None)
            base = find_inverse(base, reducing, z)
            exponent = abs(exponent)

        return Poly(base**exponent, var, domain=GF(2)) % reducing


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
    reduced = [expand_expression(coefficient, z, None, reducing)%reducing for coefficient in result.all_coeffs()]
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
    return Poly([int(x) for x in list(bin(integer)[2:])], z, domain=GF(2))

def integerify(polynomial):
    return sum([coeff*2**i for i,coeff in enumerate(polynomial.all_coeffs()[::-1])])


def substitute(polynomial, substitution):
    result = poly(0, z, domain=GF(2))
    for i, coeff in enumerate(polynomial.all_coeffs()[::-1]):
        result += Poly(coeff, z, domain=GF(2)) * Poly(substitution**i, z, domain=GF(2))
    return result


def find_inverse(polynomial, reducing, var):
    inv, _, gcd = gf_gcdex(polynomial.all_coeffs(), reducing.all_coeffs(), 2, ZZ)
    if gcd != [1]:
        print(polynomial)
        assert(False)
    return Poly(inv, var, domain=GF(2))


def find_all_powers(element, reducing):
    result = dict()
    for i in range(0, 15):
        power = Poly(element**i, z, domain=GF(2)) % reducing
        if tuple(power.all_coeffs()) in result:
            continue
        result[tuple(power.all_coeffs())] = i
    return result


def find_all_roots(polynomial, alpha, reducing):
    roots = []
    for i in range(1,16):
        root = substitute(polynomial, alpha**i) % reducing
        if root == 0:
            roots.append(alpha**i % reducing)
    return roots


def encode_bch(bits, generator):
    plaintext = Poly(bits, x, domain=GF(2))
    encoded = plaintext * generator
    return encoded.all_coeffs()


def find_error_locator(syndromes, generator, reducing, t):
    for i in range(t):
        nu = t-i  # Number of errors
        syndrome_matrix = Matrix(nu, nu, lambda a,b: syndromes[a+b])
        detection = syndrome_matrix.det() % reducing
        if detection == 0:
            continue
        
        syndrome_vector = Matrix(nu, 1, lambda a,_: -syndromes[nu+a])
        augmented = syndrome_matrix.col_insert(nu, syndrome_vector)
        locator = augmented.rref(pivots=False).col(nu)
        result = []
        for row in range(nu):
            result.append(expand_expression(locator[row], z, generator, reducing))

        return result


def find_error_pos(locator, alpha, generator, reducing):
    locator_poly = Poly(1, x, domain=GF(2)[z])
    for i, lambda_i in enumerate(locator[::-1], start=1):
        locator_poly += Poly(lambda_i%reducing, x, domain=GF(2)[z]) * Poly(x**i, x, domain=GF(2)[z])
    roots = find_all_roots(locator_poly, alpha, reducing)

    alpha_powers = find_all_powers(alpha, reducing)
    result = []
    for root in roots:
        inverse = find_inverse(root, reducing, z) % reducing
        inverse_coefficients = inverse.all_coeffs()
        result.append(alpha_powers[tuple(inverse_coefficients)])
    return result


def decode_correct_code(encoded, generator):
    decoded, _ = div(encoded, generator)
    result = decoded.all_coeffs()
    result = [0]*(7-len(result))+result
    return result


def find_syndromes(encoded, alpha, reducing, t):
    syndromes = []
    for i in range(1,2*t+1):
        syndrome = substitute(encoded, alpha**i)
        syndrome %= reducing
        syndromes.append(syndrome)
    return syndromes


def decode_bch(bits, generator, alpha, reducing, t):
    encoded = Poly(bits, x, domain=GF(2))

    syndromes = find_syndromes(encoded, alpha, reducing, t)
    if all((syndrome == 0 for syndrome in syndromes)):
        return decode_correct_code(encoded, generator)

    locator = find_error_locator(syndromes, generator, reducing, t)
    errors = find_error_pos(locator, alpha, generator, reducing)

    for error in errors:
        encoded += Poly(x**error, x, domain=GF(2))
    return decode_correct_code(encoded, generator)


def main():
    m = 3
    n = 2**m - 1
    t = 2

    reducing = Poly(z**4 + z + 1, domain=GF(2))
    alpha = Poly(z, z, domain=GF(2))
    generator = find_generator(alpha, reducing, t)

    correct = [random.choice((0,1)) for _ in range(7)]
    print("Message:", correct)
    encoded = encode_bch(correct, generator)

    # Test all 1-bit errors
    for i in range(len(encoded)):
            error = encoded[:i]+[1-encoded[i]]+encoded[i+1:]
            corrected = decode_bch(error, generator, alpha, reducing, t)
            if corrected != correct:
                print("Mistake:", corrected, "!=", correct)
                return
    print("All 1-bit errors corrected!")

    # Test all 2-bit errors
    for i in range(len(encoded)-1):
        for j in range(i+1,len(encoded)):
            error = encoded[:i]+[1-encoded[i]]+encoded[i+1:j]+[1-encoded[j]]+encoded[j+1:]
            corrected = decode_bch(error, generator, alpha, reducing, t)
            if corrected != correct:
                print("Mistake:", corrected, "!=", correct)
                return
    print("All 2-bit errors corrected!")


if __name__ == "__main__":
    main()
