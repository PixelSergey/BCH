"""BCH encoding and decoding project.

University of Oxford
Information Theory
MT24 mini-project
"""

import argparse
import random
import sys
import warnings
from sympy import GF, ZZ, Matrix, Poly, div, Symbol, Add, Mul, Pow, Integer, degree
from sympy.abc import x, z
from sympy.polys.galoistools import gf_gcdex, gf_lcm, gf_irreducible

warnings.filterwarnings("ignore", category=DeprecationWarning)


class BCH:
    """BCH code properties and functionality.

    On initialisation, constructs the required parameters for the BCH code.
    Encoding can be called using the `encode` function,
    and error-corrected decoding using the `decode` function.

    Attributes:
        m: The exponent of the Galois field size
        n: The codeword length
        t: The number of errors the code is able to correct
        c: The number of checksum bits in the code
        k: The number of data bits that can be encoded
        primitive: The primitive polynomial of the Galois field
        alpha: A primitive element of the Galois field
        generator: The generator polynomial used by the BCH code
    """

    def __init__(self, m, t, primitive=None):
        """Initialises the BCH code with the given parameters.

        Args:
            m: The exponent of the Galois field size. The codeword length is m**2-1
            t: The number of errors to correct 
            primitive: An integer representing a primitive polynomial (optional)       
        """
        self.m = m
        self.n = 2**m - 1
        self.t = t

        if not primitive:
            self.primitive = self.find_primitive()
        else:
            self.primitive = Poly([int(x) for x in primitive], z, domain=GF(2))
        self.alpha = Poly(z, z, domain=GF(2))
        self.generator = self.find_generator()

        self.c = degree(self.generator)
        self.k = self.n - self.c


    def find_primitive(self):
        """Find a primitive polynomial for the Galois field
        
        Returns:
            A primitive polynomial for the Galois field
        """
        while True:
            irreducible = Poly(gf_irreducible(self.m, 2, ZZ), z, domain=GF(2))
            for i in range(1, self.n):
                test_poly = Poly(z**i - 1, z, domain=GF(2))
                _, remainder = div(test_poly, irreducible)
                if remainder == 0:
                    break
            else:
                return irreducible


    def find_generator(self):
        """Find a generator polynomial for the Galois field
        
        Returns:
            A generator polynomial based on the primitive polynomial and alpha
        """
        generator = Poly(1, x, domain=GF(2))
        for i in range(1, 2*self.t):
            current = self.find_minimal_polynomial(self.alpha**i)
            generator = Poly(gf_lcm(generator.all_coeffs(), current.all_coeffs(), 2, ZZ), x, domain=GF(2))
        return generator


    def find_minimal_polynomial(self, element):
        """Find a minimal polynomial (mod the primitive polynomial) for a given element.
        
        Args:
            element: The element to find a minimal polynomial for. Usually a power of alpha.
        
        Returns:
            A minimal polynomial for the element
        """
        seen = set()
        i = 0
        result = Poly(1, x)
        while True:
            current = Poly(element**(2**i), z) % self.primitive
            if current in seen:
                break
            result *= Poly(x, x) - current.set_domain(ZZ)
            seen.add(current)
            i += 1

        result = Poly(result, x, domain=GF(2)[z])
        reduced = [self.expand_expression(coefficient)%self.primitive for coefficient in result.all_coeffs()]
        assert(all(coefficient in (0,1) for coefficient in reduced))

        result = Poly(reduced, x, domain=GF(2))
        return result


    def expand_expression(self, expression):
        """Expands an Expr type by evaluating all operations in the Galois field

        Args:
            expression: An Expr type to expand
        
        Returns:
            a Poly representing the evaluated expression
        """
        if isinstance(expression, Symbol):
            return Poly(expression, z, domain=GF(2)) % self.primitive

        if isinstance(expression, Integer):
            return Poly(expression, z, domain=GF(2)) % self.primitive

        if isinstance(expression, Mul):
            result = Poly(1, z, domain=GF(2))
            for term in expression.args:
                result *= self.expand_expression(term)
            return result % self.primitive

        if isinstance(expression, Add):
            result = Poly(0, z, domain=GF(2))
            for term in expression.args:
                result += self.expand_expression(term)
            return result % self.primitive

        if isinstance(expression, Pow):
            base, exponent = expression.args
            base = self.expand_expression(base)
            exponent = int(exponent)
            if exponent < 0:
                assert self.primitive is not None
                base = self.find_inverse(base)
                exponent = abs(exponent)

            return Poly(base**exponent, z, domain=GF(2)) % self.primitive


    def find_inverse(self, polynomial):
        """Finds the inverse of a polynomial modulo the primitive polynomial

        Args:
            polynomial: The polynomial to find the inverse for
        
        Returns:
            The inverse of the polynomial
        """
        inv, _, gcd = gf_gcdex(polynomial.all_coeffs(), self.primitive.all_coeffs(), 2, ZZ)
        assert gcd == [1]
        return Poly(inv, z, domain=GF(2))


    def substitute(self, polynomial, substitution):
        """Substitute a polynomial into the variables of another.

        Args:
            polynomial: The polynomial to substitute into. This will have its variables replaced.
            substitution: The polynomial to insert
        
        Returns:
            The evaluated expression `polynomial(substitution(z))`
        """
        result = Poly(0, z, domain=GF(2))
        for i, coeff in enumerate(polynomial.all_coeffs()[::-1]):
            result += Poly(coeff, z, domain=GF(2)) * Poly(substitution**i, z, domain=GF(2))
        return result


    def find_all_powers(self, element):
        """Find all powers of an element in the Galois field.
        
        This is useful for looking up powers based on an expression later on.

        Args:
            element: The element to find powers for, usually alpha

        Returns:
            A dict containing the coefficients of polynomials as keys
            and exponents of `element` as values
        """
        result = {}
        for i in range(0, self.n):
            power = Poly(element**i, z, domain=GF(2)) % self.primitive
            if tuple(power.all_coeffs()) in result:
                continue
            result[tuple(power.all_coeffs())] = i
        return result


    def find_all_roots(self, polynomial):
        """Find all the roots of a polynomial in the Galois field.
        
        Args:
            polynomial: The desired polynomial
        
        Returns:
            A list of roots in the Galois field.
            These are composed of powers of `alpha` reduced modulo the primitive polynomial.
        """
        roots = []
        for i in range(1,self.n+1):
            root = self.substitute(polynomial, self.alpha**i) % self.primitive
            if root == 0:
                roots.append(self.alpha**i % self.primitive)
        return roots


    def encode(self, bits):
        """Encode a message using the generated BCH code.

        Args:
            bits: A list containing the bits to encode
        
        Returns:
            A list containing the bits of the codeword
        """
        normalised = self.fill_data(bits, self.k)

        data = Poly(normalised, x, domain=GF(2))
        encoded = data * self.generator
        return self.fill_data(encoded.all_coeffs(), self.n)


    def decode(self, bits):
        """Decode a codeword using the generated BCH code.

        Args:
            bits: A list containing the bits of the codeword
        
        Returns:
            A list containing the bits of the decoded message, corrected for errors
        """
        normalised = self.fill_data(bits, self.n)
        encoded = Poly(normalised, x, domain=GF(2))

        syndromes = self.find_syndromes(encoded)
        if all((syndrome == 0 for syndrome in syndromes)):
            return self.decode_correct_code(encoded)

        locator = self.find_error_locator(syndromes)
        errors = self.find_error_pos(locator)

        for error in errors:
            encoded += Poly(x**error, x, domain=GF(2))

        return self.decode_correct_code(encoded)


    def find_syndromes(self, encoded):
        """Find the syndromes of a codeword

        Args:
            A polynomial representing any codeword
        
        Returns:
            A list of syndromes of the polynomial.
            If no errors have occurred, all syndromes are zero.
        """
        syndromes = []
        for i in range(1,2*self.t+1):
            syndrome = self.substitute(encoded, self.alpha**i)
            syndrome %= self.primitive
            syndromes.append(syndrome)
        return syndromes


    def find_error_locator(self, syndromes):
        """Find the error locator polynomial given syndromes.

        Args:
            syndromes: A list of syndromes obtained from `find_syndromes`

        Returns:
            An error locator vector,
            where the entries are coefficients of the error locator polynomial
        """
        for i in range(self.t):
            nu = self.t-i  # Number of errors
            syndrome_matrix = Matrix(nu, nu, lambda a,b: syndromes[a+b].as_expr())
            detection = Poly(syndrome_matrix.det(), z, domain=GF(2)) % self.primitive
            if detection == 0:
                continue

            syndrome_vector = Matrix(nu, 1, lambda a,_: -syndromes[nu+a].as_expr())
            augmented = syndrome_matrix.col_insert(nu, syndrome_vector)
            locator = augmented.rref(pivots=False).col(nu)
            result = []
            for row in range(nu):
                result.append(self.expand_expression(locator[row]))

            return result


    def find_error_pos(self, locator):
        """Find the error positions based on the error locator vector.
        
        Args:
            locator: The error locator vector obtained from `find_error_locator`
        
        Returns:
            A list of positions where errors have occurred in the codeword
        """
        locator_poly = Poly(1, x, domain=GF(2)[z])
        for i, lambda_i in enumerate(locator[::-1], start=1):
            locator_poly += Poly(lambda_i%self.primitive, x, domain=GF(2)[z]) * Poly(x**i, x, domain=GF(2)[z])
        roots = self.find_all_roots(locator_poly)

        alpha_powers = self.find_all_powers(self.alpha)
        result = []
        for root in roots:
            inverse = self.find_inverse(root) % self.primitive
            inverse_coefficients = inverse.all_coeffs()
            result.append(alpha_powers[tuple(inverse_coefficients)])
        return result


    def fill_data(self, data, length):
        """Prepend a list of bits with zeroes if the length is too short
        
        Args:
            data: The list of bits to complete with zeroes
            length: The desired length of the list

        Returns:
            A list containing the original data prepended with zeroes.
            If the length of the list is longer than the desired length, an error is thrown.
        """
        assert len(data) <= length
        return [0]*(length-len(data)) + data


    def decode_correct_code(self, encoded):
        """Given a codeword that is known to be correct, decode it.
        
        Args:
            encoded: A polynomial which has had its errors corrected

        Returns:
            A list containing the decoded codeword
        """
        decoded, _ = div(encoded, self.generator)
        result = decoded.all_coeffs()
        return self.fill_data(result, self.k)


def test(bch):
    """
    Code to test BCH functionality.
    Generates a BCH code, sends a random message, and tries to correct every possible error.
    """

    correct = [random.choice((0,1)) for _ in range(bch.k)]
    print("Message:", correct)
    encoded = bch.encode(correct)

    # Test all 1-bit errors
    for i, bit in enumerate(encoded):
        error = encoded[:i] + [1-bit] + encoded[i+1:]
        corrected = bch.decode(error)
        assert corrected == correct

    print("All 1-bit errors corrected!")

    # Test all 2-bit errors
    for i, bit1 in enumerate(encoded[:-1]):
        for j, bit2 in enumerate(encoded[i+1:], start=i+1):
            error = encoded[:i] + [1-bit1] + encoded[i+1:j] + [1-bit2] + encoded[j+1:]
            corrected = bch.decode(error)
            assert corrected == correct

    print("All 2-bit errors corrected!")


def main():
    parser = argparse.ArgumentParser(prog="BCH")
    parser.add_argument("-m", "--exponent", required=True)
    parser.add_argument("-t", "--errors-corrected", required=True)
    parser.add_argument("-p", "--primitive")
    data = parser.add_mutually_exclusive_group()
    data.add_argument("-e", "--encode")
    data.add_argument("-d", "--decode")
    data.add_argument("-x", "--test", action="store_true")
    args = parser.parse_args(sys.argv[1:])

    bch = BCH(int(args.exponent), int(args.errors_corrected), primitive=args.primitive)
    if args.test:
        test(bch)

    elif args.encode:
        message = args.encode.encode("ascii")
        bitstring = [int(a) for a in "".join([bin(character)[2:].zfill(7) for character in message])]
        padding_length = -len(bitstring) % bch.k
        bitstring += [0]*padding_length
        chunks = [bitstring[bch.k*i:bch.k*(i+1)] for i in range((len(bitstring) + bch.k - 1) // bch.k)]
        result = [bch.encode(chunk) for chunk in chunks]
        print("".join(["".join([str(bit) for bit in chunk]) for chunk in result]))
        if not args.primitive:
            print("primitive:", "".join(str(x) for x in bch.primitive.all_coeffs()))

    elif args.decode:
        message = args.decode
        bitstring = [int(a) for a in message]
        chunks = [bitstring[bch.n*i:bch.n*(i+1)] for i in range((len(bitstring) + bch.n - 1) // bch.n)]
        decoded = [bch.decode(chunk) for chunk in chunks]
        decoded_string = "".join(["".join([str(bit) for bit in chunk]) for chunk in decoded])
        padding_length = len(decoded_string) % 7
        decoded_string = decoded_string[:-padding_length] if padding_length>0 else decoded_string
        decoded_chunks = [decoded_string[7*i:7*(i+1)] for i in range((len(decoded_string) + 7 - 1) // 7)]
        result = bytes([int(chunk, 2) for chunk in decoded_chunks])
        print(result.decode("ascii"))


if __name__ == "__main__":
    main()
