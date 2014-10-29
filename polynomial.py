#!/usr/bin/env python
#
# Copyright (c) CNRS/LAAS 2014
# Author: Florent Lamiraux
#
# This file is part of symbolic-computation.
# symbolic-computation is free software: you can redistribute it
# and/or modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# symbolic-computation is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Lesser Public License for more details.  You should have
# received a copy of the GNU Lesser General Public License along with
# symbolic-computation.  If not, see
# <http://www.gnu.org/licenses/>.

from __future__ import division
from fractions import Fraction
from collections import Iterable

def isscalar (x):
    if isinstance (x, Fraction): return True
    if isinstance (x, int): return True
    return False

factorial = [1, 1,]
def fact (n):
    if n < len (factorial):
        return factorial [n]
    else:
        prev = fact (n-1)
        res = prev * n
        factorial.append (res)
        return res

combination = {(0,0):1}
def C (i,j):
    x = int (i); y = int (j)
    if x == 0:
        return 1
    if (x, int(j)) in combination.keys ():
        return combination [(i,j)]
    else:
        combination [(x, y)] = C (x,y)
        return combination [(x, y)]


class Expansion (object):
    """
    Polynomial of one variable
    """
    x = "x"
    maxDegree = 20

    def __init__ (self, *args):
        if len (args) == 1:
            if isinstance(args [0], Iterable):
                self.coefficients = args [0]
            else:
                self.coefficients = [Fraction (args [0])]
        else:
            raise TypeError ("Expecting a Fraction or a Expansion")

    def truncate (self):
        self.coefficients = self.coefficients [:self.maxDegree]
        return self

    def degree (self):
        return len (self.coefficients) - 1

    def __radd__ (self, other):
        return Expansion (other) + self

    def __rsub__ (self, other):
        return Expansion (other) - self

    def __sub__ (self, other):
        try:
            length = max (len (self.coefficients), len (other.coefficients))
            coeff2 = other.coefficients + \
                (length - len (other.coefficients)) * [0,]
        except AttributeError:
            length = max (len (self.coefficients), 1)
            coeff2 = [other,] + (length - 1) * [0,]

        coeff1 = self.coefficients + (length - len (self.coefficients)) * [0,]
        coefficients = map (lambda x: x[0]-x[1], zip (coeff1, coeff2))
        return Expansion (coefficients)

    def __neg__ (self):
        return Expansion (map (lambda x: -x, self.coefficients))

    def __add__ (self, other):
        try:
            length = max (len (self.coefficients), len (other.coefficients))
            coeff2 = other.coefficients + \
                (length - len (other.coefficients)) * [0,]
        except AttributeError:
            length = max (len (self.coefficients), 1)
            coeff2 = [other,] + (length - 1) * [0,]

        coeff1 = self.coefficients + (length - len (self.coefficients)) * [0,]
        coefficients = map (lambda x: x[0]+x[1], zip (coeff1, coeff2))
        return Expansion (coefficients)

    def __rmul__ (self, other):
        return Expansion (map (lambda x: other*x, self.coefficients))

    def __mul__ (self, other):
        try:
            coefficients = (self.degree () + other.degree () + 1) * [0,]
            deg_x = 0
            for x in self.coefficients:
                deg_y = 0
                for y in other.coefficients:
                    coefficients [deg_x + deg_y] += x*y
                    deg_y += 1
                deg_x += 1
        except AttributeError:
            coefficients = map (lambda x : x*other, self.coefficients)
        return Expansion (coefficients).truncate ()

    def __pow__ (self, i):
        x = int (i)
        if x == 0:
            return Expansion ([1,])
        else:
            return (self*(self**(i-1))).truncate ()

    def __truediv__ (self, other):
        result = Expansion (self.coefficients)
        if isscalar (other):
            result.coefficients = map (lambda x : x/other, self.coefficients)
        else:
            numerator = self.coefficients
            denominator = other.coefficients
            # Concurrently decrease degree of numerator and denominator
            while numerator [0] == 0 and denominator [0] == 0:
                del numerator [0]; del denominator [0]
            if denominator [0] == 0:
                raise ZeroDivisionError \
                    ("degree of denominator is higher than degree of numerator")
            if len (denominator) == 1:
                result /= denominator [0]
            else:
                # Make coefficient of degree zero of denominator equal to 1
                coeff = denominator [0]
                numerator = map (lambda x: x/coeff, numerator)
                denominator = [0,] + map (lambda x: x/coeff, denominator [1:])
                result = Expansion (numerator)*one_over_one_plus \
                    (Expansion (denominator))
        return result

    def __call__ (self, other):
        res = self.coefficients [0]
        for i in range (1, self.degree () + 1):
            res += self.coefficients [i] * other**i
        return res

    def __str__ (self):
        res = str (self.coefficients [0]) if self.coefficients [0] != 0 else ""
        for i in range (1, self.degree () + 1):
            if self.coefficients [i] != 0:
                coeff = self.coefficients [i] if self.coefficients [i] != 1 \
                    else ""
                monomial = self.x if i==1 else self.x + "^" + \
                    str (i)
                res += " + " + str (coeff) + monomial
        if res == "": res = "0"
        return res

class Polynomial (object):
    """
    Polynomial of n variables
    """
    @staticmethod
    def variable (i):
        return "x_{%i}"%i

    def __init__ (self, *args):
        """ Initialization by a monomial of the form
        c x_i^j with args = (c, i, j)
        variable indexing starts at 0
        """
        if len (args) == 3:
            c, i, j = args
            if j == 0:
                self.coefficients = {() : c}
            else:
                self.coefficients = {i*(0,)+(j,):c}
        elif len (args) == 0:
            self.coefficients = dict ()
        elif len (args) == 1:
            other = args [0]
            self.coefficients = other.coefficients.copy ()
        elif len (args) == 2:
            # value, exponent
            value = args [0]; exp = tuple (args [1])
            self.coefficients = {exp : value}

    def clean (self):
        for k, v in self.coefficients.items ():
            if v == 0:
                del self.coefficients [k]
        return self

    def __rmul__ (self, other):
        result = Polynomial (self)
        for k in result.coefficients.keys ():
            result.coefficients [k] *= other
        return result.clean ()

    def __radd__ (self, other):
        """
        Addition with a real number
        """
        return self + Polynomial (other, 0, 0)

    def __rsub__ (self, other):
        """
        Substraction with a real number
        """
        return Polynomial (other, 0, 0) - self

    def __add__ (self, other):
        result = Polynomial (self)
        for k, c in other.coefficients.items ():
            if k in result.coefficients.keys ():
                result.coefficients [k] += c
            else:
                result.coefficients [k] = c
        return result.clean ()

    def __sub__ (self, other):
        result = Polynomial (self)
        for k, c in other.coefficients.items ():
            if k in result.coefficients.keys ():
                result.coefficients [k] -= c
            else:
                result.coefficients [k] = -c
        return result.clean ()

    def __neg__ (self):
        coefficients = dict ()
        for k in self.coefficients.keys ():
            coefficients [k] = -self.coefficients [k]
        result = Polynomial ()
        result.coefficients = coefficients
        return result

    def __mul__ (self, other):
        result = Polynomial ()
        if isscalar (other):
            for k1, v1 in self.coefficients.items ():
                result.coefficients [k1] = v1*other
        else:
            for k1, v1 in self.coefficients.items ():
                for k2,v2 in other.coefficients.items ():
                    l1 = len (k1); l2 = len (k2)
                    l = max (l1, l2)
                    k1 = k1 + (l - l1)* (0,)
                    k2 = k2 + (l - l2)* (0,)
                    result += Polynomial (v1*v2, map (lambda x: x [0] +x [1],
                                                      zip (k1, k2)))
            result.clean ()
        return result

    def __pow__ (self, i):
        if not isinstance (i, int):
            raise TypeError ("power operator takes an int as second argument")
        if i==1:
            return Polynomial (self)
        if i==0:
            return Polynomial (1, 0, 0)
        else:
            return self*(self**(i-1))

    def __call__ (self, *args):
        result = 0
        for k,v in self.coefficients.items ():
            monomial = Polynomial (v, 0, 0)
            for ki, xi in zip (k, args):
                monomial *= xi**ki
            result += monomial
        return result

    def truncate (self, degree):
        for k in self.coefficients.keys ():
            if reduce (lambda x, y: x+y, k) > degree:
                del self.coefficients [k]

    def __str__ (self):
        result = ""
        for k in sorted (self.coefficients.keys ()):
            if self.coefficients [k] > 0:
                result += " + "
            result += str (self.coefficients [k])

            for i, ki in zip (xrange (1000000), k):
                if ki != 0:
                    result += self.variable (i)
                    if ki != 1:
                        result += "^" + str (ki)
        return result

    def __int__ (self):
        if len (self.coefficients.keys ()) == 0:
            return 0
        if self.coefficients.keys () == [()]:
            return int (self.coefficients [()])
        raise TypeError ("Cannot convert polynomial to float")

    def __float__ (self):
        if len (self.coefficients.keys ()) == 0:
            return 0.0
        if self.coefficients.keys () == [()]:
            return float (self.coefficients [()])
        raise TypeError ("Cannot convert polynomial to float")



x = Expansion ([0,1])

cos = reduce (lambda x,y : x+y, [Fraction ((-1)**i,fact(2*i))*x**(2*i) \
                                     for i in range (10)])
sin = reduce (lambda x,y : x+y, [Fraction ((-1)**i, fact(2*i+1))*x**(2*i+1) \
                                     for i in range (10)])
# 1/(1+x)
one_over_one_plus = reduce (lambda x,y : x+y, [(-1)**i*x**i for i in range (20)])
one_over_one_minus = reduce (lambda x,y : x+y, [x**i for i in range (20)])
