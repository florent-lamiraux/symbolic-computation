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

class Polynomial2 (object):
    """
    Polynomial of two variables
    """
    x0 = "x0"
    x1 = "x1"

    def __init__ (self, *args):
        """ Initialization by a monomial of the form
        c x^i y^j with
        c = coeff,
        i = pow0,
        j = pow1
        """
        if len (args) == 3:
            coeff, pow0, pow1 = args
            self.coefficients = {(pow0,pow1):coeff}
        elif len (args) == 0:
            self.coefficients = dict ()
        elif len (args) == 1:
            other = args [0]
            self.coefficients = other.coefficients.copy ()

    def clean (self):
        for k, v in self.coefficients.items ():
            if v == 0:
                del self.coefficients [k]
        return self

    def __rmul__ (self, other):
        result = Polynomial2 (self)
        for k in result.coefficients.keys ():
            result.coefficients [k] *= other
        return result.clean ()

    def __radd__ (self, other):
        """
        Addition with a real number
        """
        return self + Polynomial2 (other, 0, 0)

    def __rsub__ (self, other):
        """
        Substraction with a real number
        """
        return Polynomial2 (other, 0, 0) - self

    def __add__ (self, other):
        result = Polynomial2 (self)
        for k, c in other.coefficients.items ():
            if k in result.coefficients.keys ():
                result.coefficients [k] += c
            else:
                result.coefficients [k] = c
        return result.clean ()

    def __sub__ (self, other):
        result = Polynomial2 (self)
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
        result = Polynomial2 ()
        result.coefficients = coefficients
        return result
        
    def __mul__ (self, other):
        result = Polynomial2 ()
        if isscalar (other):
            for k1, v1 in self.coefficients.items ():
                result.coefficients [k1] = v1*other
        else:
            for k1, v1 in self.coefficients.items ():
                for k2,v2 in other.coefficients.items ():
                    result += Polynomial2 (v1*v2, k1[0]+k2[0],k1[1]+k2[1])
            result.clean ()
        return result

    def __pow__ (self, i):
        if not isinstance (i, int):
            raise TypeError ("power operator takes an int as second argument")
        if i==1:
            return Polynomial2 (self)
        if i==0:
            return Polynomial2 (1, 0, 0)
        else:
            return self*(self**(i-1))

    def __call__ (self, x0, x1):
        result = 0
        for k,v in self.coefficients.items ():
            result += v*x0**k[0]*x1**k[1]
        return result

    def truncate (self, degree):
        for k in self.coefficients.keys ():
            if reduce (lambda x, y: x+y, k) > degree:
                del self.coefficients [k]

    @staticmethod
    def comp (k1, k2):
        if k1 [0] + k1 [1] < k2 [0] + k2 [1] : return -1
        elif k1 [0] + k1 [1] == k2 [0] + k2 [1] :
            if k1 [0] < k2 [0] : return -1
            elif k1 [0] > k2 [0] : return 1
            else : return 0
        else: return 1

    def __str__ (self):
        result = ""
        for k in sorted (self.coefficients.keys (), Polynomial2.comp):
            if self.coefficients [k] > 0:
                result += " + "
            if self.coefficients [k] != 1 or k[0] == k[1] == 0:
                result += str (self.coefficients [k])
            if k [0] != 0:
                result += self.x0
                if k [0] != 1:
                    result += "^" + str (k[0])
            if k [1] != 0:
                result += "." + self.x1
                if k [1] != 1:
                    result += "^" + str (k[1])
        return result

x = Expansion ([0,1])

cos = reduce (lambda x,y : x+y, [Fraction ((-1)**i,fact(2*i))*x**(2*i) \
                                     for i in range (10)])
sin = reduce (lambda x,y : x+y, [Fraction ((-1)**i, fact(2*i+1))*x**(2*i+1) \
                                     for i in range (10)])
# 1/(1+x)
one_over_one_plus = reduce (lambda x,y : x+y, [(-1)**i*x**i for i in range (20)])
one_over_one_minus = reduce (lambda x,y : x+y, [x**i for i in range (20)])
p1 = Polynomial2 (1, 0, 0)
p2 = Polynomial2 (1, 1, 0)
p3 = Polynomial2 (1, 0, 1)

alpha_r = Polynomial2 (1,1,1)
alpha_r.x0 = "alpha"
alpha_r.x1 = "r"
