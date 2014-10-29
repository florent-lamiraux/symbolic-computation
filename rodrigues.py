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
from math import sqrt
from fractions import Fraction
from polynomial import Expansion, Polynomial, isscalar
import numpy as np
import matplotlib.pyplot as pl

class RodriguesLike (object):
    """
    Represent an expression of 3x3 matrix under the form
                                              2
    I + alpha (||r||) [r]  + alpha (||r||) [r]
     3       1           x                    x

    where alpha1 and alpha2 are polynomials, r = (r , r , r ) is a 3 vector
                                                   0   1   2
    and [r]  is the antisymmetric matrix defined by
           x

           [  0  -r   r  ]
           [       3   2 ]
     [r] = [  r   0  -r  ]
           [   3       1 ]
           [ -r   r   0  ]
           [   2   1     ]
    """
    alpha1 = 0
    alpha2 = 0
    #    3          2
    # [r]  = - ||r||  [r]
    #    x               x
    r2 = Expansion ([0,0,1])
    def __init__ (self, *args) :
        if len (args) == 2:
            self.alpha1 = args [0]
            self.alpha2 = args [1]
        elif len (args) == 1:
            self.alpha1 = args [0].alpha1
            self.alpha2 = args [0].alpha2

    def __mul__ (self, other) :
        result = RodriguesLike ()
        s1 = self.alpha1
        s2 = self.alpha2
        o1 = other.alpha1
        o2 = other.alpha2
        result.alpha1 = o1 + s1 - self.r2*(s1*o2+s2*o1)
        result.alpha2 = o2 + s2 + s1*o1 - self.r2*(s2*o2)
        return result

    def __call__ (self, *args):
        if isscalar (self.alpha1):
            alpha1 = self.alpha1
        else:
            alpha1 = self.alpha1 (*args)
        if isscalar (self.alpha2):
            alpha2 = self.alpha2
        else:
            alpha2 = self.alpha2 (*args)
        return RodriguesLike (alpha1, alpha2)

    def __str__ (self):
        result = "I + (%s) [r]_x + (%s) [r]_x^2"%(self.alpha1, self.alpha2)
        return result

if __name__ == '__main__':
    from polynomial import cos, sin, one_over_one_plus, one_over_one_minus
    r = Expansion ([0,1])
    Jinv = RodriguesLike ((1-cos)/r**2, (r - sin)/r**3)
    beta2 = (2*(1-cos)-r*sin)/(2*r**2*(1-cos))
    J = RodriguesLike (Fraction (-1,2), beta2)

    def variable (self, i):
        if i == 0:
            return "alpha"
        if i == 1:
            return "r"

    Polynomial.variable = variable
    alpha = Polynomial (1,0,1)
    r = Polynomial (1,1,1)

    x = Expansion ([0,1])
    RodriguesLike.r2 = r**2
    Jinv = RodriguesLike (alpha*((1-cos)/x**2)(alpha*r),
                          alpha**2*((x-sin)/x**3)(alpha*r))
    J = J (r)
    H1 = Jinv*J
    E_alphaT = RodriguesLike (-alpha*((sin/x)(alpha*r)),
                              alpha**2*(((1-cos)/x**2) (alpha*r)))

    H2 = E_alphaT*Jinv*J

    # square norm of difference with linear combination of velocities
    p = (alpha*r*H1.alpha1)**2 + (alpha*r**2*H1.alpha2)**2
    x = []; y = []
    for i in range (101):
        theta = Fraction (i,100)
        x.append (float (theta))
        y.append (sqrt (p (theta, Fraction (3,2))))

    x = np.array (x)
    y = np.array (y)

    fig = pl.figure ()
    ax = fig.add_subplot ('111')
    ax.plot (x, y)
    pl.show ()

