#!/usr/bin/env python

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

    def __init__ (self, coefficients):
        try:
            self.coefficients = [float (coefficients)]
        except TypeError:
            self.coefficients = coefficients
            self.variable = "x"

    def degree (self):
        return len (self.coefficients) - 1

    def __radd__ (self, other):
        return Expansion (other) + self

    def __add__ (self, other):
        length = max (len (self.coefficients), len (other.coefficients))
        coeff1 = self.coefficients + (length - len (self.coefficients)) * [0,]
        coeff2 = other.coefficients + (length - len (other.coefficients)) * [0,]
        coefficients = map (lambda x: x[0]+x[1], zip (coeff1, coeff2))
        return Expansion (coefficients)

    def __rmul__ (self, other):
        return Expansion (map (lambda x: other*x, self.coefficients))

    def __mul__ (self, other):
        coefficients = (self.degree () + other.degree () + 1) * [0.,]
        deg_x = 0
        for x in self.coefficients:
            deg_y = 0
            for y in other.coefficients:
                coefficients [deg_x + deg_y] += x*y
                deg_y += 1
            deg_x += 1
        return Expansion (coefficients)

    def __pow__ (self, i):
        x = int (i)
        if x == 0:
            return Expansion ([1.,])
        else:
            return self*(self**(i-1))

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
                monomial = self.variable if i==1 else self.variable + "^" + \
                    str (i)
                res += " + " + str (coeff) + monomial
        if res == "": res = "0"
        return res

class Polynomial2 (object):
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

    def __mul__ (self, other):
        result = Polynomial2 ()
        for k1, v1 in self.coefficients.items ():
            for k2,v2 in other.coefficients.items ():
                result += Polynomial2 (v1*v2, k1[0]+k2[0],k1[1]+k2[1])
        return result.clean ()

    def __pow__ (self, i):
        if not isinstance (i, int):
            raise TypeError ("power operator takes an int as first argument")
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

cos = reduce (lambda x,y : x+y, [((-1)**i/float(fact(2*i)))*(x**(2*i)) for i in range (10)])
sin = reduce (lambda x,y : x+y, [((-1)**i/float(fact(2*i+1)))*(x**(2*i+1)) for i in range (10)])

p1 = Polynomial2 (1, 0, 0)
p2 = Polynomial2 (1, 1, 0)
p3 = Polynomial2 (1, 0, 1)

alpha_r = Polynomial2 (1,1,1)
alpha_r.x0 = "alpha"
alpha_r.x1 = "r"
