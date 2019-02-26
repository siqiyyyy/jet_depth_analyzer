from math import hypot, pi, sqrt, fabs
def deltaPhi(a,b):
    dphi = abs(float(a) - float(b))
    if (dphi <= pi): return dphi
    else: return 2*pi - dphi
        

def deltaR(a,b):
    dphi = deltaPhi(a,b)
    return hypot((a.eta()-b.eta()),dphi)

