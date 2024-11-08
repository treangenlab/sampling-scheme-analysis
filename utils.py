import math
import numpy as np
from typing import Iterable, Union, Any
from collections import Counter, defaultdict
import itertools
from functools import lru_cache
from itertools import chain, zip_longest
import pickle, random
import copy

from copy import deepcopy
import pickle as pck
from functools import cache
from sympy import totient, isprime, divisors
from sympy.ntheory.residue_ntheory import mobius
from fractions import Fraction


# This name is weird?
def de_bruijn_invBWT(n, k):
    """de Bruijn sequence for alphabet k
    and subsequences of length n.
    """
    # Two kinds of alphabet input: an integer expands
    # to a list of integers as the alphabet..
    if isinstance(k, int):
        alphabet = list(map(str, range(k)))
    else:
        # While any sort of list becomes used as it is
        alphabet = k
        k = len(k)

    a = [0] * k * n
    sequence = []

    def db(t, p):
        if t > n:
            if n % p == 0:
                sequence.extend(a[1 : p + 1])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, k):
                a[t] = j
                db(t + 1, t)

    db(1, 1)
    return "".join(alphabet[i] for i in sequence)


def necklace(n, sigma=2):
    return int((sum(totient(d) * (sigma ** (n // d)) for d in divisors(n))) // n)


def aperiodic_necklaces(n, sigma=2):
    return int((sum(mobius(d) * sigma ** (n // d) for d in divisors(n))) // n)


def aperiodic_bound(w, k, sigma=2, sketch_size=1):
    if (w + k) * np.log2(sigma) > 31:
        return ragnar_ceil_LB(w, k)

    sampled = 0
    for d in divisors(w + k):
        primitive_n = aperiodic_necklaces(d, sigma)
        slack = math.ceil(sketch_size * d / w)
        sampled += primitive_n * slack
    return Fraction(sampled, (sigma ** (w + k)))


def aperiodic_bound_suff(w, k, sigma=2, sketch_size=1):
    ret = aperiodic_bound(w, k, sigma, sketch_size)
    if k % w != 1:
        ret = max(ret, aperiodic_bound(w, math.ceil(k / w) * w + 1, sigma, sketch_size))
    return ret


def ragnar_ceil_LB(w, k):
    return Fraction(math.ceil((w + k) / w), (w + k))


def ragnar_WABI_LB(w, k):
    return 1.5 / (w + k - 0.5)


def min_under_rot(s):
    mv = s
    for i in range(len(s)):
        s = s[-1] + s[:-1]
        if s < mv:
            mv = s
    return mv


def get_rots(s):
    rots = list()
    mv = s
    while len(rots) == 0 or s != mv:
        rots.append(s)
        s = s[1:] + s[0]
    return rots


@cache
def get_necklaces(v, sigma=2):
    nodes = (
        "".join(x)
        for x in itertools.product(list(str(x) for x in range(sigma)), repeat=v)
    )
    lwords = set()
    for vmer in nodes:
        lwords.add(min(vmer := (vmer[-1] + vmer[:-1]) for _ in range(len(vmer))))

    necklaces = defaultdict(list)
    for lw in lwords:
        word = lw
        for i in range(len(lw)):
            necklaces[lw].append(word)
            word = word[1:] + word[0]
            # word = word[-1] + word[:-1]
            if word == lw:
                break
    return necklaces
