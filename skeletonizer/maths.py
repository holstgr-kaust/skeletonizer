"""
    Skeletonizer: Python Cell Morphology Analysis and Construction Toolkit

    KAUST, BESE, Neuro-Inspired Computing Project
    (c) 2014-2015. All rights reserved.
"""
"""
    Skeletonize maths module.
"""

import math
import operator


def vlogger(func):
    def inner(*args, **kwargs):
        print("Called: %s(%s, %s)" % (func.__name__, args, kwargs))
        result = func(*args, **kwargs)
        print(" result: %s, len:%s\n" % (str(result), str(vlength(result))))
        return result
    return inner

def square(x):
    return x * x

def distance_squared(v1, v2):
    return sum(map(lambda x, y: square(x - y), v1, v2))

def distance_squared(v1, v2):
    return sum(map(lambda x, y: square(x - y), v1, v2))

def distance(v1, v2):
    return math.sqrt(distance_squared(v1, v2))


def vlength(vect):
    return math.sqrt(sum(map(lambda v: square(v), vect)))

def vmuls3(v, x):
    return (v[0]*x, v[1]*x, v[2]*x)

def vdivs3(v, x):
    return (v[0]/x, v[1]/x, v[2]/x)

def vadds3(v, x):
    return (v[0]+x, v[1]+x, v[2]+x)

def vsubs3(v, x):
    return (v[0]-x, v[1]-x, v[2]-x)

def vadd3(v1, v2):
    return (v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2])

def vsub3(v1, v2):
    return (v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])

def vmin3(v1, v2):
    return (min(v1[0],v2[0]), min(v1[1],v2[1]), min(v1[2],v2[2]))

def vmax3(v1, v2):
    return (max(v1[0],v2[0]), max(v1[1],v2[1]), max(v1[2],v2[2]))


def vnormalize3(v):
    m = vlength(v)
    assert(m != 0)
    return vdivs3(v, m)

def vnormalize_zero3(v):
    m = vlength(v)
    return vnormalize3(v) if m != 0 else v


def v3_to_aabb(v1, v2):
    """
    Creates an AABB (Axis Aligned Bounding Box) from two maximally extreme points of the box.
    :param v1: pos of first corner.
    :param v2: pos of second corner, maximally distant to v1.
    :return: pair of (min, max) vectors describing AABB.
    """
    return (vmin3(v1, v2), vmax3(v1, v2))

def adjust_aabb(aabb, n):
    """
    Creates an adjusted AABB where each side is moved out from, or closer to, the centre.
    :param aabb: The source AABB (Axis Aligned Bounding Box).
    :param n: scaler of the amount to grow / reduce each side (positive values grow; negative values shrink AABB).
    :return: Adjusted AABB.
    """
    return v3_to_aabb(vsubs3(aabb[0], n), vadds3(aabb[1], n))

def inside_aabb(aabb, v):
    """
    Tests if a point is inside the AABB.
    :param aabb: The AABB (Axis Aligned Bounding Box).
    :param v: position vector.
    :return: True if v is inside, False otherwise.
    """
    inside_min = all(map(lambda (aabbn, vn): aabbn < vn, zip(aabb[0], v)))
    inside_max = all(map(lambda (aabbn, vn): aabbn > vn, zip(aabb[1], v)))
    return inside_min and inside_max


#@vlogger
def vadjust_offset_length3(v, centre, min_length):
    """
    Returns a vector offset as if centre point was moved to zero; and, at least as long as min_length.
    """
    nv = vsub3(v, centre)
    m = vlength(nv)
    return nv if m > min_length else vmuls3(vnormalize_zero3(nv), min_length)



