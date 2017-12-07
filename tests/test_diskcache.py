from __future__ import absolute_import, print_function
from moment_polytopes.disk_cache import disk_cache


@disk_cache
def slow(x, y=1):
    slow.callcount += 1
    return x + 2 * y


slow.callcount = 0


def test_disk_cache():
    # initialize test
    slow.cache = {}
    assert slow.callcount == 0

    # check call counts (we are quite pessimistic with regards to the different ways of calling a function with the same parameters)
    assert slow(3) == 5
    assert slow.callcount == 1

    assert slow(3, 1) == 5
    assert slow.callcount == 2

    assert slow(3, 2) == 7
    assert slow.callcount == 3

    assert slow(3, y=2) == 7
    assert slow.callcount == 4

    assert slow(y=2, x=3) == 7
    assert slow.callcount == 5
