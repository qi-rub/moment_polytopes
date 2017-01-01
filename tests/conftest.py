from __future__ import absolute_import, print_function
import pytest
import moment_polytopes.disk_cache


def pytest_addoption(parser):
    # --disable-disk-cache disables the disk cache
    parser.addoption(
        '--disable-disk-cache', action='store_true', help='disable disk cache')

    # --algorithm controls algorithm used to test Ressayre elements
    parser.addoption(
        "--algorithms",
        "--algorithm",
        help="select algorithm(s) used to test Ressayre elements (separate by commas)"
    )


def pytest_configure(config):
    if config.getoption('--disable-disk-cache'):
        print('Disk cache disabled')
        moment_polytopes.disk_cache.DISABLED = True


def pytest_generate_tests(metafunc):
    if 'algorithm' in metafunc.fixturenames:
        algorithms = metafunc.config.getoption('algorithms')

        # no algorithms specified? use default
        if algorithms:
            algorithms = algorithms.split(',')
            metafunc.parametrize('algorithm', algorithms, ids=algorithms)
        else:
            metafunc.parametrize('algorithm', [None], ids=['default'])
