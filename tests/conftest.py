from __future__ import absolute_import, print_function
import pytest
import moment_polytopes.disk_cache


def pytest_addoption(parser):
    # add --disable-disk-cache option
    parser.addoption(
        '--disable-disk-cache', action='store_true', help='disable disk cache')


def pytest_configure(config):
    if config.getoption('--disable-disk-cache'):
        print('Disk cache disabled')
        moment_polytopes.disk_cache.DISABLED = True
