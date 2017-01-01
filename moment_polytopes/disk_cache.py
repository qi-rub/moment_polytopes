from __future__ import absolute_import, print_function
import logging, os.path
import cPickle as pickle
from functools import wraps
from . import __name__ as package_name
from . import __version__ as package_version

__all__ = ["DISK_CACHE_DIR", "disk_cache"]

logger = logging.getLogger(__name__)

#:Cache directory used by the :func:`disk_cache` decorator.
DISK_CACHE_DIR = os.path.expanduser('~/.cache/%s-%s' %
                                    (package_name, package_version))
logger.debug('Using directory %s', DISK_CACHE_DIR)

#:Globally disable disk cache.
DISABLED = False


def disk_cache(f):
    """
    Decorator that wraps a function such that calls are only evaluated once and otherwise retrieve from an on-disk cache.
    """
    MISSING = object()  # dummy object representing missing values in the cache
    fname = "%s.%s" % (f.__module__, f.__name__)  # absolute function name

    @wraps(f)
    def wrapper(*args, **kwargs):
        if DISABLED:
            return f(*args, **kwargs)

        key = (args, tuple(sorted(kwargs.items())))

        # lazily load cache
        if wrapper.cache is None:
            logger.debug('Loading cache for %s() from %s' %
                         (fname, wrapper.cache_path))
            if os.path.exists(wrapper.cache_path):
                wrapper.cache = pickle.load(open(wrapper.cache_path, 'rb'))
            else:
                wrapper.cache = {}

        # look up result in cache
        result = wrapper.cache.get(key, MISSING)
        if result is not MISSING:
            return result

        # call function and update cache
        logger.debug('Cache miss for %s() with args %r and kwargs %r' %
                     (fname, args, kwargs))
        result = f(*args, **kwargs)
        wrapper.cache[key] = result

        # save cache to disk
        pickle.dump(wrapper.cache, open(wrapper.cache_path, 'wb'))

        return result

    # create cache dir if does not exist
    if not os.path.exists(DISK_CACHE_DIR):
        os.mkdir(DISK_CACHE_DIR)

    # configure paths
    wrapper.cache_path = os.path.join(DISK_CACHE_DIR, "%s.p" % fname)
    wrapper.cache = None

    return wrapper
