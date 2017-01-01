Installation
============

This package requires SageMath_ 7.5 or higher. It can be installed as follows:

.. code-block:: bash

  sage -i lrslib
  sage -pip install git+git://github.com/catch22/moment_polytopes

SSL Troubleshooting
~~~~~~~~~~~~~~~~~~~

If on the second line you get an error message of the type ``pip is configured with locations that require TLS/SSL, however the ssl module in Python is not available``, please run the following and retry:

.. code-block:: bash

  sage -i openssl
  sage -f python2
