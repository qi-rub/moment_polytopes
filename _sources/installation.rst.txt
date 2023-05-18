Installation
============

This package requires `SageMath <http://www.sagemath.org>`_ 9.8 or higher. It can be installed as follows:

.. code-block:: bash

  sage -pip install git+git://github.com/qi-rub/moment_polytopes

To test your installation, download the :download:`three_qubits.py <../examples/three_qubits.py>` example and run it as follows:

.. code-block:: bash

  sage three_qubits.py


Troubleshooting
---------------

If on the second line you get an error message saying that "pip is configured with locations that require TLS/SSL, however the ssl module in Python is not available", please run the following and retry:

.. code-block:: bash

  sage -i openssl
  sage -f python2


Mathematica Integration
-----------------------

`Wolfram Mathematica <https://www.wolfram.com/mathematica/>`_ contains some clever heuristics for evaluating determinants of polynomial matrices, and we provide the ``mathematica`` algorithm in :func:`moment_polytopes.ressayre_tester` etc. to leverage its functionality.

To use it, Mathematica needs to be installed and the ``math`` executable has to be available in the current ``PATH``. Run ``print(mathematica._install_hints())`` at the ``sage`` prompt for further information on how to set up SageMath's Mathematica integration.
