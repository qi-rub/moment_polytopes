Welcome!
========

This is a `SageMath <http://www.sagemath.org>`_ package for computing moment polytopes associated with finite-dimensional representations of compact and connected Lie groups.

Introduction
------------

The following code solves the one-body quantum marginal problem for three qubits::

  import moment_polytopes

  # compute three-qubit moment cone in H-representation
  three_qubits = (2, 2, 2)
  hrepr = moment_polytopes.qmp_hrepr_irred(three_qubits)
  print '%s facets' % len(hrepr.ieqs)

  # convert to V-representation
  vrepr = hrepr.vrepr()
  print '%s vertices' % len(vrepr.vertices)

We have used this package to compute moment cones associated with many other quantum marginal problems.


Installation
------------

This package requires SageMath_ 7.5 or higher. It can be installed as follows:

.. code-block:: bash

  sage -i lrslib
  sage -pip install git+git://github.com/catch22/moment_polytopes

SSL Troubleshooting
~~~~~~~~~~~~~~~~~~~

If on the second line you get an error message of the type ``pip is configured with locations that require TLS/SSL, however the ssl module in Python is not available.``, please run the following and retry:

.. code-block:: bash

  sage -i openssl
  sage -f python2

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Miscellaneous Pages
-------------------

.. toctree::
   :maxdepth: 2

   changes
   license
