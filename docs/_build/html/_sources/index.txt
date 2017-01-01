Welcome!
========

This is a `SageMath <http://www.sagemath.org>`_ package for computing moment polytopes associated with finite-dimensional representations of compact and connected Lie groups, based on the algorithm proposed by `Vergne and Walter (2014) <https://arxiv.org/abs/1410.8144>`_.


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

We have used this package to compute moment cones associated with many other quantum marginal problems (TODO: add link).

You can install this package as follows:

.. code-block:: bash

  sage -i lrslib
  sage -pip install git+git://github.com/catch22/moment_polytopes

See :doc:`installation` for further detail and troubleshooting.

Citation
--------

If you use this software, please consider citing our article:

.. code-block:: bibtex

  @article{moment_polytopes,
    author  = {Vergne, M. and Walter, M.},
    title   = {Inequalities for moment cones of finite-dimensional representations},
    year    = {2014},
    journal = {to appear in Journal of Symplectic Geometry},
    eprint  = {1410.8144},
    note    = {Software available at \url{https://github.com/catch22/moment_polytopes}.},
  }

Documentation
-------------

.. toctree::
   :maxdepth: 1

   installation

.. toctree::
   :maxdepth: 2

   reference
   changes
   license
