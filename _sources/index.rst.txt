.. title:: moment_polytopes

Welcome!
========

This is a `SageMath <http://www.sagemath.org>`_ package for computing moment polytopes associated with finite-dimensional representations of compact and connected Lie groups, based on the algorithm proposed by `Vergne and Walter (2014) <https://arxiv.org/abs/1410.8144>`_ (see also `Walter (2014) <https://arxiv.org/abs/1410.6820>`_ for further detail).


Introduction
------------

The following code solves the one-body quantum marginal problem for three qubits:

.. literalinclude:: ../examples/three_qubits.py
   :language: python

You can install the latest version of this package as follows:

.. code-block:: bash

  sage -i lrslib
  sage -pip install git+git://github.com/qi-rub/moment_polytopes --upgrade

Now download the :download:`three_qubits.py <../examples/three_qubits.py>` example and run it via ``sage three_qubits.py``.
See :doc:`installation` for further detail and troubleshooting.

We have used this package to compute moment polytopes associated with various :doc:`qmp`.


Documentation
-------------

.. toctree::
   :maxdepth: 2

   installation
   qmp
   reference
   changes
   license


Citation
--------

If you use this software, please consider citing our article:

.. code-block:: bibtex

  @article{moment_polytopes,
    author  = {Vergne, M. and Walter, M.},
    title   = {Inequalities for moment cones of finite-dimensional representations},
    journal = {Journal of Symplectic Geometry},
    year    = {2017},
    volume  = {15},
    number  = {4},
    pages   = {1209--1250},
    eprint  = {1410.8144},
    doi     = {10.4310/JSG.2017.v15.n4.a8},
    note    = {Software available at \url{https://qi-rub.github.io/moment_polytopes/}.},
  }
