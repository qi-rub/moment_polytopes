Reference
=========

.. module:: moment_polytopes

Ressayre Elements
-----------------

This module contains functionality for computing moment polytopes of arbitrary finite-dimensional representations of compact, connected Lie groups (see :func:`is_ressayre` and :class:`Representation` below).

.. autodata:: DEFAULT_FAILURE_PROBABILITY
.. autoclass:: RessayreTester
   :members:
.. autofunction:: ressayre_tester
.. autofunction:: is_ressayre
.. autofunction:: is_admissible
.. autofunction:: c_candidates

Representations
---------------

.. autoclass:: Representation
   :members:
.. autofunction:: weyl_module
.. autofunction:: external_tensor_product
.. autofunction:: positive_weyl_chamber_hrepr
.. autofunction:: is_dual_root_primitive
.. autofunction:: dual_root_primitive

Quantum Marginal Problem
------------------------

This module contains functionality for computing moment polytopes associated with the pure-state marginal problem, i.e., the :math:`\times_i GL(d_i)`-representation :math:`\bigotimes_i \mathbb C^{d_i}`.

.. currentmodule:: moment_polytopes.qmp
.. autofunction:: hrepr
.. autofunction:: vrepr
.. autofunction:: H_AB_dominant
.. autofunction:: H_ABC_dominant
.. autofunction:: H_dominant_admissible
.. autofunction:: H_candidates
.. autofunction:: H_ressayre
.. autofunction:: facet_normal_form
.. autofunction:: ieqs_wo_perms
.. autofunction:: vertices_wo_perms

.. currentmodule:: moment_polytopes

Polyhedra
---------

.. autoclass:: HRepr
   :members:
.. autoclass:: VRepr
   :members:
.. autoclass:: LrsError
   :members:

Combinatorics
-------------

.. autofunction:: rect_tableaux
.. autofunction:: cubicle
.. autofunction:: cubicle_tableaux
.. autofunction:: is_dominant
.. autofunction:: is_extremal_edge
.. autofunction:: is_extremal_edge_ieq
.. autofunction:: extremal_edges
.. autofunction:: perms_of_length
.. autofunction:: length_tuples
.. autofunction:: is_shuffle
.. autofunction:: is_antishuffle
.. autofunction:: shuffles
.. autofunction:: antishuffles
.. autofunction:: perm_action
.. autoclass:: StabilizerGroup
   :members:

Third-Party Data
----------------

This module contains inequalities of moment polytopes computed by others.
We use them to verify that our implementation is correct.

.. currentmodule:: moment_polytopes.third_party
.. autodata:: KLYACHKO_FERMI_SCENARIOS
.. autofunction:: klyachko_fermi_hrepr
.. autodata:: KLYACHKO_QMP_SCENARIOS
.. autodata:: KLYACHKO_GOOD_QMP_SCENARIOS
.. autofunction:: klyachko_qmp_hrepr

Internals
---------

.. currentmodule:: moment_polytopes
.. autodata:: __version__

.. currentmodule:: moment_polytopes.disk_cache
.. autodata:: DISK_CACHE_DIR
.. autofunction:: disk_cache

.. currentmodule:: moment_polytopes.utils
.. autofunction:: dim_affine_hull
