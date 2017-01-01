Reference
=========

.. module:: moment_polytopes

Lie Group Representations
-------------------------

.. autoclass:: Representation
   :members:
.. autofunction:: weyl_module
.. autofunction:: external_tensor_product
.. autofunction:: positive_weyl_chamber_hrepr
.. autofunction:: is_dual_root_primitive
.. autofunction:: dual_root_primitive

Ressayre Elements
-----------------

.. autodata:: DEFAULT_FAILURE_PROBABILITY
.. autoclass:: RessayreTester
   :members:
.. autofunction:: ressayre_tester
.. autofunction:: is_ressayre

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

Internals
---------

.. autodata:: __version__

.. currentmodule:: moment_polytopes.disk_cache

.. autodata:: DISK_CACHE_DIR
.. autofunction:: disk_cache
