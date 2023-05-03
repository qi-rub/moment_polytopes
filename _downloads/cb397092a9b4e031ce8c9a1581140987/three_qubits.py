import moment_polytopes
import logging

# enable logging
logging.basicConfig(level=logging.DEBUG)

# compute three-qubit moment polytope in H-representation
three_qubits = (2, 2, 2)
hrepr = moment_polytopes.qmp.hrepr(three_qubits)
print(f"{len(hrepr.ieqs)} facets")

# convert to V-representation
vrepr = hrepr.vrepr()
print(f"{len(vrepr.vertices)} vertices")
