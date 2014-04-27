csc6870-river
=============

CSC 6870 River simulation using palabos

Team Memebers:
Ryan Price
Brock Jackman
Chad Fisher

    Objective:
    Simulate a river flow through certain obstacles.  These obstacles include a small island,
    and a bridge support.  The river will start with a certain amount of predefined water, and the water
    will be allowed to flow through the river (i.e. the river.stl mesh file).


It's Working!

The mesh needs to have a reservior at its +x+y corner. The reservoir edge
needs to be at least 25% high and and enclose a 10% x 10% area in the mesh's
corner. (Sizes relative to total size of the mesh). The end of the flow should be
at the -x-y corner of the mesh, and the drain surface must reach within 12% of the
lowest point in the mesh.

stream.stl satisfies these conditions

The inlet pressure is a bit wonky, but values should be slightly above unity.
The pressure does not relate to physical units, and the apparent inflow speed
will depend on the size of the lattice so you'll have to play with it.
