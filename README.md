Extending OpenFOAM extended stencils for periodic boundary conditions

Centred stencil (X=cell, I=face)

  X X X X
  X X|X X
  X X X X

Upwinded version consists of two stencils:

X X X X
X X X|X
X X X X

and

    X X X X
    X|X X X
    X X X X
