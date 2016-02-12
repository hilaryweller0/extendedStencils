 face stencil:
- now has two components:
    - untransformed part. Uses globalIndex numbering into cells and
      boundary faces.
    - transformed part. Uses labelPair from globalIndexAndTransform.

We need to know from the stencil which on is the owner data and which the
neighbour. Owner data never transformed, neighbour data might be.


If we know the patch the face is on:

- internal face:
    untransformed: (owner neighbour rest)

- boundary face:
    untransformed: (owner rest)

- coupled face
    no transform (processor patch)
        untransformed: (owner neighbour rest)
        transformed  : (rest)
    transform (cyclic patch)
        untransformed: (owner -1 rest)
        transformed  : (neighbour rest)


TODO:
- remove globalTransforms member data
- move nullIndex transformed data to untransformed
- test positional transformation (collectData?)
