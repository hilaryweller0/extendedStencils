Stencil now has two components:
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
    no transform (see below)
        untransformed: (owner neighbour rest)
        transformed  : (rest)
    transform (cyclic patch)
        untransformed: (owner rest)
        transformed  : (neighbour rest)


To check for no transform:

    const globalIndexAndTransform& globalTransforms =
        mesh.globalData().globalTransforms();


    label patchI = ...
    const polyPatch& pp = mesh.boundaryMesh()[patchI];

    if (pp.coupled())
    {
        const labelPair& transSign =
            globalTransforms.patchTransformSign()[patchI];

        if (transSign == globalTransforms.nullTransformIndex())
        {
            // No transform, neighbour data is untransformedElements[1]
        }
        else
        {
            // Transform, neighbour data is transformedElements[0]
        }
    }
    else
    {
        // Uncoupled, no neighbour.
    }
