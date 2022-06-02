/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "splitCyclicFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"
#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(splitCyclicFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, splitCyclicFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::splitCyclicFvPatch::makeWeights(fvsPatchScalarField& w) const
{
    const splitCyclicFvPatch& nbrPatch = neighbFvPatch();

    const scalarField deltas(nf()&fvPatch::delta());
    const scalarField nbrDeltas(nbrPatch.nf()&nbrPatch.fvPatch::delta());

    forAll(deltas, facei)
    {
        scalar di = deltas[facei];
        scalar dni = nbrDeltas[facei];

        w[facei] = dni/(di + dni);
    }
}


// Make patch face - neighbour cell distances
void Foam::splitCyclicFvPatch::makeDeltaCoeffs(fvsPatchScalarField& dc) const
{
    vectorField d = delta();
    vectorField n = nf();

    forAll(d, facei)
    {
        // Stabilised form for bad meshes.  HJ, 24/Aug/2011
        dc[facei] = 1.0/max(n[facei] & d[facei], 0.05*mag(d[facei]));
    }
}


Foam::tmp<Foam::vectorField> Foam::splitCyclicFvPatch::delta() const
{
    const vectorField patchD(fvPatch::delta());
    const vectorField nbrPatchD(neighbFvPatch().fvPatch::delta());

    tmp<vectorField> tpdv(new vectorField(patchD.size()));
    vectorField& pdv = tpdv();

    // To the transformation if necessary
    if (parallel())
    {
        forAll(patchD, facei)
        {
            vector ddi = patchD[facei];
            vector dni = nbrPatchD[facei];

            pdv[facei] = ddi - dni;
        }
    }
    else
    {
        forAll(patchD, facei)
        {
            vector ddi = patchD[facei];
            vector dni = nbrPatchD[facei];

            pdv[facei] = ddi - transform(forwardT()[0], dni);
        }
    }

    return tpdv;
}


Foam::tmp<Foam::labelField> Foam::splitCyclicFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::splitCyclicFvPatch::transfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    notImplemented("cyclicFvPatch::transfer");
    tmp<labelField> tpnf(new labelField(this->size()));
    return tpnf;
}


Foam::tmp<Foam::labelField> Foam::splitCyclicFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return neighbFvPatch().patchInternalField(iF);
}


// ************************************************************************* //
