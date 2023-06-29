/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "splitCyclicAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(splitCyclicAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        AMGInterfaceField,
        splitCyclicAMGInterfaceField,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::splitCyclicAMGInterfaceField::splitCyclicAMGInterfaceField
(
    const AMGInterface& AMGCp,
    const lduInterfaceField& fineInterface
)
:
    AMGInterfaceField(AMGCp, fineInterface),
    cyclicInterface_(refCast<const splitCyclicAMGInterface>(AMGCp)),
    doTransform_(false),
    rank_(0)
{
    const splitCyclicLduInterfaceField& p =
        refCast<const splitCyclicLduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::splitCyclicAMGInterfaceField::~splitCyclicAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::splitCyclicAMGInterfaceField::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{
    Info << "In splitCyclicAMGInterfaceField::updateInterfaceMatrix" << endl;
    Info << cyclicInterface_.neighbPatchID() << endl;
    Info << cyclicInterface_.neighbPatch().faceCells().size() << endl;

    // Get neighbouring field
    scalarField pnf
    (
        cyclicInterface_.neighbPatch().interfaceInternalField(psiInternal)
    );

    Info << pnf.size() << " " << coeffs.size() << " " << cyclicInterface_.faceCells().size() << endl;

    Info << "In splitCyclicAMGInterfaceField::updateInterfaceMatrix 1" << endl;

    const unallocLabelList& faceCells = cyclicInterface_.faceCells();

    transformCoupleField(pnf, cmpt);

    if (switchToLhs)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}


// ************************************************************************* //
