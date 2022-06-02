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

#include "splitCyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "SubField.H"
#include "unitConversion.H"
#include "transform.H"
#include "matchPoints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(splitCyclicPolyPatch, 0);

    // int splitPolyPatch::disallowGenericsplitPolyPatch
    // (
    //     debug::debugSwitch("disallowGenericsplitPolyPatch", 0)
    // );

    //defineRunTimeSelectionTable(splitCyclicPolyPatch, word);
    //defineRunTimeSelectionTable(splitCyclicPolyPatch, dictionary);

    addToRunTimeSelectionTable(polyPatch, splitCyclicPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, splitCyclicPolyPatch, dictionary);

    template<>
    const char* NamedEnum<splitCyclicPolyPatch::transformType, 5>::names[] =
    {
        "unknown",
        "rotational",
        "translational",
        "coincidentFullMatch",
        "noOrdering"
    };

    const NamedEnum<splitCyclicPolyPatch::transformType, 5>
        splitCyclicPolyPatch::transformTypeNames;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::splitCyclicPolyPatch::findMaxArea
(
    const pointField& points,
    const faceList& faces
)
{
    label maxI = -1;
    scalar maxAreaSqr = -GREAT;

    forAll(faces, faceI)
    {
        scalar areaSqr = magSqr(faces[faceI].normal(points));

        if (areaSqr > maxAreaSqr)
        {
            maxAreaSqr = areaSqr;
            maxI = faceI;
        }
    }
    return maxI;
}


void Foam::splitCyclicPolyPatch::calcTransforms()
{
    if (size())
    {
        // Half0
        const splitCyclicPolyPatch& half0 = *this;
        vectorField half0Areas(half0.size());
        forAll(half0, facei)
        {
            half0Areas[facei] = half0[facei].normal(half0.points());
        }

        // Half1
        const splitCyclicPolyPatch& half1 = neighbPatch();
        vectorField half1Areas(half1.size());
        forAll(half1, facei)
        {
            half1Areas[facei] = half1[facei].normal(half1.points());
        }

        calcTransforms
        (
            half0,
            half0.faceCentres(),
            half0Areas,
            half1.faceCentres(),
            half1Areas
        );
    }
}


void Foam::splitCyclicPolyPatch::calcTransforms
(
    const primitivePatch& half0,
    const pointField& half0Ctrs,
    const vectorField& half0Areas,
    const pointField& half1Ctrs,
    const vectorField& half1Areas
)
{
    if (debug && master())
    {
        fileName casePath(boundaryMesh().mesh().time().path());
        {
            fileName nm0(casePath/name()+"_faces.obj");
            Pout<< "splitCyclicPolyPatch::calcTransforms : Writing " << name()
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, half0, half0.points());
        }
        const splitCyclicPolyPatch& half1 = neighbPatch();
        {
            fileName nm1(casePath/half1.name()+"_faces.obj");
            Pout<< "splitCyclicPolyPatch::calcTransforms : Writing " << half1.name()
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, half1, half1.points());
        }
        {
            OFstream str(casePath/name()+"_to_" + half1.name() + ".obj");
            label vertI = 0;
            Pout<< "splitCyclicPolyPatch::calcTransforms :"
                << " Writing coupled face centres as lines to " << str.name()
                << endl;
            forAll(half0Ctrs, i)
            {
                const point& p0 = half0Ctrs[i];
                str << "v " << p0.x() << ' ' << p0.y() << ' ' << p0.z() << nl;
                vertI++;
                const point& p1 = half1Ctrs[i];
                str << "v " << p1.x() << ' ' << p1.y() << ' ' << p1.z() << nl;
                vertI++;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }


    // Some sanity checks

    if (half0Ctrs.size() != half1Ctrs.size())
    {
        FatalErrorIn
        (
            "splitCyclicPolyPatch::calcTransforms()"
        )   << "For patch " << name()
            << " there are " << half0Ctrs.size()
            << " face centres, for the neighbour patch " << neighbPatch().name()
            << " there are " << half1Ctrs.size()
            << exit(FatalError);
    }

    /*if (transform() != neighbPatch().transform())
    {
        FatalErrorIn
        (
            "splitCyclicPolyPatch::calcTransforms()"
        )   << "Patch " << name()
            << " has transform type " << transformTypeNames[transform()]
            << ", neighbour patch " << neighbPatchName()
            << " has transform type "
            << neighbPatch().transformTypeNames[neighbPatch().transform()]
            << exit(FatalError);
    }*/

    // Calculate transformation tensors

    if (half0Ctrs.size() > 0)
    {
        vectorField half0Normals(half0Areas.size());
        vectorField half1Normals(half1Areas.size());

        scalar maxAreaDiff = -GREAT;
        label maxAreaFacei = -1;

        forAll(half0, facei)
        {
            scalar magSf = mag(half0Areas[facei]);
            scalar nbrMagSf = mag(half1Areas[facei]);
            scalar avSf = (magSf + nbrMagSf)/2.0;

            if (magSf < ROOTVSMALL && nbrMagSf < ROOTVSMALL)
            {
                // Undetermined normal. Use dummy normal to force separation
                // check. (note use of sqrt(VSMALL) since that is how mag
                // scales)
                half0Normals[facei] = point(1, 0, 0);
                half1Normals[facei] = half0Normals[facei];
            }
            else
            {
                scalar areaDiff = mag(magSf - nbrMagSf)/avSf;

                if (areaDiff > maxAreaDiff)
                {
                    maxAreaDiff = areaDiff;
                    maxAreaFacei = facei;
                }

                if (areaDiff > polyPatch::matchTol_())
                {
                    FatalErrorIn
                    (
                        "splitCyclicPolyPatch::calcTransforms()"
                    )   << "face " << facei
                        << " area does not match neighbour by "
                        << 100*areaDiff
                        << "% -- possible face ordering problem." << endl
                        << "patch:" << name()
                        << " my area:" << magSf
                        << " neighbour area:" << nbrMagSf
                        << " matching tolerance:" << polyPatch::matchTol_()
                         << endl
                        << "Mesh face:" << start()+facei
                        << " fc:" << half0Ctrs[facei]
                        << endl
                        << "Neighbour fc:" << half1Ctrs[facei]
                        << endl
                        << "If you are certain your matching is correct"
                        << " you can increase the 'matchTolerance' setting"
                        << " in the patch dictionary in the boundary file."
                        << endl
                        << "Rerun with cyclic debug flag set"
                        << " for more information." << exit(FatalError);
                }
                else
                {
                    half0Normals[facei] = half0Areas[facei] / magSf;
                    half1Normals[facei] = half1Areas[facei] / nbrMagSf;
                }
            }
        }


        // Print area match
        if (debug)
        {
            Pout<< "splitCyclicPolyPatch::calcTransforms :"
                << " patch:" << name()
                << " Max area error:" << 100*maxAreaDiff << "% at face:"
                << maxAreaFacei << " at:" << half0Ctrs[maxAreaFacei]
                << " coupled face at:" << half1Ctrs[maxAreaFacei]
                << endl;
        }


        // Calculate transformation tensors

        if (transform() == ROTATIONAL)
        {
            // Calculate using the given rotation axis and centre. Do not
            // use calculated normals.
            vector n0 = findFaceMaxRadius(half0Ctrs);
            vector n1 = -findFaceMaxRadius(half1Ctrs);
            n0 /= mag(n0) + VSMALL;
            n1 /= mag(n1) + VSMALL;

            if (debug)
            {
                scalar theta = radToDeg(acos(n0 & n1));

                Pout<< "splitCyclicPolyPatch::calcTransforms :"
                    << " patch:" << name()
                    << " Specified rotation :"
                    << " n0:" << n0 << " n1:" << n1
                    << " swept angle: " << theta << " [deg]"
                    << endl;
            }

            // Extended tensor from two local coordinate systems calculated
            // using normal and rotation axis
            const tensor E0
            (
                rotationAxis_,
                (n0 ^ rotationAxis_),
                n0
            );
            const tensor E1
            (
                rotationAxis_,
                (-n1 ^ rotationAxis_),
                -n1
            );
            const tensor revT(E1.T() & E0);

            const_cast<tensorField&>(forwardT()) = tensorField(1, revT.T());
            const_cast<tensorField&>(reverseT()) = tensorField(1, revT);
            const_cast<vectorField&>(separation()).setSize(0);
            //const_cast<boolList&>(collocated()) = boolList(1, false);
        }
        else
        {
            scalarField half0Tols
            (
                polyPatch::matchTol_()
               *calcFaceTol
                (
                    half0,
                    half0.points(),
                    static_cast<const pointField&>(half0Ctrs)
                )
            );

            calcTransformTensors
            (
                static_cast<const pointField&>(half0Ctrs),
                static_cast<const pointField&>(half1Ctrs),
                half0Normals,
                half1Normals,
                half0Tols
                //polyPatch::matchTol_(),
                //transform()
            );


            if (transform() == TRANSLATIONAL)
            {
                if (debug)
                {
                    Pout<< "splitCyclicPolyPatch::calcTransforms :"
                        << " patch:" << name()
                        << " Specified separation vector : "
                        << separationVector_ << endl;
                }

                // Check that separation vectors are same.
                const scalar avgTol = average(half0Tols);
                if
                (
                    mag(separationVector_ + neighbPatch().separationVector_)
                  > avgTol
                )
                {
                    WarningIn
                    (
                        "splitCyclicPolyPatch::calcTransforms()"
                    )   << "Specified separation vector " << separationVector_
                        << " differs by that of neighbouring patch "
                        << neighbPatch().separationVector_
                        << " by more than tolerance " << avgTol << endl
                        << "patch:" << name()
                        << " neighbour:" << neighbPatchName()
                        << endl;
                }


                // Override computed transform with specified.
                if
                (
                    separation().size() != 1
                 || mag(separation()[0] - separationVector_) > avgTol
                )
                {
                    WarningIn
                    (
                        "splitCyclicPolyPatch::calcTransforms()"
                    )   << "Specified separationVector " << separationVector_
                        << " differs from computed separation vector "
                        << separation() << endl
                        << "This probably means your geometry is not consistent"
                        << " with the specified separation and might lead"
                        << " to problems." << endl
                        << "Continuing with specified separation vector "
                        << separationVector_ << endl
                        << "patch:" << name()
                        << " neighbour:" << neighbPatchName()
                        << endl;
                }

                // Set tensors
                const_cast<tensorField&>(forwardT()).clear();
                const_cast<tensorField&>(reverseT()).clear();
                const_cast<vectorField&>(separation()) = vectorField
                (
                    1,
                    separationVector_
                );
                //const_cast<boolList&>(collocated()) = boolList(1, false);
            }
        }
    }
}


void Foam::splitCyclicPolyPatch::getCentresAndAnchors
(
    const primitivePatch& pp0,
    const primitivePatch& pp1,

    pointField& half0Ctrs,
    pointField& half1Ctrs,
    pointField& anchors0,
    scalarField& tols
) const
{
    // Get geometric data on both halves.
    half0Ctrs = pp0.faceCentres();
    //anchors0 = getAnchorPoints(pp0, pp0.points(), transform());
    anchors0 = getAnchorPoints(pp0, pp0.points());
    half1Ctrs = pp1.faceCentres();

    if (debug)
    {
        Pout<< "splitCyclicPolyPatch::getCentresAndAnchors :"
            << " patch:" << name() << nl
            << "half0 untransformed faceCentres (avg) : "
            << gAverage(half0Ctrs) << nl
            << "half1 untransformed faceCentres (avg) : "
            << gAverage(half1Ctrs) << endl;
    }

    if (half0Ctrs.size())
    {
        switch (transform())
        {
            case ROTATIONAL:
            {
                vector n0 = findFaceMaxRadius(half0Ctrs);
                vector n1 = -findFaceMaxRadius(half1Ctrs);
                n0 /= mag(n0) + VSMALL;
                n1 /= mag(n1) + VSMALL;

                if (debug)
                {
                    scalar theta = radToDeg(acos(n0 & n1));

                    Pout<< "splitCyclicPolyPatch::getCentresAndAnchors :"
                        << " patch:" << name()
                        << " Specified rotation :"
                        << " n0:" << n0 << " n1:" << n1
                        << " swept angle: " << theta << " [deg]"
                        << endl;
                }

                // Extended tensor from two local coordinate systems calculated
                // using normal and rotation axis
                const tensor E0
                (
                    rotationAxis_,
                    (n0 ^ rotationAxis_),
                    n0
                );
                const tensor E1
                (
                    rotationAxis_,
                    (-n1 ^ rotationAxis_),
                    -n1
                );
                const tensor revT(E1.T() & E0);

                // Rotation
                forAll(half0Ctrs, faceI)
                {
                    half0Ctrs[faceI] =
                        Foam::transform
                        (
                            revT,
                            half0Ctrs[faceI] - rotationCentre_
                        )
                      + rotationCentre_;
                    anchors0[faceI] =
                        Foam::transform
                        (
                            revT,
                            anchors0[faceI] - rotationCentre_
                        )
                      + rotationCentre_;
                }

                break;
            }
            case TRANSLATIONAL:
            {
                // Transform 0 points.

                if (debug)
                {
                    Pout<< "splitCyclicPolyPatch::getCentresAndAnchors :"
                        << " patch:" << name()
                        << "Specified translation : " << separationVector_
                        << endl;
                }

                // Note: getCentresAndAnchors gets called on the slave side
                // so separationVector is owner-slave points.

                half0Ctrs -= separationVector_;
                anchors0 -= separationVector_;
                break;
            }
            default:
            {
                // Assumes that cyclic is rotational. This is also the initial
                // condition for patches without faces.

                // Determine the face with max area on both halves. These
                // two faces are used to determine the transformation tensors
                label max0I = findMaxArea(pp0.points(), pp0);
                vector n0 = pp0[max0I].normal(pp0.points());
                n0 /= mag(n0) + VSMALL;

                label max1I = findMaxArea(pp1.points(), pp1);
                vector n1 = pp1[max1I].normal(pp1.points());
                n1 /= mag(n1) + VSMALL;

                if (mag(n0 & n1) < 1-polyPatch::matchTol_())
                {
                    if (debug)
                    {
                        Pout<< "splitCyclicPolyPatch::getCentresAndAnchors :"
                            << " patch:" << name()
                            << " Detected rotation :"
                            << " n0:" << n0 << " n1:" << n1 << endl;
                    }

                    // Rotation (around origin)
                    const tensor revT(rotationTensor(n0, -n1));

                    // Rotation
                    forAll(half0Ctrs, faceI)
                    {
                        half0Ctrs[faceI] = Foam::transform
                        (
                            revT,
                            half0Ctrs[faceI]
                        );
                        anchors0[faceI] = Foam::transform
                        (
                            revT,
                            anchors0[faceI]
                        );
                    }
                }
                else
                {
                    // Parallel translation. Get average of all used points.

                    const point ctr0(sum(pp0.localPoints())/pp0.nPoints());
                    const point ctr1(sum(pp1.localPoints())/pp1.nPoints());

                    if (debug)
                    {
                        Pout<< "splitCyclicPolyPatch::getCentresAndAnchors :"
                            << " patch:" << name()
                            << " Detected translation :"
                            << " n0:" << n0 << " n1:" << n1
                            << " ctr0:" << ctr0 << " ctr1:" << ctr1 << endl;
                    }

                    half0Ctrs += ctr1 - ctr0;
                    anchors0 += ctr1 - ctr0;
                }
                break;
            }
        }
    }

    // Calculate typical distance per face
    tols = polyPatch::matchTol_()*calcFaceTol(pp1, pp1.points(), half1Ctrs);
}


Foam::vector Foam::splitCyclicPolyPatch::findFaceMaxRadius
(
    const pointField& faceCentres
) const
{
    // Determine a face furthest away from the axis

    const vectorField n((faceCentres - rotationCentre_) ^ rotationAxis_);

    const scalarField magRadSqr(magSqr(n));

    label faceI = findMax(magRadSqr);

    if (debug)
    {
        Info<< "findFaceMaxRadius(const pointField&) : patch: " << name() << nl
            << "    rotFace  = " << faceI << nl
            << "    point    = " << faceCentres[faceI] << nl
            << "    distance = " << Foam::sqrt(magRadSqr[faceI])
            << endl;
    }

    return n[faceI];
}


void Foam::splitCyclicPolyPatch::initAddressing()
{
    polyPatch::initAddressing();
}


void Foam::splitCyclicPolyPatch::calcAddressing()
{
    polyPatch::calcAddressing();
}


void Foam::splitCyclicPolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}


void Foam::splitCyclicPolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();

    calcTransforms();
}


void Foam::splitCyclicPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);
}


void Foam::splitCyclicPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    calcTransforms();
}


void Foam::splitCyclicPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();
}


void Foam::splitCyclicPolyPatch::updateMesh()
{
    polyPatch::updateMesh();
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}


void Foam::splitCyclicPolyPatch::syncOrder() const
{}


void Foam::splitCyclicPolyPatch::initOrder(const primitivePatch& pp) const
{
    if (master())
    {
        // Save patch for use in non-owner side ordering. Equivalent to
        // processorPolyPatch using OPstream.
        ownerPatchPtr_.reset
        (
            new primitivePatch
            (
                pp,
                pp.points()
            )
        );
    }
}


bool Foam::splitCyclicPolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    if (debug)
    {
        Pout<< "order : of " << pp.size()
            << " faces of patch:" << name()
            << " neighbour:" << neighbPatchName()
            << endl;
    }
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    if (transform() == NOORDERING)
    {
        // No faces, nothing to change.
        return false;
    }

    if (master())
    {
        // Do nothing (i.e. identical mapping, zero rotation).
        // See comment at top.
        forAll(faceMap, patchFaceI)
        {
            faceMap[patchFaceI] = patchFaceI;
        }

        return false;
    }
    else
    {
        // Get stored geometry from initOrder invocation of owner.
        const primitivePatch& pp0 = neighbPatch().ownerPatchPtr_();

        // Get geometric quantities
        pointField half0Ctrs, half1Ctrs, anchors0;
        scalarField tols;
        getCentresAndAnchors
        (
            pp0,
            pp,

            half0Ctrs,
            half1Ctrs,
            anchors0,
            tols
        );

        if (debug)
        {
            Pout<< "half0 transformed faceCentres (avg)   : "
                << gAverage(half0Ctrs) << nl
                << "half1 untransformed faceCentres (avg) : "
                << gAverage(half1Ctrs) << endl;
        }

        // Geometric match of face centre vectors
        bool matchedAll = matchPoints
        (
            half1Ctrs,
            half0Ctrs,
            tols,
            true,
            faceMap
        );

        if (!matchedAll || debug)
        {
            // Dump halves
            fileName nm0
            (
                boundaryMesh().mesh().time().path()
               /neighbPatch().name()+"_faces.obj"
            );
            Pout<< "cyclicPolyPatch::order : Writing neighbour"
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, pp0, pp0.points());

            fileName nm1
            (
                boundaryMesh().mesh().time().path()
               /name()+"_faces.obj"
            );
            Pout<< "cyclicPolyPatch::order : Writing my"
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, pp, pp.points());

            OFstream ccStr
            (
                boundaryMesh().mesh().time().path()
               /name() + "_faceCentres.obj"
            );
            Pout<< "cyclicPolyPatch::order : "
                << "Dumping currently found cyclic match as lines between"
                << " corresponding face centres to file " << ccStr.name()
                << endl;

            // Recalculate untransformed face centres
            //pointField rawHalf0Ctrs =
            //    calcFaceCentres(half0Faces, pp.points());
            label vertI = 0;

            forAll(half1Ctrs, i)
            {
                if (faceMap[i] != -1)
                {
                    // Write edge between c1 and c0
                    const point& c0 = half0Ctrs[faceMap[i]];
                    const point& c1 = half1Ctrs[i];
                    writeOBJ(ccStr, c0, c1, vertI);
                }
            }
        }

        if (!matchedAll)
        {
            SeriousErrorIn
            (
                "cyclicPolyPatch::order"
                "(const primitivePatch&, labelList&, labelList&) const"
            )   << "Patch:" << name() << " : "
                << "Cannot match vectors to faces on both sides of patch"
                << endl
                << "    Perhaps your faces do not match?"
                << " The obj files written contain the current match." << endl
                << "    Continuing with incorrect face ordering from now on!"
                << endl;

                return false;
        }


        // Set rotation.
        forAll(faceMap, oldFaceI)
        {
            // The face f will be at newFaceI (after morphing) and we want its
            // anchorPoint (= f[0]) to align with the anchorpoint for the
            // corresponding face on the other side.

            label newFaceI = faceMap[oldFaceI];

            const point& wantedAnchor = anchors0[newFaceI];

            rotation[newFaceI] = getRotation
            (
                pp.points(),
                pp[oldFaceI],
                wantedAnchor,
                tols[oldFaceI]
            );

            if (rotation[newFaceI] == -1)
            {
                SeriousErrorIn
                (
                    "cyclicPolyPatch::order(const primitivePatch&"
                    ", labelList&, labelList&) const"
                )   << "in patch " << name()
                    << " : "
                    << "Cannot find point on face " << pp[oldFaceI]
                    << " with vertices "
                    << IndirectList<point>(pp.points(), pp[oldFaceI])()
                    << " that matches point " << wantedAnchor
                    << " when matching the halves of processor patch " << name()
                    << "Continuing with incorrect face ordering from now on!"
                    << endl;

                return false;
            }
        }

        ownerPatchPtr_.clear();

        // Return false if no change neccesary, true otherwise.

        forAll(faceMap, faceI)
        {
            if (faceMap[faceI] != faceI || rotation[faceI] != 0)
            {
                return true;
            }
        }

        return false;
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::splitCyclicPolyPatch::splitCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm),
    neighbPatchName_(word::null),
    neighbPatchID_(-1),
    rotationAxis_(vector::zero),
    rotationCentre_(point::zero),
    separationVector_(vector::zero),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::splitCyclicPolyPatch::splitCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    neighbPatchName_(dict.lookupOrDefault("neighbourPatch", word::null)),
    coupleGroup_(dict),
    neighbPatchID_(-1),
    rotationAxis_(vector::zero),
    rotationCentre_(point::zero),
    separationVector_(vector::zero),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    if (neighbPatchName_ == word::null && !coupleGroup_.valid())
    {
        FatalIOErrorIn
        (
            "cyclicPolyPatch::cyclicPolyPatch\n"
            "(\n"
            "    const word& name,\n"
            "    const dictionary& dict,\n"
            "    const label index,\n"
            "    const polyBoundaryMesh& bm\n"
            ")",
            dict
        )   << "No \"neighbourPatch\" provided." << endl
            << "Is your mesh uptodate with split cyclics?" << endl
            << "Run foamUpgradeCyclics to convert mesh and fields"
            << " to split cyclics." << exit(FatalIOError);
    }

    if (neighbPatchName_ == name)
    {
        FatalIOErrorIn("cyclicPolyPatch::cyclicPolyPatch(..)", dict)
            << "Neighbour patch name " << neighbPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    switch (transform())
    {
        case ROTATIONAL:
        {
            dict.lookup("rotationAxis") >> rotationAxis_;
            dict.lookup("rotationCentre") >> rotationCentre_;

            scalar magRot = mag(rotationAxis_);
            if (magRot < SMALL)
            {
                FatalIOErrorIn("cyclicPolyPatch::cyclicPolyPatch(..)", dict)
                    << "Illegal rotationAxis " << rotationAxis_ << endl
                    << "Please supply a non-zero vector."
                    << exit(FatalIOError);
            }
            rotationAxis_ /= magRot;

            break;
        }
        case TRANSLATIONAL:
        {
            dict.lookup("separationVector") >> separationVector_;
            break;
        }
        default:
        {
            // no additional info required
        }
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}



Foam::splitCyclicPolyPatch::splitCyclicPolyPatch
(
    const splitCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    neighbPatchName_(pp.neighbPatchName()),
    coupleGroup_(pp.coupleGroup_),
    neighbPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    separationVector_(pp.separationVector_),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    if (pp.neighbPatchName() == name())
    {
        FatalErrorIn("cyclicPolyPatch::cyclicPolyPatch(..)")
            << "Neighbour patch name " << neighbPatchName()
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::splitCyclicPolyPatch::~splitCyclicPolyPatch()
{
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::splitCyclicPolyPatch::neighbPatchName() const
{
    if (neighbPatchName_.empty())
    {
        // Try and use patchGroup to find samplePatch and sampleRegion
        label patchID = coupleGroup_.findOtherPatchID(*this);

        neighbPatchName_ = boundaryMesh()[patchID].name();
    }
    return neighbPatchName_;
}


Foam::label Foam::splitCyclicPolyPatch::neighbPatchID() const
{
    if (neighbPatchID_ == -1)
    {
        neighbPatchID_ = this->boundaryMesh().findPatchID(neighbPatchName());

        if (neighbPatchID_ == -1)
        {
            FatalErrorIn("cyclicPolyPatch::neighbPatchID() const")
                << "Illegal neighbourPatch name " << neighbPatchName()
                << endl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a cyclic
        const splitCyclicPolyPatch& nbrPatch = refCast<const splitCyclicPolyPatch>
        (
            this->boundaryMesh()[neighbPatchID_]
        );

        if (nbrPatch.neighbPatchName() != name())
        {
            WarningIn("cyclicPolyPatch::neighbPatchID() const")
                << "Patch " << name()
                << " specifies neighbour patch " << neighbPatchName()
                << endl << " but that in return specifies "
                << nbrPatch.neighbPatchName()
                << endl;
        }
    }
    return neighbPatchID_;
}


void Foam::splitCyclicPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    if (!neighbPatchName_.empty())
    {
        os.writeKeyword("neighbourPatch") << neighbPatchName_
            << token::END_STATEMENT << nl;
    }
    coupleGroup_.write(os);
    switch (transform_)
    {
        case ROTATIONAL:
        {
            os.writeKeyword("rotationAxis") << rotationAxis_
                << token::END_STATEMENT << nl;
            os.writeKeyword("rotationCentre") << rotationCentre_
                << token::END_STATEMENT << nl;
            break;
        }
        case TRANSLATIONAL:
        {
            os.writeKeyword("separationVector") << separationVector_
                << token::END_STATEMENT << nl;
            break;
        }
        case NOORDERING:
        {
            break;
        }
        default:
        {
            // no additional info to write
        }
    }
}


// ************************************************************************* //
