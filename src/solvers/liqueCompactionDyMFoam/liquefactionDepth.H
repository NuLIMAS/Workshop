scalar clmin = 0;
scalar clmax = min(mesh.C().component(vector::Z)).value();
scalar csmin = min(mesh.C().component(vector::Z)).value();
scalar csmax = min(mesh.C().component(vector::Z)).value();

if (min(Cl.internalField()) ==1)
{
    clmin = min(mesh.C().component(vector::Z)).value();
    clmax = min(mesh.C().component(vector::Z)).value();
    csmin = min(mesh.C().component(vector::Z)).value();
    csmax = min(mesh.C().component(vector::Z)).value();
}
else
{
    if (max(Cl.internalField()) ==0)
    {
        clmin = max(mesh.C().component(vector::Z)).value();
        clmax = max(mesh.C().component(vector::Z)).value();
        csmin = min(mesh.C().component(vector::Z)).value();
        csmax = min(mesh.C().component(vector::Z)).value();
    }
    else
    {
        clmin = 0;
        clmax = min(mesh.C().component(vector::Z)).value();
        csmax= min(mesh.C().component(vector::Z)).value();
        csmin = 0;
        forAll(Cl, cellI)
        {
            if ( Cl[cellI] == 1)
            {
                clmin = min(clmin,mesh.C()[cellI].component(vector::Z));
            }
            else //if ( Cl[cellI] == 0)
            {
                clmax = max(clmax,mesh.C()[cellI].component(vector::Z));
            }
            if ( Cs[cellI] == 1)
            {
                csmax = max(csmax,mesh.C()[cellI].component(vector::Z));
            }
            else //if ( Cs[cellI] == 0)
            {
                csmin = min(csmin,mesh.C()[cellI].component(vector::Z));
            }

            }
    Info << "Partial liquefaction" << endl;
    }
}

std::ofstream file;
file.open ("liquefactionDepth.txt", std::ofstream::out | std::ofstream::app);
file << runTime.timeName() << " " << clmin << " " << clmax <<  " "  << csmin << " " << csmax  << std::endl << "\n";
file.close();