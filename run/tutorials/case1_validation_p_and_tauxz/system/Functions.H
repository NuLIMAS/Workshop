#include "$FOAM_CASE/caseInfo"

functions
(

    Profiles
    {
        type                   sets;
        functionObjectLibs     ("libsampling.so");
        setFormat              raw;
        interpolationScheme    cellPoint;
        outputControl          outputTime;

        sets
        (
            Profile1
            {
                type        uniform;
                axis        z;
                start       ($xm 0 $z1);
                end         ($xm 0 $z2);
                nPoints     200;
            }

            Profile2
            {
                type        uniform;
                axis        z;
                start       ($xq 0 $z1);
                end         ($xq 0 $z2);
                nPoints     200;
            }
        );
        fields
        (
            p
            tauXZ
            U
        );
    }

);

