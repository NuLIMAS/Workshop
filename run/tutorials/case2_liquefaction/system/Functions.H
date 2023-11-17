#include "$FOAM_CASE/caseInfo"

functions
(
    Probes
    {
        type               probes;
        functionObjectLibs ("libsampling.so");
        setFormat          gnuplot;
        probeLocations
        (
            ($xm 0.0 -0.085)
            ($xm 0.0 -0.185)
            ($xm 0.0 -0.24)
            ($xm 0.0 -0.3)
            ($xm 0.0 -0.4)
        );
        fields
        (
            pE
            p
            liqueCriteria
        );
    }

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
                nPoints     100;
            }
        );
        fields
        (
            pE
            liqueCriteria
        );
    }

    fieldAverage
    {
        type                   fieldAverage;
        functionObjectLibs     ("libfieldFunctionObjects.so");
        setFormat              raw;
        outputControl          outputTime;
        fields
        (
            tauXZ
            {
                mean 	 yes;
                prime2Mean 	 yes;
                base 	 time;
            }
        );
    }
);

