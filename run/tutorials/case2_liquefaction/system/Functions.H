functions
(
    Probes
    {
        type               probes;
        functionObjectLibs ("libsampling.so");
        setFormat          gnuplot;
        probeLocations
        (
            (1.59 0.0 -0.085)
            (1.59 0.0 -0.185)
            (1.59 0.0 -0.24)
            (1.59 0.0 -0.3)
            (1.59 0.0 -0.4)
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
                start       (1.59 0 0);
                end         (1.59 0 -0.4);
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

