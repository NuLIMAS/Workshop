functions
(
    Probes
    {
        type               probes;
        functionObjectLibs ("libsampling.so");
        setFormat          gnuplot;
        probeLocations
        (
            ($xm 0 -0.05)
            ($xm 0 -0.10)
            ($xm 0 -0.15)
            ($xm 0 -0.18)
        );
        fields
        (
            p
            pMean
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
            p
            {
                mean         yes;
                prime2Mean   no;
                base         time;
            }
            tauXZ
            {
                mean         yes;
                prime2Mean   yes;
                base         time;
            }
        );
    }


);

