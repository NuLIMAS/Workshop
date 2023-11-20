functions
(
    Probes
    {
        type               probes;
        functionObjectLibs ("libsampling.so");
        setFormat          gnuplot;
        probeLocations     $probes;
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

