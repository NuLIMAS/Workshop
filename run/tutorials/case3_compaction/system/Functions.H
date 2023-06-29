functions
(
    Probes 
    {
        type               probes;
        functionObjectLibs ("libsampling.so");
        setFormat          gnuplot;
        probeLocations
	(
          (1 0 -0.050000)
          (1 0 -0.100000)
          (1 0 -0.150000)
          (1 0 -0.180000)

	);
	fields
	(
	    p
	    pMean
	    liqueCriteria
	);
    }

    Profiles
    {
        type                   sets;
	functionObjectLibs     ("libsampling.so");
	setFormat              raw;
	interpolationScheme    cellPointFace;
	outputControl          outputTime;
	sets
	(
            Profile1
            {
                type        uniform;
                axis        z;
                start       (1 0 0);
                end         (1 0 -0.2);
                nPoints     20;
            }
            Profile2
            {
                type        uniform;
                axis        z;
                start       (1 0 -0.0);
                end         (1 0 -0.2);
                nPoints     200;
            }
            Profile3
            {
                type        uniform;
                axis        x;
                start       (0    0 -0.2);
                end         (2 0 -0.2);
                nPoints     200;
            }
            
        );
        fields
        (
            p
            tauXZ
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
	        mean 	 yes;
		prime2Mean 	 no;
		base 	 time;
	    }
	    tauXZ
 	    {
	        mean 	 yes;
		prime2Mean 	 yes;
		base 	 time;
	    }
	);
    }


);

