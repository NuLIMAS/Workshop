functions
(
    Probes 
    {
        type               probes;
        functionObjectLibs ("libsampling.so");
        setFormat          gnuplot;
        probeLocations
	(
	    (61.0375 0 -0.01) 
	    ( 61.0375 0 -10) 
	    (61.0375 0 -20) 
	    (61.0375 0 -40) 
	    ( 61.0375 0 -49.99)
	);
	fields
	(
	    pE
            p	 
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
                start       (61.0375 0 0);
                end         (61.0375 0 -50);
                nPoints     200;
            }
        );
        fields
        (
            p
            D
            tauXZ
        );
    }
	
    flux
    {
        type        divFlux;
        phiName     phi;
        functionObjectLibs ("libutilityFunctionObjects.so");
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
		prime2Mean 	 yes;
		base 	 time;
	    }
	    tauXZ
 	    {
	        mean 	 yes;
		prime2Mean 	 yes;
		base 	 time;
	    }
	    phi 
 	    {
	        mean 	 yes;
		prime2Mean 	 yes;
		base 	 time;
	    }
 
	);
    }


);

