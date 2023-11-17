functions
(
    Probes
    {
        type               probes;
        functionObjectLibs ("libsampling.so");
        setFormat          gnuplot;
        probeLocations
	    (
            ($xm 0 -0.050000)
            ($xm 0 -0.100000)
            ($xm 0 -0.150000)
            ($xm 0 -0.180000)
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

