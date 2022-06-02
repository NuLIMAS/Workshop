functions
(
	 Probes 
	{
		type               probes;
		functionObjectLibs ("libsampling.so");
		setFormat          gnuplot;
		probeLocations
		(
			 (0.45 0 0) (0.45 0 -0.05) (0.45 0 -0.10) (0.45 0 -0.14) (0.45 0 -0.17)
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
				start       (0.45 0 0);
				end         (0.45 0 -0.17);
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
		);
	}
);

