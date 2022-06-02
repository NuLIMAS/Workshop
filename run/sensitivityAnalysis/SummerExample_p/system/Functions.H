functions
(
	 Probes 
	{
		type               probes;
		functionObjectLibs ("libsampling.so");
		setFormat          gnuplot;
		probeLocations
		(
			 (87 0 0) ( 87 0 -0.2) (87 0 -0.4) (87 0 -0.75) ( 87 0 -1.0)
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
				start       (84 0 0.05);
				end         (84 0 -1);
				nPoints     100;
			}

		

		);
		fields
		(
			 pE
			 liqueCriteria
		);
	}


);

