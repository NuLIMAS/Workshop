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
				start       (0.45 0 0);
				end         (0.45 0 -0.17);
				nPoints     200;
			}

		

		);
		fields
		(
			 E
			 sigma0
		);
	}

);

