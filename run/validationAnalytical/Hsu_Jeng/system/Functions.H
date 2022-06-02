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
				start       (61.0375 0 0);
				end         (61.0375 0 -50);
				nPoints     200;
			}

			Profile2
			{
				type        uniform;
				axis        z;
				start       (30.51875 0 0);
				end         (30.51875 0 -50);
				nPoints     200;
			}
		

		);
		fields
		(
			 p
			 U
			 tauXZ
		);
	}


);

