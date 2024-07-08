functions
(
	 Probes 
	{
		type               probes;
		functionObjectLibs ("libsampling.so");
		setFormat          gnuplot;
		probeLocations     $probes;
		fields             $field;
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

