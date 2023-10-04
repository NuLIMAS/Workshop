    if (runTime.outputTime())
    {
        sigma0 = -(sigmaD.component(symmTensor::XX) + sigmaD.component(symmTensor::YY) + sigmaD.component(symmTensor::ZZ))/3;
        runTime.write(); 
    }
