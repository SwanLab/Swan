function T_trained=computeT_NN(mesh,R,T_NN,pol_deg)
    nPoints   = size(mesh.coord,1);
    rRep      = repmat(R,nPoints,1);
    ndim      = mesh.ndim;
    dataInput = [rRep, mesh.coord];
    dataFull  = Data.buildModel(dataInput,pol_deg);
    T_all     = T_NN.computeOutputValues(dataFull);

    nF        = size(T_all,2) / mesh.ndim; % num coarse fun
    T_tmp     = reshape(T_all.', ndim, nF, nPoints);
    T_tmp = permute(T_tmp, [2 1 3]); 
    T_trained = reshape(T_tmp, nF, ndim*nPoints).';
end