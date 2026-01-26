function T_trained=computeT_Hybrid(basis,R,Q_NN,pol_deg)
    rFull = Data.buildModel(R,pol_deg);
    qNN=Q_NN.computeOutputValues(rFull).';
    T_NN=basis*qNN;
    T_trained=reshape(T_NN,[],8);
end