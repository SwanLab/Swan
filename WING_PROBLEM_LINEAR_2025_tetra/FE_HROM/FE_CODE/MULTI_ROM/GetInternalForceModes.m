function  [SNAPforceS,MSG] = GetInternalForceModes(BASES,DATAIN,MSG,BdomRED,Wdom,DATAROM,Bdom)


BasisS = BASES.STRESSES.U ;
MSG{end+1} = ['Number of stress modes  = ',num2str(size(BasisS,2))];
ngaus = length(Wdom);
nstrain = size(BdomRED,1)/length(Wdom);


    SNAPforceS = BasisRED_BasisStress(BdomRED,BasisS,nstrain,ngaus) ;
    %--------------------------------------------------------------

