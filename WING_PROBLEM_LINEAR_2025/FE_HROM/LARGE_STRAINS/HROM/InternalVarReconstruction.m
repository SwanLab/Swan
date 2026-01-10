function [OTHER_output,DATAinpGID,RECONS_PK2stress] = InternalVarReconstruction(DATAinpGID,OTHER_output,ECMdata,DATAHROM,BasisStwo)

% RECONSTRUCTION Internal Variables 
% ---------------------------------
DATAinpGID.OPERreconstr.PK2STRESS.BASIS =  BasisStwo ;
OTHER_output = DefaultField(OTHER_output,'BasisStwoZ',[]) ;
if isempty(OTHER_output.BasisStwoZ)
    if  isempty(ECMdata.setPoints)
        BasisStwoZ =  InterpolationGaussVariablesECM(BasisStwo,ECMdata,DATAHROM.MESH.ngaus_STRESS,DATAHROM.MESH.nstrain) ;
    else
        setIndices = small2large(ECMdata.setPoints,DATAHROM.MESH.nstrain) ;
        BasisStwoZ = BasisStwo(setIndices,:) ;
    end
else
    BasisStwoZ = OTHER_output.BasisStwoZ ;
end

DATAinpGID.OPERreconstr.PK2STRESS.coeff =  (BasisStwoZ'*BasisStwoZ)\BasisStwoZ' ;
RECONS_PK2stress.BASIS = BasisStwo ;
RECONS_PK2stress.coeff = DATAinpGID.OPERreconstr.PK2STRESS.coeff ;