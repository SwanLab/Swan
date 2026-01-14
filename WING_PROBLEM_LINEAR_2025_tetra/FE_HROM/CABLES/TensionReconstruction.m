function [OTHER_output,DATAinpGID,RECONS_TENSIONV] = TensionReconstruction(DATAinpGID,OTHER_output,ECMdata,DATAHROM,BasisTENSIONV)

if nargin == 0
    load('tmp.mat')
end

% RECONSTRUCTION  TENSIONV
% -----------------------------
DATAinpGID.OPERreconstr.TENSIONV.BASIS =  BasisTENSIONV ;
OTHER_output = DefaultField(OTHER_output,'BasisTENSIONV_z',[]) ;
if isempty(OTHER_output.BasisTENSIONV_z)
    if  ~isfield(ECMdata,'setPoints') || isempty(ECMdata.setPoints)
        BasisTENSIONV_z =  InterpolationGaussVariablesECM(BasisTENSIONV,ECMdata,DATAHROM.MESH.ngaus_STRESS,DATAHROM.MESH.ndim) ;
    else
        setIndices = small2large(ECMdata.setPoints,DATAHROM.MESH.ndim) ;
        BasisTENSIONV_z = BasisTENSIONV(setIndices,:) ;
    end
else
    BasisTENSIONV_z = OTHER_output.BasisTENSIONV_z ;
end

DATAinpGID.OPERreconstr.TENSIONV.coeff =  (BasisTENSIONV_z'*BasisTENSIONV_z)\BasisTENSIONV_z' ;
RECONS_TENSIONV.BASIS = BasisTENSIONV ;
RECONS_TENSIONV.coeff = DATAinpGID.OPERreconstr.TENSIONV.coeff ;