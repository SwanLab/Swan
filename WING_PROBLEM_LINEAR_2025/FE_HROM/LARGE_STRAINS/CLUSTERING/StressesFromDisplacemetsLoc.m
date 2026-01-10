function [SNAPstressSTWOproj_LOC,SNAPstressPonePROJ_LOC] = StressesFromDisplacemetsLoc(IND_LOCAL_CLUSTER,icluster,BasisU_cluster,SNAPdisp,DOFr,DOFl,OPERFE,MATPRO,DATA)
if nargin == 0
    load('tmp.mat')
end



if isstruct(SNAPdisp)
    coeff_FEast = BasisU_cluster.coeff{icluster}'*SNAPdisp.DOFl.coeff(:,IND_LOCAL_CLUSTER) ;
    dL_ast =  BasisU_cluster.coeff{icluster}*coeff_FEast;
    dL = BasisU_cluster.BASIS*dL_ast ;
    
    coeff_FEast =  SNAPdisp.DOFr.coeff(:,IND_LOCAL_CLUSTER) ;
    dR = SNAPdisp.DOFr.BASIS*coeff_FEast;
    
else
    coeff = BasisU_cluster{icluster}'*SNAPdisp(DOFl,IND_LOCAL_CLUSTER) ;
    % ------------------------------------------------
    dL = BasisU_cluster{icluster}*coeff; % Snapshot displacements (local)
    dR = SNAPdisp(DOFr,IND_LOCAL_CLUSTER) ;
end


ndof = size(dL,1)+size(dR,1) ;
d = zeros(ndof,size(dL,2)) ;
d(DOFl,:)  = dL ;
d(DOFr,:)  = dR ;



%
% 2. Deformation gradient at all Gauss points
FgradST = OPERFE.Bst*d + repmat(OPERFE.IDENTITY_F,1,size(d,2)) ;
% 3. Green-Lagrante strains at all Gauss points
GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
% 4. 2nd Piola-Kirchhoff stresses at all Gauss Points
SNAPstressSTWOproj_LOC = zeros(size(GLSTRAINS)) ;
for isnap = 1:size(GLSTRAINS,2)
    [SNAPstressSTWOproj_LOC(:,isnap) ]= PK2stress_Constitutive_Model(GLSTRAINS(:,isnap),MATPRO,DATA,FgradST(:,isnap)) ;
end
% 5. 1st Piola-Kirchhoff stresses at all Gauss Points
SNAPstressPonePROJ_LOC = PK1stress(SNAPstressSTWOproj_LOC,FgradST,DATA.MESH.ndim) ;