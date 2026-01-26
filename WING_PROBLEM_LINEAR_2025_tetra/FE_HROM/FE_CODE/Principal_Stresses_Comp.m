function [PRINCIPAL_STRESSES stressVONMISES ] =  Principal_Stresses_Comp(DATAOUT,MATERIAL) ; 
% Computing principal stresses for both fibers and matrix
%dbstop('4')
if nargin == 0 
    load('tmp.mat')
end

nstress = 6 ; 
ngaus =  length(DATAOUT.stress)/nstress;  
nelem = length(DATAOUT.MaterialType) ; 
ngausE = ngaus/nelem ; 
%%% AVERAGE STRESSES ON ELEMENTS (GIVEN STRESSES AT ALL GAUSS POINTS)
stressEL = Stress_Elements_Average(DATAOUT.wSTs,DATAOUT.stress,ngausE,nstress,nelem) ; 

% STRESSm = zeros(3,3,nelem) ; 
% IND = cell(6,1) ; 
% for i=1:nstress
%     IND{i} = i:nstress:nstress*nelem ; 
% end

% Correspondence between indices (Voigt)
STRESSm(1,1,:) =  stressEL(1,:) ; 
STRESSm(2,2,:) =  stressEL(2,:) ; 
STRESSm(3,3,:) =  stressEL(3,:) ; 
STRESSm(2,3,:) =  stressEL(4,:) ; 
STRESSm(3,2,:) =  stressEL(4,:) ; 
STRESSm(1,3,:) =  stressEL(5,:) ; 
STRESSm(3,1,:) =  stressEL(5,:) ; 
STRESSm(1,2,:) =  stressEL(6,:) ; 
STRESSm(2,1,:) =  stressEL(6,:) ;  

%% Principal stresses (vectorized form)
PSTRESS = eig3(STRESSm) ;
%% 

%% VON MISES CRITERIA 
sigmaVONMISES = 0.5*((PSTRESS(1,:)-PSTRESS(2,:)).^2 +...
    (PSTRESS(1,:)-PSTRESS(3,:)).^2 + (PSTRESS(3,:)-PSTRESS(2,:)).^2) ; 
sigmaVONMISES = sqrt(sigmaVONMISES) ; 

% Principal stresses fibers 

elemFIB = find(DATAOUT.MaterialType == MATERIAL.FIBER.INDEX) ;  % Fiber elements
elemMAT= find(DATAOUT.MaterialType == MATERIAL.MATRIX.INDEX) ;  % Matrix elements 

% Corresponding gauss points FIBER 
%gaussFIB = GaussPointsAssociatedElements(ngausE,elemFIB) ; 
% Corresponding gauss points MATRIX 
%gaussMAT = GaussPointsAssociatedElements(ngausE,elemMAT) ; 



PRINCIPAL_STRESSES.FIBER = PSTRESS(:,elemFIB) ; 
PRINCIPAL_STRESSES.MATRIX = PSTRESS(:,elemMAT) ; 
stressVONMISES.FIBER = sigmaVONMISES(elemFIB) ; 
stressVONMISES.MATRIX = sigmaVONMISES(elemMAT) ; 
