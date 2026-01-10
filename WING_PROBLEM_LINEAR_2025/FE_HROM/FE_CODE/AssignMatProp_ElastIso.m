function [celasglo celasgloINV typePROBLEM] = AssignMatProp_ElastIso(ndim,MATERIAL,nelem,MaterialType,DATA)

%dbstop('4')
if nargin == 0
    load('tmp.mat')
    DATA.typePROBLEM = 'pstress' ; 
elseif nargin == 4
    DATA = [] ; 
end

DATA =DefaultField(DATA,'typePROBLEM','3D')  ; 
typePROBLEM = DATA.typePROBLEM;  %'pstress'/'pstrain'/'3D';  Plane stress/ plane strain problem
if ndim==2
    nstrain = 3;
 
else
    nstrain = 6 ;
end

celasglo = zeros(nstrain,nstrain,nelem) ;  % Global array of elasticity matrices
celasgloINV = zeros(6,6,nelem) ;  % Global array of compliance matrices (3D)

NAMES_MAT = fieldnames(MATERIAL) ;
nmat = length(NAMES_MAT);

if isempty(MaterialType)
    MaterialType = ones(nelem,1) ;
end
%dbstop('27')
for imat = 1:nmat
    
    % Elastic properties 
    E = MATERIAL.(NAMES_MAT{imat}).E   ; %
    nu = MATERIAL.(NAMES_MAT{imat}).nu; % Poisson's coefficient
    Gshear = MATERIAL.(NAMES_MAT{imat}).Gshear; % Poisson's coefficient
    
    % Rotated fiber ? 
    if isfield(MATERIAL.(NAMES_MAT{imat}),'ANGLE')
        ANGLE = MATERIAL.(NAMES_MAT{imat}).ANGLE ; 
    else
        ANGLE = 0 ; 
    end
    
    switch DATA.typePROBLEM
        case '3D'
            [celasINV3D celas] = Compliance_ElasticityANGLE(E,nu,Gshear,typePROBLEM,ANGLE)  ;
        otherwise
            [celasINV3D celas] = Compliance_Elasticity(E,nu,Gshear,typePROBLEM)  ;
            
    end
    
    INDEXmat =   MATERIAL.(NAMES_MAT{imat}).INDEX  ;
    
    elemMAT = find(MaterialType == INDEXmat) ;
    
    
    
    for eloc=1:length(elemMAT)
        e = elemMAT(eloc) ;
        celasglo(:,:,e) = celas ;
        celasgloINV(:,:,e) = celasINV3D ;
    end
    
    
end