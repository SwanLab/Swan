function [celasglo celasgloINV typePROBLEM densGLO,DATA] = ...
    AssignMatProp_ElastIso(ndim,MATERIAL,nelem,MaterialType,DATA)

%dbstop('4')
if nargin == 0
    load('tmp1.mat')
elseif nargin == 4
    DATA = [] ; 
end

EXIST_DENS = 0 ;,
DATA = DefaultField(DATA,'typePROBLEM','pstress');  %'pstress'/'pstrain'/'3D';  Plane stress/ plane strain problem
DATA = DefaultField(DATA,'StrainStressWith4Components',1) ;
typePROBLEM = DATA.typePROBLEM ; 
if ndim==2 && DATA.StrainStressWith4Components == 0
    nstrain = 3;
elseif ndim==2 && DATA.StrainStressWith4Components == 1
    nstrain = 4 ;
else
    nstrain = 6 ;
    typePROBLEM ='3D' ;
end

DATA.ISPSTRAIN = 0; ;

celasglo = zeros(nstrain,nstrain,nelem) ;  % Global array of elasticity matrices
celasgloINV = zeros(6,6,nelem) ;  % Global array of compliance matrices (3D)
densGLO = ones(nelem,1) ;  % Global array of compliance matrices (3D)


NMAT_fe = length(unique(MaterialType));
nmat = length(MATERIAL.PLY) ;

if NMAT_fe > nmat
    for imat = nmat:NMAT_fe
        MATERIAL.PLY(imat) =  MATERIAL.PLY(nmat)   ;
    end
    nmat = NMAT_fe ;
end


if isempty(MaterialType)
    MaterialType = ones(nelem,1) ;
end


ISPSTRAIN = 0 ;


%dbstop('30')
for imat = 1:nmat
    % Load elasticity tensor
    MATLOC = MATERIAL.PLY(imat)  ;
    MATLOC=DefaultField(MATLOC,'NAMEWS',[]) ;
    MATLOC = DefaultField(MATLOC,'typePROBLEM','3D') ; 
    typePROBLEM = MATLOC.typePROBLEM  ;
    StrainStressWith4Components = 0 ;
    switch typePROBLEM   % Change 14-july-2019
        case '3D'
            if ndim ==2
                error('Set typePROBLEM to either pstress or pstrain')
            end
            
        case 'pstrain'
            StrainStressWith4Components = DATA.StrainStressWith4Components ;
            DATA.ISPSTRAIN = 1;
        case 'pstress'
            if DATA.StrainStressWith4Components == 1
                error('Incompatible option. ')
            end
    end
    
    if isempty(MATLOC.NAMEWS)
        E = MATLOC.E ;
        nu = MATLOC.nu ;
        G = MATLOC.G ;
   
        [celasINV3D Celas] = Compliance_Elasticity(E,nu,G,typePROBLEM,StrainStressWith4Components) ;
        
        DATA.StrainStressWith4Components = StrainStressWith4Components;
        if isfield(MATLOC,'densCOMP')
            dens = MATLOC.densCOMP ; % Density
            EXIST_DENS  =1;
        end
        
    else
        nameWS =  MATLOC.NAMEWS  ;
        SS =  load(nameWS) ;
        fff =fieldnames(SS) ;
        Celas = [] ;
        for iii = 1:length(fff)
            if strcmp(fff{iii},'Celas') |  strcmp(fff{iii},'celas')
                Celas = SS.(fff{iii}) ;
            end
        end
        if isempty(Celas)
            error('The variable containing the elasticity matrix should be named Celas ')
        end
        
         [celasINV3D Celas] = Compliance_ElasticityGIVEN(Celas,typePROBLEM,StrainStressWith4Components) ; 
        
        %%
        % density
        if isfield(SS,'densCOMP')==1 
            dens = SS.densCOMP ;
            EXIST_DENS = 1 ;
        elseif  isfield(MATLOC,'densCOMP')  
            dens = MATLOC.densCOMP ;
            EXIST_DENS = 1 ;
        else
            %  warning('Density is not defined')
            EXIST_DENS = 0 ;
        end
        
    end
    
    % Angle
    if  ndim == 3
        MATLOC = DefaultField(MATLOC,'ANGLE',0) ;
        a =  MATLOC.ANGLE  ;
        % Transformation of coordinates
        T = RotationMatrix(a) ;
        celasPLY = T*Celas*T' ;
    else
        celasPLY = Celas ;
    end
    
    elemMAT = find(MaterialType == imat) ;
    if ndim ==3
        inv_celasPLY = inv(celasPLY) ;
    else
        inv_celasPLY = 0 ;
    end
    
    for eloc=1:length(elemMAT)
        e = elemMAT(eloc) ;
        celasglo(:,:,e) = celasPLY ;
        
        celasgloINV(:,:,e) = inv_celasPLY ;
        
    end
    if EXIST_DENS==1
        densGLO(elemMAT) = dens ;
    end
    
    %     % Elastic properties (isotropic)
    %     E = MATERIAL.(NAMES_MAT{imat}).E   ; %
    %     nu = MATERIAL.(NAMES_MAT{imat}).nu; % Poisson's coefficient
    %     Gshear = MATERIAL.(NAMES_MAT{imat}).Gshear; % Poisson's coefficient
    %     [celasINV3D celas] = Compliance_Elasticity(E,nu,Gshear,typePROBLEM)  ;
    %
    %     INDEXmat =   MATERIAL.(NAMES_MAT{imat}).INDEX  ;
    %
    %     elemMAT = find(MaterialType == INDEXmat) ;
    %
    %
    %
    
    
    
end


%%%

