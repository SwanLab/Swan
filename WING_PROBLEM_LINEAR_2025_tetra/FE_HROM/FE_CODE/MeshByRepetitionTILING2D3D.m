function [COORrecons,CNrecons,Materials,dFACES,DOFs_to_INCLUDE] = ...
    MeshByRepetitionTILING2D3D(COORdom,CNrefNEW,NODESfaces,DATAONLINE,MaterialType,DATAINM,nMATdomain,...
    DOMAINS_TO_INCLUDE)

%-------------------------------------------------------------------
% CONSTRUCTING THE MATRIX OF COORDINATES, CONNECTIVITIES AND MATERIAL
% LIBRARY
% -------------------------------------------------------------------
%dbstop('9')
if nargin == 0
    load('tmp.mat')
end

DATAONLINE = DefaultField(DATAONLINE,'NdomZ',1) ;
DATAONLINE = DefaultField(DATAONLINE,'NdomY',1) ;
DATAONLINE = DefaultField(DATAONLINE,'NdomX',1) ;




DATAINM = DefaultField(DATAINM,'TypeUnitCell','HEXAHEDRA') ;


switch  DATAINM.TypeUnitCell
    case  'HEXAHEDRA'
        NdomX = DATAONLINE.NdomX ;
        NdomY = DATAONLINE.NdomY ;
        NdomZ = DATAONLINE.NdomZ ;
        NDOM_dir = {NdomX,NdomY,NdomZ} ;
        f1NOD = cell(1,3) ;f2NOD = cell(1,3) ;
        f1NOD{1} = [NODESfaces{1} ]; % Nodes face xMIN
        f2NOD{1} = [NODESfaces{3} ]; % Nodes face xMAX
        f1NOD{2} = [NODESfaces{2} ]; % Nodes face yMIN
        f2NOD{2} = [NODESfaces{4} ]; % Nodes face yMAX
        f1NOD{3} = [NODESfaces{6} ]; % Nodes face zMIN
        f2NOD{3} = [NODESfaces{5} ]; % Nodes face zMAX
    case 'HEXAG_2D_SQUARE'
        NdomX = DATAONLINE.NdomX ;
        NdomY = DATAONLINE.NdomY ;
        NDOM_dir = {NdomX,NdomY} ;
        f1NOD = cell(1,2) ;f2NOD = cell(1,2) ;
        f1NOD{1} = [NODESfaces{1} ]; % Nodes face xMIN
        f2NOD{1} = [NODESfaces{3} ]; % Nodes face xMAX
        f1NOD{2} = [NODESfaces{2} ]; % Nodes face yMIN
        f2NOD{2} = [NODESfaces{4} ]; % Nodes face yMAX
        
    otherwise
        error('Option  not implemented')
end




%DATAIN.PlotSeparatedDomains = 0.01 ;
DATAINM =  DefaultField(DATAINM,'PlotSeparatedDomains',0) ;

SEPARATION = DATAINM.PlotSeparatedDomains ;
dFACES = zeros(length(NDOM_dir),1) ;
DOFs_to_INCLUDE = [] ;
%%% LOOP OVER DIRECTIONS
if isempty(DOMAINS_TO_INCLUDE)
    COORrecons = [COORdom];
    CNrecons = [CNrefNEW]  ;
    Materials = MaterialType ;
    
    for idir = 1:length(NDOM_dir)
        nDOM = NDOM_dir{idir} ;
        translationVECTOR = COORdom(f2NOD{idir}(1),:)-COORdom(f1NOD{idir}(1),:);
        dFACES(idir) = norm(COORdom(f2NOD{idir}(1),:)-COORdom(f1NOD{idir}(1),:));
        TRANSLATION =0 ;
        %  nmat = length(unique(Materials))/nMATdomain ;
        for e = 2:nDOM
            % Reference point (local)
            TRANSLATION =  TRANSLATION + (1+SEPARATION)*translationVECTOR ;
            COORloc=  COORdom + repmat(TRANSLATION,size(COORdom,1),1) ;
            COORrecons = [COORrecons; COORloc] ;
            CNloc = CNrefNEW +(e-1)*size(COORdom,1) ;
            CNrecons = [CNrecons ; CNloc] ;
            NewMaterial = MaterialType +  (e-1)*nMATdomain ;
            Materials =[Materials; NewMaterial]  ;
            
        end
        
        COORdom = COORrecons ;
        CNrefNEW = CNrecons ;
        MaterialType = Materials ;
        nMATdomain = nMATdomain*nDOM ;
    end
else
    idir = 1 ;
    nDOM = NDOM_dir{idir} ;
    translationVECTOR = COORdom(f2NOD{idir}(1),:)-COORdom(f1NOD{idir}(1),:);
    dFACES(idir) = norm(COORdom(f2NOD{idir}(1),:)-COORdom(f1NOD{idir}(1),:));
    TRANSLATION =0 ;
    
    COORrecons = [] ; %[COORdom];
    CNrecons = [] ;%[CNrefNEW]  ;
    Materials = [] ;% MaterialType ;
    DOFs_to_INCLUDE = [] ;
    nDOFS = prod(size(COORdom)) ;
    DOFS = (1:nDOFS)' ;
    %  nmat = length(unique(Materials))/nMATdomain ;
    for I = 1:length(DOMAINS_TO_INCLUDE)
        e = DOMAINS_TO_INCLUDE(I) ;
        NewDOFS = DOFS + (e-1)*nDOFS ;
        DOFs_to_INCLUDE =[DOFs_to_INCLUDE ; NewDOFS] ;
        % Reference point (local)
        TRANSLATION =  (e-1)*translationVECTOR ;
        COORloc=  COORdom + repmat(TRANSLATION,size(COORdom,1),1) ;
        COORrecons = [COORrecons; COORloc] ;
        CNloc = CNrefNEW +(I-1)*size(COORdom,1) ;
        CNrecons = [CNrecons ; CNloc] ;
        NewMaterial = MaterialType +  (e-1)*nMATdomain ;
        Materials =[Materials; NewMaterial]  ;
        
    end
    
    COORdom = COORrecons ;
    CNrefNEW = CNrecons ;
    MaterialType = Materials ;
    %nMATdomain = nMATdomain*nDOM ;
    
    
    
end