function [INDSEL,DATAOUTdistorsion] = PointsToIncludeECM(DATA_GENGAUSS,MESH,Nst,DATA)
% Copy of 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/SVDlibrary/GeneralizedGauss/RestrictedDomainForECMpoints.m
% This function determine the indexes of the Gauss Points that will be
% considered as candidate points in the ECM algorithm
% For cartesian domains, the  domain of interest is determined by
% DATA_GENGAUSS.RESTRICTED_DOMAIN_SELECTION_POINTS
% For general domains, use RESTRICTED_DOMAIN_SELECTION_BYLAYERS >0
% JAHO, 13-JAN-2020
if nargin == 0
    load('tmp1.mat')
end


% Measuring the "DISTORSION" of each element 

 
DATAOUTdistorsion = DistorsionElements(MESH.COOR,MESH.CN,MESH.TypeElement)  ;
 

INDSEL = [] ;
%DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'RESTRICTED_DOMAIN_SELECTION_POINTS',[] ) ;
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'RESTRICTED_DOMAIN_SELECTION_BYLAYERS',0 ) ;
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'WHICH_BOUNDARIES_TO_EXCLUDE',[]) ;

nto = size(Nst,1) ;
ndim = size(MESH.COOR,2) ;
ngausTOT = nto/ndim ;
nelem = size(MESH.CN,1) ;
ngausELEM = ngausTOT/nelem ;
CN = MESH.CN ; 
if DATA_GENGAUSS.RESTRICTED_DOMAIN_SELECTION_BYLAYERS >0   
      %  [NODES_FACES,NODES_LINES] = NodesFacesLinesGID(DATANAMES) ;
        NODES_FACES = MESH.NODES_FACES(:) ;
        NODES_FACES = NODES_FACES(DATA_GENGAUSS.WHICH_BOUNDARIES_TO_EXCLUDE) ;        
        NODESbnd = cell2mat(NODES_FACES) ;
        NODESbnd = NODESbnd(:) ;     
    % Elements sharing one of these nodes
    ElementsSharing = [] ;
    for innode = 1:size(MESH.CN,2)
        IndShared= ismember(MESH.CN(:,innode), NODESbnd) ;
        ElementsSharingLOC =  find(IndShared==1) ;
        ElementsSharing = [ElementsSharing; ElementsSharingLOC] ;
    end
    
    ElementsSharing = unique(ElementsSharing) ;
    % Corresponding INDEX elements
    % ------------------------------ 
    ElementsSharing_previousLAYER = ElementsSharing ;
    % More than 1 layer
    ilayer = 1;
    while  ilayer <= DATA_GENGAUSS.RESTRICTED_DOMAIN_SELECTION_BYLAYERS-1
        % Nodes pertaining to ElementsSharing
        NodesLAYER = MESH.CN(ElementsSharing_previousLAYER,:) ;
        NodesLAYER = unique(NodesLAYER);
        % Elements connected to NodesLAYER
        ElementsSharingNEW = [] ;
        for innode = 1:size(CN,2)
            [dummy,ElementsSharingLOC   ]= intersect(MESH.CN(:,innode), NodesLAYER) ;
            ElementsSharingNEW = [ElementsSharingNEW; ElementsSharingLOC] ;
        end
        
        ElementsSharingNEWlayer = setdiff(ElementsSharingNEW,ElementsSharing) ;
        ElementsSharing_previousLAYER = ElementsSharingNEWlayer ;
        ElementsSharing = [ElementsSharing; ElementsSharingNEWlayer] ;
        ilayer = ilayer+1;
    end
    ElementsExclude = unique(ElementsSharing) ;  
    
    IndGaussExclude = small2large(ElementsExclude,ngausELEM) ;
    
    INDSEL= 1:(ngausELEM*nelem) ;
    INDSEL(IndGaussExclude) = [] ;
else
    ElementsExclude = [] ;
    
end

% METHOD BASED ON DISTORSION
% ***********************************+
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'EXCLUDE_DEGENERATED_ELEMENTS_IN_PERCENTAGE',0) ;

if DATA_GENGAUSS.EXCLUDE_DEGENERATED_ELEMENTS_IN_PERCENTAGE > 0
    AllElements = 1:size(MESH.CN,1) ;
    VariableDEG = DATAOUTdistorsion.distorsionANGLE ;
    [VariableDEGsort,IND1 ]=  sort(VariableDEG,'descend');
    
    PERCENTAGE_EXCLUDE= DATA_GENGAUSS.EXCLUDE_DEGENERATED_ELEMENTS_IN_PERCENTAGE  ;
    NEXCLUDE= PERCENTAGE_EXCLUDE/100*size(MESH.CN,1) ;
    NEXCLUDE = ceil(NEXCLUDE) ;
    disp(['Excluding all elements with distorsion angle higher than (degrees) = ',num2str(VariableDEGsort(NEXCLUDE))]) ;
    
    
    ElementsSharing = IND1(1:NEXCLUDE) ;
    DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'LAYERS_DEGENERATED_ELEMENTS_TO_ELIMINATE',0) ;
    
    ElementsSharing_previousLAYER = ElementsSharing ;
    % More than 1 layer
    ilayer = 1;
    while  ilayer <= DATA_GENGAUSS.LAYERS_DEGENERATED_ELEMENTS_TO_ELIMINATE
        % Nodes pertaining to ElementsSharing
        NodesLAYER = MESH.CN(ElementsSharing_previousLAYER,:) ;
        NodesLAYER = unique(NodesLAYER);
        % Elements connected to NodesLAYER
        ElementsSharingNEW = [] ;
        for innode = 1:size(CN,2)
            [dummy,ElementsSharingLOC   ]= intersect(CN(:,innode), NodesLAYER) ;
            ElementsSharingNEW = [ElementsSharingNEW; ElementsSharingLOC] ;
        end
        
        ElementsSharingNEWlayer = setdiff(ElementsSharingNEW,ElementsSharing) ;
        ElementsSharing_previousLAYER = ElementsSharingNEWlayer ;
        ElementsSharing = [ElementsSharing; ElementsSharingNEWlayer] ;
        ilayer = ilayer+1;
    end
    ElementsExcludeDIST = unique(ElementsSharing) ;
    
    
    %disp(DATAOUTdistorsion.MSG )
    
    ElementsExclude = unique([ElementsExclude;ElementsExcludeDIST]) ;
    
    IndGaussExclude = small2large(ElementsExclude,ngausELEM) ;
    
    INDSEL= 1:(ngausELEM*nelem) ;
    INDSEL(IndGaussExclude) = [] ;
    
else
    DATAOUTdistorsion.MSG = '' ;
end


% CANDIDATE POINTS EXCLUDED FROM INFO IN FILE 
% -------------------------------------------

DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'ListElementsExclude_fromGID',[]) ; 
if ~isempty(DATA_GENGAUSS.ListElementsExclude_fromGID)
     ListElementsToExclude = load(DATA_GENGAUSS.ListElementsExclude_fromGID) ; 
    ListElementsToExclude = ListElementsToExclude(:,1) ; 
    ListGaussToExclude = small2large(ListElementsToExclude,ngausELEM) ; 
    INDSEL  = setdiff(INDSEL,ListGaussToExclude) ; 
end
 

 
