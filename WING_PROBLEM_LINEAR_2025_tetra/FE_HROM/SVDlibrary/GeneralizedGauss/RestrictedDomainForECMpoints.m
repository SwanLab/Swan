function [INDSEL,DATAOUTdistorsion] = RestrictedDomainForECMpoints(DATA_GENGAUSS,COOR,Nst,CN,CONNECTb)
% This function determine the indexes of the Gauss Points that will be
% consdiered as candidate points in the ECM algorithm
% For cartesian domains, the  domain of interest is determined by
% DATA_GENGAUSS.RESTRICTED_DOMAIN_SELECTION_POINTS
% For general domains, use RESTRICTED_DOMAIN_SELECTION_BYLAYERS >0
% JAHO, 12-APril-2020(30th day of COVID-19 confinement )
if nargin == 0
    load('tmp1.mat')
end

DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'TypeElement',[]) ;
if ~isempty(DATA_GENGAUSS.TypeElement)
    DATAOUTdistorsion = DistorsionElements(COOR,CN,DATA_GENGAUSS.TypeElement)  ;
end

INDSEL = [] ;
%DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'RESTRICTED_DOMAIN_SELECTION_POINTS',[] ) ;
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'RESTRICTED_DOMAIN_SELECTION_BYLAYERS',0 ) ;
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'WHICH_BOUNDARIES_TO_EXCLUDE',[]) ;



% if DATA_GENGAUSS.RESTRICTED_DOMAIN_SELECTION_BYLAYERS >0
%     DATA_GENGAUSS.RESTRICTED_DOMAIN_SELECTION_POINTS = [];
% end

% if ~isempty(DATA_GENGAUSS.RESTRICTED_DOMAIN_SELECTION_POINTS)
%     % CARTESIAN DOMAIN (obsolete)
%     error('Obsolete option...')
%     COORg = COOR' ;
%     xGAUSSall = Nst*COORg(:) ;
%
%     COORg = reshape(xGAUSSall,ndim,[])' ;
%
%     ndim = 2;
%     for idim = 1:ndim
%         DOM =  DATA_GENGAUSS.RESTRICTED_DOMAIN_SELECTION_POINTS{idim} ;
%         INDMIN = find((COORg(:,idim)-DOM(1))>=0) ;
%         INDMAX = find((COORg(:,idim)-DOM(2))<=0) ;
%         INDSEL{idim} = intersect(INDMIN,INDMAX) ;
%     end
%     INDSEL = intersect(INDSEL{1},INDSEL{2});
% end

nto = size(Nst,1) ;
ndim = size(COOR,2) ;
ngausTOT = nto/ndim ;
nelem = size(CN,1) ;
ngausELEM = ngausTOT/nelem ;
if DATA_GENGAUSS.RESTRICTED_DOMAIN_SELECTION_BYLAYERS >0
    
    if iscell(CONNECTb)
        % Domain decomposition problems
        if   ~isempty(DATA_GENGAUSS.WHICH_BOUNDARIES_TO_EXCLUDE)  %&& DATA_GENGAUSS.WHICH_BOUNDARIES_TO_EXCLUDE ~=0
            CONNECTb = CONNECTb(DATA_GENGAUSS.WHICH_BOUNDARIES_TO_EXCLUDE) ;
        end
        CONNECTb = cell2mat(CONNECTb');
        % Boundary nodes
        NODESbnd = unique(CONNECTb(:)) ;
    else
        DATANAMES = DATA_GENGAUSS.NameFileMesh(1:end-4);
        
        [NODES_FACES,NODES_LINES] = NodesFacesLinesGID(DATANAMES) ;
        NODES_FACES = NODES_FACES(:) ;
        NODES_FACES = NODES_FACES(DATA_GENGAUSS.WHICH_BOUNDARIES_TO_EXCLUDE) ;
        
        NODESbnd = cell2mat(NODES_FACES) ;
        NODESbnd = NODESbnd(:) ;
        % We have to read the nodes of each surface/line
        
    end
    
    %   NODESbnd = [NODESbndINI;NODESbnd ];
    
    
    
    
    
    % Elements sharing one of these nodes
    ElementsSharing = [] ;
    for innode = 1:size(CN,2)
        IndShared= ismember(CN(:,innode), NODESbnd) ;
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
        NodesLAYER = CN(ElementsSharing_previousLAYER,:) ;
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
    ElementsExclude = unique(ElementsSharing) ;
    
    
    
    
    
    IndGaussExclude = small2large(ElementsExclude,ngausELEM) ;
    
    INDSEL= 1:(ngausELEM*nelem) ;
    INDSEL(IndGaussExclude) = [] ;
else
    ElementsExclude = [] ;
    
end

DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'EXCLUDE_DEGENERATED_ELEMENTS_IN_PERCENTAGE',0) ;

if DATA_GENGAUSS.EXCLUDE_DEGENERATED_ELEMENTS_IN_PERCENTAGE > 0
    AllElements = 1:size(CN,1) ;
    VariableDEG = DATAOUTdistorsion.distorsionANGLE ;
    [VariableDEGsort,IND1 ]=  sort(VariableDEG,'descend');
    
    PERCENTAGE_EXCLUDE= DATA_GENGAUSS.EXCLUDE_DEGENERATED_ELEMENTS_IN_PERCENTAGE  ;
    NEXCLUDE= PERCENTAGE_EXCLUDE/100*size(CN,1) ;
    NEXCLUDE = ceil(NEXCLUDE) ;
    DATAOUTdistorsion.MSG = ['Excluding all elements with distorsion angle higher than (degrees) = ',num2str(VariableDEGsort(NEXCLUDE))] ;
    
    
    ElementsSharing = IND1(1:NEXCLUDE) ;
    DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'LAYERS_DEGENERATED_ELEMENTS_TO_ELIMINATE',0) ;
    
    ElementsSharing_previousLAYER = ElementsSharing ;
    % More than 1 layer
    ilayer = 1;
    while  ilayer <= DATA_GENGAUSS.LAYERS_DEGENERATED_ELEMENTS_TO_ELIMINATE
        % Nodes pertaining to ElementsSharing
        NodesLAYER = CN(ElementsSharing_previousLAYER,:) ;
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
    
    
    disp(DATAOUTdistorsion.MSG )
    
    ElementsExclude = unique([ElementsExclude;ElementsExcludeDIST]) ;
    
    IndGaussExclude = small2large(ElementsExclude,ngausELEM) ;
    
    INDSEL= 1:(ngausELEM*nelem) ;
    INDSEL(IndGaussExclude) = [] ;
    
else
    DATAOUTdistorsion.MSG = '' ;
end
