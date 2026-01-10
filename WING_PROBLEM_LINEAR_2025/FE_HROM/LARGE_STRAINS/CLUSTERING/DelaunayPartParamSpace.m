function [INDEX_CLUSTERS_WITH_OVERLAPPING,INDEX_CLUSTERS_NO_OVERLAPPING,DELAUNAY_PARAM_SPACE] =...
    DelaunayPartParamSpace(DATAoffline,MATRIX_POINTS_SPACE_PARAMETER) ;

if nargin == 0
    load('tmp1.mat')
    DATAoffline.PARTITION_PARAMETER{1} = linspace(0,1,5) ;
    DATAoffline.PARTITION_PARAMETER{2} = linspace(0,1,5) ;
end
% Limits cartesian domain
xMIN = min(MATRIX_POINTS_SPACE_PARAMETER) ;
xMAX = max(MATRIX_POINTS_SPACE_PARAMETER) ;

% Normalization points parameter domain
Xnorm = zeros(size(MATRIX_POINTS_SPACE_PARAMETER)) ;
for idim = 1:length(xMIN)
    if  xMAX(idim)-xMIN(idim) == 0
        error('Program this case')
    else
        Xnorm(:,idim) = (MATRIX_POINTS_SPACE_PARAMETER(:,idim)-xMIN(idim))/(xMAX(idim)-xMIN(idim)) ;
    end
end

DATAoffline = DefaultField(DATAoffline,'REGULAR_MESH_Delaunay',0) ;
DATAoffline = DefaultField(DATAoffline,'TrainingPointsAreVertex',0) ;


% 1) Creating a regular grid
%*****************************+
if DATAoffline.TrainingPointsAreVertex == 0
    ninp = size(Xnorm,2) ;
    DATAoffline = DefaultField(DATAoffline,'PARTITION_PARAMETER',[]) ;
    if  isempty(DATAoffline.PARTITION_PARAMETER)
        nclusters = DATAoffline.NCLUSTERS_BASIS_DISP ;
        ndiv = ceil((nclusters)^(1/ninp))*ones(1,ninp) ;
        Xinp = cell(1,ninp) ;
        for iparaminp = 1:ninp
            Xinp{iparaminp} = linspace(0,1,ndiv) ;
        end
    else
        Xinp = cell(1,ninp) ;
        for iparaminp = 1:ninp
            Xinp{iparaminp} =   DATAoffline.PARTITION_PARAMETER{iparaminp}  ;
        end
        [aaaa,ndiv] = cellfun(@size,Xinp)  ;
    end
    
    XGRID = cell(size(Xinp)) ;
    [XGRID{:}]  = ndgrid(Xinp{:}) ;
    % Coordinates
    % ************
    COOR = cell(1,ninp) ;
    for iparaminp = 1:ninp
        COOR{iparaminp} = XGRID{iparaminp}(:) ;
    end
    COOR = cell2mat(COOR) ;
    if ninp == 2
        figure(127)
        hold on
        xlabel('Param 1') ;
        ylabel('Paaram 2') ;
        % plot(COOR(:,1),COOR(:,2),'k*') ;
    end
    
    
else
    COOR = [];
end

% **********************************************************

%2) TRAINING PARAMETERS CLOSEST TO COOR
% ************************************+

% for iii = 1:size(COOR,1)
%     [Idx,D] = knnsearch(MATRIX_POINTS_SPACE_PARAMETER,COOR(iii,:))  ;
%     dist = bsxfun(@minus,MATRIX_POINTS_SPACE_PARAMETER',COOR(iii,:)')' ;
%     dist2  =sum(dist.^2,2);
%     [aaa,bbb]  =min(dist2) ;
%     text(COOR(iii,1),COOR(iii,2),num2str(iii))
%     plot(MATRIX_POINTS_SPACE_PARAMETER(Idx,1),MATRIX_POINTS_SPACE_PARAMETER(Idx,2),'ro')
% end


if DATAoffline.REGULAR_MESH_Delaunay == 0
    if ~isempty(COOR)
        [Idx,D] = knnsearch(Xnorm,COOR)  ;
        
        Idx = unique(Idx) ;
        COORnew = Xnorm(Idx,:) ;
    else
        
        DATAoffline = DefaultField(DATAoffline,'FractionPointsPointsAreVertex',1);
        
        if DATAoffline.FractionPointsPointsAreVertex == 1
            [COORnew,IND_AAA,IND_BBB ]= unique(Xnorm,'rows') ;
        else
            [COORnew,IND_AAA_1,IND_BBB ]= unique(Xnorm,'rows') ;
            rng default
            npointsALL = size(COORnew,1) ;
            npontsSELECT = ceil(DATAoffline.FractionPointsPointsAreVertex*npointsALL) ;
            IND_AAA_2 =  randi(npointsALL,[1,npontsSELECT]) ;
            
            IND_AAA_2 = unique(IND_AAA_2) ;
            
            COORnew = COORnew(IND_AAA_2,:) ;
            IND_AAA = [] ;
        end
    end
else
    
    COORnew = COOR ;
    
end

% if ninp == 2
%   %  plot(COORnew(:,1),COORnew(:,2),'rx') ;
% end

% for iii = 1:size(MATRIX_POINTS_SPACE_PARAMETER,1)
%     plot(MATRIX_POINTS_SPACE_PARAMETER(iii,1),MATRIX_POINTS_SPACE_PARAMETER(iii,2),'y.')
%     text(MATRIX_POINTS_SPACE_PARAMETER(iii,1),MATRIX_POINTS_SPACE_PARAMETER(iii,2),num2str(iii)) ;
% end




% 3) Construct a Delaunay triangulation of COORnew
DATAoffline = DefaultField(DATAoffline,'IS1DparamSPACE',0);
%[UU,VV,SS] = SVDT(COORnew') ;
if  DATAoffline.IS1DparamSPACE == 1
    % 1D parametric space
    
    nCOOR = sum(COORnew.^2,2);
    [~,IDX] = sort(nCOOR) ;
    CN = [IDX(1:end-1),IDX(2:end)] ;
else
    CN =  delaunayn(COORnew) ;
    %DT = delaunayTriangulation(COORnew) ;
    
    %CN = DT.ConnectivityList ;
    %
    
end


% Centroid
CENTROID = zeros(size(CN,1),size(COORnew,2)) ;
for inode = 1:size(CN,2)
    CENTROID = CENTROID +  COORnew(CN(:,inode),:)/size(CN,2) ;
end


if size(CN,2) == 3
    if size(COORnew,2) == 2
        triplot(CN,COORnew(:,1),COORnew(:,2));
    end
    
    
    %     for ielem = 1:size(CN,1)
    %         text(CENTROID(ielem,1),CENTROID(ielem,2),num2str(ielem)) ;
    %     end
    %
    %     for inode = 1:size(COORnew,1)
    %         dd= text(COORnew(inode,1),COORnew(inode,2),num2str(inode)) ;
    %         set(dd,'FontSize',8,'Color',[1,0,0]) ;
    %     end
    
end


% SNAPSHOTS ASSOCIATED TO EACH CLUSTER
% ************************************
nclusters = size(CN,1) ;
% nnodeE = size(CN,2) ;


INCLUDE_ALL_NEIG =0 ;
if INCLUDE_ALL_NEIG == 1
    [InvCNmatrix, ElemNode, NeighboringElements, ElemShared]= InverseConnectMatrix(CN) ;
else
    
    % This should be done in a vectorized fashion
    [InvCNmatrix, ElemNode, NeighboringElements, ElemShared]= InverseConnectMatrix(CN) ;
    if size(COORnew,2) == 2
        
        if size(CN,2) ==3
            NeighboringElements = NeigSharingEdgeDELAN(CN,InvCNmatrix) ;
        elseif size(CN,2) ==2
            NeighboringElements = NeigSharingEdgeDELAN_1d(CN,InvCNmatrix) ;
            
        else
            error('Option not implemented')
        end
        
        
        
    else
        
        error('Not implemented')
    end
    
    
end



INDEX_CLUSTERS_NO_OVERLAPPING = cell(1,nclusters) ;
INDEX_CLUSTERS_WITH_OVERLAPPING = cell(1,nclusters) ;

%for icluster = 1:nclusters
if ~isempty(COOR)
    IND_clusters = tsearchn(COORnew, CN , Xnorm);
    
    
    for icluster = 1:nclusters
        IndLocal = find(IND_clusters == icluster)  ;
        
        
        if DATAoffline.REGULAR_MESH_Delaunay == 0
            IndLocal = [IndLocal;CN(icluster,:)' ] ;
        end
        
        INDEX_CLUSTERS_NO_OVERLAPPING{icluster}  = unique(IndLocal);
    end
    
else
    % Points are vertices
    for  icluster = 1:nclusters
        IND_local = CN(icluster,:) ;
        if ~isempty(IND_AAA)
            INDEX_CLUSTERS_NO_OVERLAPPING{icluster}  =  IND_AAA(IND_local)  ;
        else
            INDEX_CLUSTERS_NO_OVERLAPPING{icluster}  =  IND_AAA_1(IND_AAA_2(IND_local))  ;
        end
        
    end
    
end
% NEIGHBORS

for icluster = 1:nclusters
    IndLocal = INDEX_CLUSTERS_NO_OVERLAPPING{icluster}  ;
    CLUSTER_NEIG = NeighboringElements(icluster,:) ;
    for jclusterLOC = 1:length(CLUSTER_NEIG)
        jcluster =CLUSTER_NEIG(jclusterLOC) ;
        if jcluster ~=0
            IndLocal = [IndLocal;INDEX_CLUSTERS_NO_OVERLAPPING{jcluster}] ;
        end
        
    end
    INDEX_CLUSTERS_WITH_OVERLAPPING{icluster} = unique(IndLocal)  ;
    
end


[nsnap_nooverlap, ~] = cellfun(@size,INDEX_CLUSTERS_NO_OVERLAPPING) ;

[nsnap_withoverla, ~] = cellfun(@size,INDEX_CLUSTERS_WITH_OVERLAPPING) ;



% Reverse the normalization centroids
%CENTROIDnew = zeros(size(CENTROID)) ;
% POINTS_SIMPLEX = zeros(size(COORnew)) ;
%
% for idim = 1:length(xMIN)
%     CENTROIDnew(:,idim) =     xMIN(idim) +  CENTROID(:,idim)*(xMAX(idim)-xMIN(idim)) ;
%     POINTS_SIMPLEX(:,idim) =     xMIN(idim) +  COORnew(:,idim)*(xMAX(idim)-xMIN(idim)) ;
%
% end

DELAUNAY_PARAM_SPACE.COOR =  COORnew ;
DELAUNAY_PARAM_SPACE.CN =  CN ;
DELAUNAY_PARAM_SPACE.InvCNmatrix = InvCNmatrix ;
DELAUNAY_PARAM_SPACE.BOUNDARIES_MAX = xMAX ;
DELAUNAY_PARAM_SPACE.BOUNDARIES_MIN= xMIN ;
DELAUNAY_PARAM_SPACE.CENTROID= CENTROID ;



%
% if ninp == 2
%
%     %CENTROID = (COORnew(CN(:,1),:) +COORnew(CN(:,2),:)  + COORnew(CN(:,3),:))/3;
%     %     iclusterPLOT = 720 ;
%     %     Cl2Plot = [iclusterPLOT,NeighboringElements(iclusterPLOT,:)] ;
%     %     for i = 1:length(Cl2Plot)
%     %         icluster = Cl2Plot(i) ;
%     %         text(CENTROID(icluster,1),CENTROID(icluster,2),num2str(icluster))
%     %     end
%
% end



%end




%



