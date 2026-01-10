function [DATAIN] =    ShowSpatialSVDPlot(COORg,PHI_OLD,W,DATALOC,DATAIN,coorECM,z,w )

if nargin == 0
    load('tmp1.mat')
end

ndim  = size(COORg,2) ;

% CARTESIAN MESH
% --------------
xx = cell(1,ndim) ;
TOLloc = 1e-4; % TOLERANCE for determining which points are aligned with those  with y =  ymin
LENGTHSall= zeros(ndim,1); %cellfun

INDICES_ECM_POINTS = zeros(size(coorECM)) ;

for idim = 1:ndim
    xxSORT = unique(sort(COORg(:,idim))) ;
    dxxSORT = diff(xxSORT)  ;
    dx = max(dxxSORT) ;
    dxxSORT = dxxSORT/dx ;
    IROW = find(dxxSORT> TOLloc) ;
    xxLOC = [xxSORT(IROW); xxSORT(end) ];
    xx{idim} = xxLOC;
    LENGTHSall(idim) = length(xxLOC) ;
    INDICES_ECM_POINTS(:,idim) = dsearchn(xxLOC,coorECM(:,idim)) ;
    
end


if prod(LENGTHSall) == size(COORg,1)
    
    % GRID
    
    xINT = cell(1,idim) ;
    if ndim == 2
        [xINT{1}, xINT{2}]= meshgrid(xx{1},xx{2})  ;  % Matrices xINT and yINT for interpolation
        COORsort = [xINT{1}(:),xINT{2}(:)];
    else
        [xINT{1}, xINT{2},xINT{3}]= meshgrid(xx{1},xx{2},xx{3})  ;
        % Notice that   xx{1} --> x direction,  xx{2} --> y direction,  xx{3} --> z direction,
        % Accordingly, xINT{i} is    ny x nx x nz
        COORsort = [] ;
        nnnodesTOT = prod(size(xINT{1})) ;
        for idim = 1:ndim
            COORloc = reshape(xINT{idim},nnnodesTOT,[]) ;
            COORsort =  [COORsort, COORloc] ;
            % Note that, in reshaping, we order the points as follows
        end
    end
    % Coordinates ---> Sorted in a cartesian fashion
    % Now we seek a set of indices INDnear such that COORg(INDnear,:) =
    % COORsort
    dt = delaunayTriangulation(COORg);  % Del. triangulation new coordinates
    INDnear = nearestNeighbor(dt, COORsort);  % Nearest points old coordinates
    [~,INDnearINVERSE ]= sort(INDnear) ;
    zNEW = INDnearINVERSE(z) ;
    %wNEW = w(INDnearINVERSE(z) ) ;
    Wnew = W(INDnear) ;
    %------------------------------------
    % Therefore
    % PHI --> Matrix with the functions we wish to integrate
    PHI =  PHI_OLD(INDnear,:) ; %
    % Check that everything is consistent
    % Exact integral
    EXACT = PHI'*Wnew;
    % Approx. integral
    APPROX = PHI(zNEW,:)'*w;
    ERROR_chech = norm(EXACT-APPROX)/norm(EXACT)*100 ;   % Just for checking
    % DATAIN.MSGPRINT{end+1}  = ['Int. Error after changing indexes = ',num2str(ERROR_chech),'%'] ;
    % disp(DATAIN.MSGPRINT{end})
    
    
    
    % Loop over columns of PHI (one column, one function)
    TOL_SVD = 1e-4 ;
    
    
    
    for iCOL = 1:size(PHI,2)
        PHIloc = PHI(:,iCOL) ;  % Vector to be integrated --> It contains NX x NY x NZ elements
        LABELloc = ['SVD of PHI(:,',num2str(iCOL),')'] ;
            NFIG = 100 + iCOL ;
        if ndim == 2
            
            SpatialSVD_2D(xx,PHIloc,INDICES_ECM_POINTS,TOL_SVD,LABELloc,NFIG) ;           
            
        else
             SpatialSVD_3D(xx,PHIloc,INDICES_ECM_POINTS,TOL_SVD,LABELloc,NFIG) ;        
        end
        
        
    end
    
    
    
    
    
else
    DATAIN.MSGPRINT{end+1} = '-------------------------' ;
    DATAIN.MSGPRINT{end+1} = ['NON-CARTESIAN MESH. SVD spatial decomposition is not possible']
    DATAIN.MSGPRINT{end+1} = '-------------------------' ;
end




