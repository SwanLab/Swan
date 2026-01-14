function [PHI,zNEW,COORsort,Wnew,DATAIN] = OrderGaussPointsCartesian(xx,COORg,z,PHI_OLD,W,w,DATAIN)
% Change of indexes--> From the numbering provided by GId, to the numbering
% of the CARTESIAN mesh 


  xINT = cell(1,2) ;
    [xINT{1} xINT{2}]= meshgrid(xx{1},xx{2})  ;  % Matrices xINT and yINT for interpolation
    
    COORsort = [xINT{1}(:),xINT{2}(:)];  % Coordinates ---> Sorted in a cartesian fashion
    % Now we seek a set of indices INDnear such that COORg(INDnear,:) =
    % COORsort
    dt = delaunayTriangulation(COORg(:,1),COORg(:,2));  % Del. triangulation new coordinates
    INDnear = nearestNeighbor(dt, COORsort(:,1),COORsort(:,2));  % Nearest points old coordinates
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
    ERROR_chech = norm(EXACT-APPROX)/norm(EXACT)*100
    DATAIN.MSGPRINT{end+1}  = ['Int. Error after changing indexes = ',num2str(ERROR_chech),'%'] ;
    disp(DATAIN.MSGPRINT{end})