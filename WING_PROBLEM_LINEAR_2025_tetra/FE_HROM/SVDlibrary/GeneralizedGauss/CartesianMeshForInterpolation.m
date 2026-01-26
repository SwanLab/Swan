function xx = CartesianMeshForInterpolation(COORg,DATALOC,z,xxOLD)




ndim = size(COORg,2) ;



xx  =cell(1,ndim) ;

DATALOC = DefaultField(DATALOC,'SIZE_CARTESIAN_MESH',[0.05,0.05]) ;

DATALOC = DefaultField(DATALOC,'CARTESIAN_MESH_COINCIDENT_WITH_ECM_POINTS',0) ;

DATALOC = DefaultField(DATALOC,'IS_CARTESIAN_WITH_HOLES',0) ;


for idim = 1:ndim
    xMAX =max(COORg(:,idim)) ;
    xMIN =min(COORg(:,idim)) ;
    Lx = xMAX-xMIN ;
        dx  = DATALOC.SIZE_CARTESIAN_MESH(idim) ;
    
    if  DATALOC.IS_CARTESIAN_WITH_HOLES == 1
        
       xx =  xxOLD ; 
        
    elseif   DATALOC.CARTESIAN_MESH_COINCIDENT_WITH_ECM_POINTS == 0
        % Equally spaced points 
        % ---------------------
             nsteps = ceil(Lx/dx);
             xx{idim} = linspace(xMIN,xMAX,nsteps);  % Equally spaced
    else
        
        
       
        xECM = sort(unique(COORg(z,idim))) ;
        xECM = [xMIN,xECM',xMAX] ;
        xTOTAL = cell(1,length(xECM)-1) ;
        
        for itramo = 1:length(xTOTAL)
            xINI = xECM(itramo) ;
            xFIN = xECM(itramo+1) ;
            dxLOC = xFIN-xINI ;
            nsteps =     (dxLOC/dx) ;
            if nsteps <=1.1
                xLOC = [xINI,xFIN] ;
            else
                nsteps = ceil(nsteps) ;
                xLOC = linspace(xINI,xFIN,nsteps);
            end
            xTOTAL{itramo} = xLOC;
        end
        
        xTOTAL = cell2mat(xTOTAL) ;
        xx{idim} = unique(xTOTAL);
        % Now we add the coordinates of the ECM points
        
    end
    
    
end


