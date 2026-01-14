function  EstimateSingleMode(TOLERANCE_SVD,SCALE_FACTOR,...
    FACTOR_TEST,NAME_MODES,FACTOR_SCALES_ALL,TYPE,TYPE_MODE,imode)


%%%%%%%%%%%55
ModesMatrix = [] ;
COOR = {} ;
LENGTHS = [] ;
FOLDER = ['MODES',filesep] ;
for ifile = 1:length(NAME_MODES)
    load([FOLDER,NAME_MODES{ifile}],'BASES','DATA_REFMESH') ;
   % COOR{ifile} = DATA_REFMESH.COOR ;
    %     xMAX = max(COOR{ifile}(:,1)) ;
    %     xMIN = min(COOR{ifile}(:,1)) ;
    %     LENGTHS(ifile) = xMAX-xMIN ;
    Basis = BASES.(TYPE).U(:,imode) ;
    ModesMatrix = [ModesMatrix,Basis] ;
    
end



DATA.RELATIVE_SVD = 1;

[U,S,V] =RSVDT(ModesMatrix,TOLERANCE_SVD,0,0,DATA ) ;
% ---------------------------------


if  length(SCALE_FACTOR.Z) == 1
    figure(1)
    hold on
    xlabel('Factor scale Y')
    ylabel('Right Singular Value')
    COLORS = {'r*-','go-','bo-','m--','yo-'} ;
    %COLORS = {'r*','g','b','m','y','r*','go','bo','m-','yo'} ;
    
    for i=1:size(V,2)
        h(i) = plot(SCALE_FACTOR.Y ,V(:,i),COLORS{i})
        LL{i} = ['Right S. Vector =',num2str(i)]
    end
    title([TYPE,'; ',TYPE_MODE{imode} ])
    legend(h,LL)
    
    
    %  interpolation
    % -------------------------
    % Mode = sum_i  U_i S_i V_i_interpolated(e)^T
    
    Mode = zeros(size(ModesMatrix,1),1) ;
    for imodeLOC = 1:length(S)
        % Splain interpolation
        Vlocal = V(:,imodeLOC)  ;
        Vtest = spline(LENGTHS,Vlocal,eTEST)  ;
        figure(1)
        hhh(imodeLOC) =  plot(eTEST,Vtest,'Marker','x','MarkerSize',7);
        LLE{imodeLOC} = ['Interpol. R.S.V =',num2str(imodeLOC),' (e=',num2str(eTEST),')']
    end
    legend(hhh,LLE)
    
    
else
    
    %% We have to decompose the right singular vectors
    %  Modes are sorted so that
    % V(:,i) = [V(ez,ey=cte1)      (nz x 1 )
    %           V(ez,ey=cte2,      (nz x 2 )
    %           V(ez,ey =cte3      (nz x 3 )
    % ..............     ]
    % (See arrays FACTOR_SCALES_ALL.Y and FACTOR_SCALES_ALL.Z)
    % Therefore, in the following each column of V is decomposed using the
    % SVD as follows
    % V_i  = U_z*S*V_y
    
    Vall = zeros(size(S)) ;
    
    for i = 1:size(V,2)
        Vloc = V(:,i) ;
        % Now Vloc is reshaped so that it becomes a nZ x nY matrix
        nY = length( SCALE_FACTOR.Y) ;
        nZ = length(  SCALE_FACTOR.Z) ;
        Vloc = reshape(Vloc,nZ,nY)  ;
        [U_z,Sloc,V_y] =RSVDT(Vloc,TOLERANCE_SVD,0,0,DATA ) ;
        
        figure(1000 + i)
        hold on
        subplot(2,1,1)
        hold on
        title(['Global mode =',num2str(i)])
        
        xlabel('Factor scale Y')
        ylabel('Left Singular Vector')
        COLORS =  ColoresMatrix(length(Sloc)) ;% {'r*-','go-','bo-','m--','yo-'} ;
        %COLORS = {'r*','g','b','m','y','r*','go','bo','m-','yo'} ;
        h = [] ;
        LL = {} ;
        
        
        Funinterp_Y = zeros(size(Sloc)) ;
        Funinterp_Z = zeros(size(Sloc)) ;
        
        for iLOC=1:size(V_y,2)
            h(end+1) = plot(SCALE_FACTOR.Y ,V_y(:,iLOC),'Color',COLORS(iLOC,:)) ; 
            LL{end+1} = [' U_i =',num2str(iLOC)]  ;
            
            
            Funinterp_Y(iLOC) = spline(SCALE_FACTOR.Y,V_y(:,iLOC),FACTOR_TEST.Y)  ;
            
            plot(FACTOR_TEST.Y, Funinterp_Y(iLOC),'Marker','x','MarkerSize',7);
            text(FACTOR_TEST.Y, Funinterp_Y(iLOC),'INTP')
            %   LLE{imode} = ['Interpolated   ']
            
            
        end
        %        title([TYPE,'; ',TYPE_MODE{imode} ])
        legend(h,LL)
        
        
        subplot(2,1,2)
        hold on
        xlabel('Factor scale Z')
        ylabel('Right Singular VEctor')
    %    COLORS = {'r*-','go-','bo-','m--','yo-'} ;
        %COLORS = {'r*','g','b','m','y','r*','go','bo','m-','yo'} ;
        h = [] ;
        LL = {} ;
        for iLOC=1:size(U_z,2)
            h(iLOC) = plot(SCALE_FACTOR.Z ,U_z(:,iLOC),'Color',COLORS(iLOC,:))
            LL{iLOC} = [' V_i =',num2str(iLOC)] ;
            
            Funinterp_Z(iLOC) = spline(SCALE_FACTOR.Z,U_z(:,iLOC),FACTOR_TEST.Z)  ;
            
            plot(FACTOR_TEST.Z, Funinterp_Z(iLOC),'Marker','o','MarkerSize',7);
            text(FACTOR_TEST.Z, Funinterp_Z(iLOC),'INTP')
        end
        %      title([TYPE,'; ',TYPE_MODE{imode} ])
        legend(h,LL)
        
        %
        
        Vall(i) = sum((Funinterp_Z.*Funinterp_Y).*Sloc) ;
    end
    
    % Reconstructed mode
    
    MODE_RECONSTRUCTED = zeros(size(U,1),1) ;
    
    for imodeLOC = 1:size(U,2)
        MODE_RECONSTRUCTED = MODE_RECONSTRUCTED + U(:,imodeLOC)*S(imodeLOC)*Vall(imodeLOC) ; 
    end
    
    %% PLOTTING THE MODE IN GID 
    % -------------------------    
    % Reference mesh 
    % ---------------
    ifile = 1;
    load([FOLDER,NAME_MODES{ifile}],'BASES','DATA_REFMESH') ;
    COOR = DATA_REFMESH.COOR ;    
    eyREF =   FACTOR_SCALES_ALL.Y(ifile) ;
    ezREF =    FACTOR_SCALES_ALL.Z(ifile) ;     
    eY = FACTOR_TEST.Y  ;
    eZ = FACTOR_TEST.Z  ;
    COOR(:,2) = COOR(:,2)/eyREF*eY ; 
    COOR(:,3) = COOR(:,3)/ezREF*eZ ;   
    
     NAME_MODES = [cd,filesep,'Mode_',TYPE,'_',TYPE_MODE{imode}];
    DATA= [] ;
    GidPostProcessModes_dom(COOR,DATA_REFMESH.CN,DATA_REFMESH.TypeElement,MODE_RECONSTRUCTED,DATA_REFMESH.posgp,...
        NAME_MODES,DATA,TYPE);
    
    
    
    
end



