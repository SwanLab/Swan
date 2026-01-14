function GAMMA_INTERP = InterpolationMatricesDirect...
    (FACTOR_SCALES_ALL,GAMMA,FACTOR_TEST,SCALE_FACTOR,DATAINPUT)

if nargin == 0
    load('tmp.mat')
end

DATAINPUT = DefaultField(DATAINPUT,'SHOW_PLOTS_SVD',1) ;

SHOW = DATAINPUT.SHOW_PLOTS_SVD ;

% S1) Determining physical parameters (maximum 2 )
% ------------------------------------------------
DATAINPUT = DefaultField(DATAINPUT,'PARAM_INTERP',[]) ;
if isempty(DATAINPUT.PARAM_INTERP)
    PARAM = {'X','Y','Z'} ;
    IPARAM = [] ;
    FACTOR_SCALES_ALL = DefaultField(FACTOR_SCALES_ALL,'X',ones(size(FACTOR_SCALES_ALL.Y))) ;
    for iparam = 1:length(PARAM)
        DDD = diff(FACTOR_SCALES_ALL.(PARAM{iparam}));
        if any(DDD)
            IPARAM(end+1) = iparam ;
        end
    end
    
    PARAM = PARAM(IPARAM) ;
    
else
    PARAM = {'parameter'} ;
    SCALE_FACTOR.(PARAM{1}) =  DATAINPUT.PARAM_INTERP{1} ;
    FACTOR_TEST_p.(PARAM{1})  = FACTOR_TEST ;
    FACTOR_TEST = FACTOR_TEST_p ;
end

% -------------
% Reshape GAMMA
%--------------
[nrows, ncols] = size(GAMMA{1}) ;
for iproj = 1:length(GAMMA)
    GAMMA{iproj} = GAMMA{iproj}(:) ;
end
GAMMA = cell2mat(GAMMA) ;
% SVD GAMMA ---for separation of variables
DATA.RELATIVE_SVD = 1;
DATAINPUT = DefaultField(DATAINPUT,'TOLERANCE_SVD',1e-4) ;
TOLERANCE_SVD = DATAINPUT.TOLERANCE_SVD;
[U,S,V] =RSVDT(GAMMA,TOLERANCE_SVD,0,0,DATA ) ;

%-------------------------
normMATRIX = norm(GAMMA,'fro') ;


DATAINPUT = DefaultField(DATAINPUT,'TITLE',[]) ; 
DATAINPUT = DefaultField(DATAINPUT,'NFIGURE',1) ; 


if  length(PARAM) == 1
    if SHOW == 1
        figure(DATAINPUT.NFIGURE)
        hold on
        xlabel(['Factor scale ',PARAM{1}])
        ylabel('Right Singular Value')
        if ~isempty(DATAINPUT.TITLE)
            title(DATAINPUT.TITLE)
        end
        
        COLORS =  ColoresMatrix(length(S)) ;
        for i=1:size(V,2)
            h(i) = plot(SCALE_FACTOR.(PARAM{1}) ,V(:,i),'Color',COLORS(i,:));
            LL{i} = ['Right S. Vector =',num2str(i),' (S_{rel}=',num2str(S(i)/S(1),4),')'];
        end
        %    title([TYPE,'; ',TYPE_MODE{imode} ])
        legend(h,LL);
        
    end
    %  interpolation
    % -------------------------
    % Mode = sum_i  U_i S_i V_i_interpolated(e)^T
    eTEST = FACTOR_TEST.(PARAM{1}) ;
    GAMMA_INTERP = zeros(size(U,1),1) ;
    for imodeLOC = 1:length(S)
        % Splain interpolation
        Vlocal = V(:,imodeLOC)  ;
        Vtest = spline(SCALE_FACTOR.(PARAM{1}) ,Vlocal,eTEST)  ;
        if SHOW == 1
            figure(1)
            hhh(imodeLOC) =  plot(eTEST,Vtest,'Marker','x','MarkerSize',7);
            %   LLE{imodeLOC} = ['Interpol. R.S.V =',num2str(imodeLOC),' (e=',num2str(eTEST),')'] ;
        end
        GAMMA_INTERP = GAMMA_INTERP + U(:,imodeLOC)*S(imodeLOC)*Vtest ;
    end
    
    
    
    %  legend(hhh,LLE);
    
    
    GAMMA_INTERP = reshape(GAMMA_INTERP,nrows,ncols) ;
    
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
        nY = length( SCALE_FACTOR.(PARAM{1})) ;
        nZ = length(  SCALE_FACTOR.(PARAM{2})) ;
        Vloc = reshape(Vloc,nZ,nY)  ;
        [U_z,Sloc,V_y] =RSVDT(Vloc,TOLERANCE_SVD,0,0,DATA ) ;
        
        if SHOW == 1
            figure(1000 + i)
            hold on
            subplot(2,1,1)
            hold on
            title(['Global mode =',num2str(i),'  (SVALrel = ',  num2str(S(i)/S(1)),' )'])
            
            xlabel('Factor scale Y')
            ylabel('Left Singular Vector')
            COLORS =  ColoresMatrix(length(Sloc)) ;% {'r*-','go-','bo-','m--','yo-'} ;
            %COLORS = {'r*','g','b','m','y','r*','go','bo','m-','yo'} ;
            h = [] ;
            LL = {} ;
        end
        
        Funinterp_Y = zeros(size(Sloc)) ;
        Funinterp_Z = zeros(size(Sloc)) ;
        
        for iLOC=1:size(V_y,2)
            if SHOW == 1
                h(end+1) = plot(SCALE_FACTOR.Y ,V_y(:,iLOC),'Color',COLORS(iLOC,:)) ;
                LL{end+1} = [' U_i =',num2str(iLOC)]  ;
            end
            
            Funinterp_Y(iLOC) = spline(SCALE_FACTOR.Y,V_y(:,iLOC),FACTOR_TEST.Y)  ;
            if SHOW == 1
                plot(FACTOR_TEST.Y, Funinterp_Y(iLOC),'Marker','x','MarkerSize',7);
                text(FACTOR_TEST.Y, Funinterp_Y(iLOC),'INTP') ;
            end
            %   LLE{imode} = ['Interpolated   ']
            
            
        end
        %        title([TYPE,'; ',TYPE_MODE{imode} ])
        if SHOW == 1
            legend(h,LL)
        end
        
        if SHOW == 1
            subplot(2,1,2)
            hold on
            xlabel('Factor scale Z')
            ylabel('Right Singular VEctor')
            %    COLORS = {'r*-','go-','bo-','m--','yo-'} ;
            %COLORS = {'r*','g','b','m','y','r*','go','bo','m-','yo'} ;
            h = [] ;
            LL = {} ;
        end
        for iLOC=1:size(U_z,2)
            if SHOW == 1
                h(iLOC) = plot(SCALE_FACTOR.Z ,U_z(:,iLOC),'Color',COLORS(iLOC,:))
                LL{iLOC} = [' V_i =',num2str(iLOC)] ;
            end
            
            Funinterp_Z(iLOC) = spline(SCALE_FACTOR.Z,U_z(:,iLOC),FACTOR_TEST.Z)  ;
            if SHOW == 1
                plot(FACTOR_TEST.Z, Funinterp_Z(iLOC),'Marker','o','MarkerSize',7);
                text(FACTOR_TEST.Z, Funinterp_Z(iLOC),'INTP') ;
            end
        end
        %      title([TYPE,'; ',TYPE_MODE{imode} ])
        if SHOW == 1
            legend(h,LL)
        end
        
        %
        
        Vall(i) = sum((Funinterp_Z.*Funinterp_Y).*Sloc) ;
    end
    
    % Reconstructed mode
    
    MODE_RECONSTRUCTED = zeros(size(U,1),1) ;
    
    for imodeLOC = 1:size(U,2)
        MODE_RECONSTRUCTED = MODE_RECONSTRUCTED + U(:,imodeLOC)*S(imodeLOC)*Vall(imodeLOC) ;
    end
    
    
    GAMMA_INTERP = reshape(MODE_RECONSTRUCTED,nrows,ncols) ;
    
    %     %% PLOTTING THE MODE IN GID
    %     % -------------------------
    %     % Reference mesh
    %     % ---------------
    %     ifile = 1;
    %     load([FOLDER,NAME_MODES{ifile}],'BASES','DATA_REFMESH') ;
    %     COOR = DATA_REFMESH.COOR ;
    %     eyREF =   FACTOR_SCALES_ALL.Y(ifile) ;
    %     ezREF =    FACTOR_SCALES_ALL.Z(ifile) ;
    %     eY = FACTOR_TEST.Y  ;
    %     eZ = FACTOR_TEST.Z  ;
    %     COOR(:,2) = COOR(:,2)/eyREF*eY ;
    %     COOR(:,3) = COOR(:,3)/ezREF*eZ ;
    %
    %     ndof = prod(size(COOR));
    %     MODE_RECONSTRUCTED = reshape(MODE_RECONSTRUCTED,ndof,[]) ;
    %
    %     if  DATAINPUT.ESTIMATE_ALL_MODES_AT_ONCE ==0
    %         TYPE_MODE = {TYPE_MODE{DATAINPUT.imode}} ;
    %
    %     end
    %
    %     for imode = 1:length(TYPE_MODE)
    %
    %         NAME_MODES = [cd,filesep,'Mode_',TYPE,'_',TYPE_MODE{imode}];
    %         DATA= [] ;
    %
    %
    %         MODES_LOC =    MODE_RECONSTRUCTED(:,imode) ;
    %
    %
    %
    %         GidPostProcessModes_dom(COOR,DATA_REFMESH.CN,DATA_REFMESH.TypeElement,...
    %             MODES_LOC,DATA_REFMESH.posgp,...
    %             NAME_MODES,DATA,TYPE);
    %
    %     end
    %
    
    
    
end

%
%
% %%% Comparison with the modes  obtained by running the FE program
% % ----------------------------
% RUN_MODES.NAME_MODES_FILES = DATAINPUT.NAME_WS_TEST;
% FOLDER = [cd,filesep,'MODES',filesep] ;
%
% NAME_MODES   =[FOLDER,'MODES_',RUN_MODES.NAME_MODES_FILES,'.mat'] ;
%
% if ~exist(NAME_MODES ,'file')
%     SCALE_FACTOR_loc = FACTOR_TEST  ;
%     SCALE_FACTOR_loc.X = SCALE_FACTOR.X ;
%
%     RUN_FE.SCALE_FACTOR = SCALE_FACTOR_loc ;
%     LSCRIPT_fun(DATAINPUT.EXECUTABLE_FOLDER,DATAINPUT.NAME_GEOMETRY_MATERIAL_DATA,...
%         DATAINPUT.NAME_LOAD_DATA,RUN_MODES,RUN_FE) ;
% end
%
% load(NAME_MODES,'BASES','DATA_REFMESH') ;
% ModesMatrixFE = BASES.(TYPE).U ;
%
% if  DATAINPUT.ESTIMATE_ALL_MODES_AT_ONCE ==0
%     ModesMatrixFE = ModesMatrixFE(:,DATAINPUT.imode) ;
% end
%
% ERROR_APPROXIMATION = []; % Normalized vector
%
% for imodeLOC = 1:size(ModesMatrixFE,2)
%     V_FE  = ModesMatrixFE(:,imodeLOC) ;
%     V_FE = V_FE/norm(V_FE) ;
%     V_RECONSTR = MODE_RECONSTRUCTED(:,imodeLOC) ;
%     V_RECONSTR = V_RECONSTR/norm(V_RECONSTR) ;
%
%     ERROR_APPROXIMATION(imodeLOC) = norm(V_RECONSTR-V_FE)*100 ;
%
% end
%
% disp(['Error approx. (= %)',num2str(ERROR_APPROXIMATION)])
%
% if size(ModesMatrixFE,2) ==1
%     ERROR_ALL = V_RECONSTR-V_FE;
%     NAME_MODES = [cd,filesep,'Error_',TYPE,'_',TYPE_MODE{imode}];
%     GidPostProcessModes_dom(COOR,DATA_REFMESH.CN,DATA_REFMESH.TypeElement,...
%         ERROR_ALL,DATA_REFMESH.posgp,...
%         NAME_MODES,DATA,TYPE);
%
% end
%
