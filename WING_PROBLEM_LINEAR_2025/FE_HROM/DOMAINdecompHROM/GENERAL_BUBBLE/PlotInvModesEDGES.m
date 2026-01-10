function PlotInvModesEDGES(DATA,PhiDEF_b,Mintf)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/09_BACK_to_HEXAG.mlx
if nargin == 0
    load('tmp.mat')
end

INFO_EDGES = DATA.INFO_EDGES;

INFO_EDGES = DefaultField(INFO_EDGES,'xiEDGES',[]) ;
ndim = DATA.MESH.ndim;
if ndim ~=2
    error('THis is only valid for 2D')
end


% STEP 1) WEIGHTED SVD OF PhiDEF_b
TOL_ALL = 1e-1;
DATAlocc.TOL =TOL_ALL;
[PhiDEF_b,Sbs,Vbs,Mintf_CHOL] = WSVDT( PhiDEF_b,Mintf,DATAlocc) ;

REQUIRED_degree_all_edges = [] ; 

if  ~isempty(INFO_EDGES.xiEDGES)
    
    REQUIRED_degree_all_edges = zeros(length(INFO_EDGES.xiEDGES),1);
    for iedge = 1:length(INFO_EDGES.xiEDGES)
        xi = INFO_EDGES.xiEDGES{iedge} ;
        
        
        
        
        figure(981+iedge)
        hold on
        
        IND_NODES= INFO_EDGES.IndPointsBNDedge{iedge} ;
        
        IND_DOFS = small2large(IND_NODES,ndim) ;
        
        Mintf_loc = Mintf(IND_DOFS,IND_DOFS) ;
      %  Mintf_edge = Mintf_loc(1:2:end,1:2:end)  ;
        
        
        PhiDEF_b_edge = PhiDEF_b(IND_DOFS,:) ;
        
        DATALOC.RELATIVE_SVD = 1;
        
        DATALOC.TOL = TOL_ALL ;
        [UUall,SSall,VVall] = WSVDT(PhiDEF_b_edge,Mintf_loc,DATALOC) ;
        
        PhiDEF_b_edge_filter = UUall*diag(SSall) ;
        
        
        nDEG = 10;
        %         xiALL = zeros(length(xi),nDEG+1) ;
        %         for ideg = 0:nDEG
        %             xiALL(:,ideg+1) = xi.^ideg;
        %             nx = sqrt(xiALL(:,ideg+1)'*Mintf_edge*xiALL(:,ideg+1)) ;
        %              xiALL(:,ideg+1) =  xiALL(:,ideg+1)/nx ;
        %         end
        
        
        DIRECTIONS = {'X','Y'};
        TOLedge = 1e-2;
        
        
        
        REQUIRED_DEGREE = 0; 
        
        for idirection = 1:2
            
            subplot(2,1,idirection)
            hold on
            title(['Modes of PhiDEF edge, comp. (',DIRECTIONS{idirection},') versus   normalized coordinates; edge number ',num2str(iedge)])
            
            Mintf_edge = Mintf_loc(1:2:end,1:2:end)  ;
            UUloc  = PhiDEF_b_edge_filter(idirection:2:end,:) ;
            DATALOC.TOL = TOLedge ;
            [UUloc,SSS,VVV] = WSVDT(UUloc,Mintf_edge,DATALOC ) ;
            
            % SSSloc = zeros(size(SSall)) ;
            %  SSSloc(1:length(SSS)) = SSS;
            
            %         DATALOC.RELATIVE_SVD = 1;
            %         TOL = 1e-2;
            %         [UU,SS,VV] = SVDT(PhiDEF_b_edge_x,TOL,DATALOC) ;
            hh= zeros(size(SSS)) ;
            LL = cell(size(SSS)) ;
            for imode =1:length(SSS)
                hh(imode) = plot(xi,UUloc(:,imode)) ;
                LL{imode} = ['\lambda = ',num2str(SSS(imode))] ;
            end
            legend(hh,LL)
            
            %   PROJECT = xiALL'*Mintf_edge*UUloc ;
            PLOT_FUNCTIONS = 0;
            
            %  [Q,M] =  OrthPoly_ChatGPT2(xi,nDEG,PLOT_FUNCTIONS)
            %  Q = OrthPoly_ChatGPT1(xi,nDEG,Mintf_edge,PLOT_FUNCTIONS) ;
            nDEGmax = min(nDEG,length(xi)-2) ;
            ORTHOGONALITY_L2 = 0; 
            [Q,M] =  OrthPoly_ChatGPT2(xi,nDEGmax,PLOT_FUNCTIONS,ORTHOGONALITY_L2);
            %    Q = OrthPoly_ai(xi,nDEGmax,Mintf_edge,PLOT_FUNCTIONS) ;
            PROJECT = Q'*M*UUloc ;
            % What is this ?  This is the projection of the SVD modes onto
            % an orthogonal Legrende basis for polynomials.  
            % We want, for each edge, to infer, given an error threshold,
            % what should be the polynomial which describes best such SVD
            % modes
            nPROJECTION = sqrt(sum(PROJECT.^2,1)) ; 
            PROJECT_norm = bsxfun(@times,PROJECT',1./nPROJECTION')' ; 
            cSUM_proj = sqrt(cumsum(PROJECT_norm.^2,1)) ; 
            
            TOL_FILTER=  0.01;  
            
            for imodes = 1:size(cSUM_proj,2)
               [aaa,bbb] =  find(cSUM_proj(:,imodes) >= 1-TOL_FILTER) ; 
               REQUIRED_DEGREE = max(REQUIRED_DEGREE,aaa(1)-1) ; 
            end
            
            
            
            
        end
        
        REQUIRED_degree_all_edges(iedge)= max(REQUIRED_DEGREE,1) ; 
    end
    disp('***************************************************************************************++')
    disp(['REQUIRED DEGREE polynomial EDGES = ',num2str(REQUIRED_degree_all_edges')])
    disp('***************************************************************************************++')
    disp('***************************************************************************************++')
    disp(['EMPLOYED DEGREE polynomial EDGES = ',num2str(INFO_EDGES.NumberNodesPerEdge-1)])
    disp('***************************************************************************************++')
    
end