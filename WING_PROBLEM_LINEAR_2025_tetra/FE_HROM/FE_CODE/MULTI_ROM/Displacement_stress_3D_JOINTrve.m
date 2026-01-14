function [DISP3D,DISP3D_lateral,stressGLO,STRESS_DATA,DATAIN,DATA_REFMESH] ...
    = Displacement_stress_3D_JOINTrve(DATAROM,MESH2D,qDEF,qRB,DATAIN,DATA_REFMESH,DATAadd,DATARUN)


if nargin ==0
    load('tmp0.mat')
end


DATAIN = DefaultField(DATAIN,'DOMAINS_POSTPROCESS_SELECT',[]) ;
DATAIN.DOMAINS_POSTPROCESS_SELECT = DefaultField(DATAIN.DOMAINS_POSTPROCESS_SELECT,'NUMBER',[]) ;
DATAIN.DOMAINS_POSTPROCESS_SELECT = DefaultField(DATAIN.DOMAINS_POSTPROCESS_SELECT,'VARIABLE','VONMISES') ;
DATAIN = DefaultField(DATAIN,'POST_PROCESS_LATERAL_SURFACES',[]) ;
DATAIN.POST_PROCESS_LATERAL_SURFACES = DefaultField(DATAIN.POST_PROCESS_LATERAL_SURFACES,'ACTIVE',0) ;
DATAIN.POST_PROCESS_LATERAL_SURFACES = DefaultField(DATAIN.POST_PROCESS_LATERAL_SURFACES,'COARSE_MESH',1) ;
%DATAIN.POST_PROCESS_LATERAL_SURFACES = DefaultField(DATAIN.POST_PROCESS_LATERAL_SURFACES,'NAME_MESH_COARSE',[]) ;

DATAIN = DefaultField(DATAIN,'PRINT_DISTRIBUTED_FORCES',0) ; % Print in GID external traction forces (over external surfaces)
% If this option is enabled, then no inner elements are plotted (only lateral surfaces)
% Likewise, the plot is made using the finer mesh
nDOM = size(qDEF,2) ;
if DATAIN.PRINT_DISTRIBUTED_FORCES == 1
    DATAIN.DOMAINS_POSTPROCESS = [] ;
    DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER = [] ;
    NO_STRESSES = 1;
else
    DATAIN = DefaultField(DATAIN,'DOMAINS_POSTPROCESS',[]) ;
    if isempty(DATAIN.DOMAINS_POSTPROCESS)  && isempty(DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER)
        DATAIN.DOMAINS_POSTPROCESS = 1:nDOM ;
    end
    NO_STRESSES = 0;
end









DISP3D = [] ;
DISP3D_lateral = [] ;
%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%DISP3D_all = DISP3D ;




%--------------
%% STRESSES
% ---------
[stressGLO,STRESS_DATA,SELECTED_DOMAINS] = ...
    StressDom_SlicesRVE(DATAROM,MESH2D,DATAIN,qDEF,DATA_REFMESH,NO_STRESSES,DATAadd,DATARUN) ;

%% DISPLACEMENTS
% ------------------
DISP3D = ComputeDisplacementsRVE(DATAIN,DATAROM,MESH2D,qRB,DATA_REFMESH,qDEF,SELECTED_DOMAINS) ;


% DISPLACEMENT LATERAL SURFACES
% --------------------------------
[DISP3D,DATAIN ]= DisplacementLateralSurfacesRVE(DATAIN,qDEF,DATA_REFMESH,DISP3D,SELECTED_DOMAINS,qRB,DATAROM,MESH2D) ;
%b
%----------------------------------------------------------------------------------

DATAIN.SELECTED_DOMAINS = SELECTED_DOMAINS ;

end

% --------------------------------



function  DISP3D = ComputeDisplacementsRVE(DATAIN,DATAROM,MESH2D,qRB,DATA_REFMESH,qDEF,SELECTED_DOMAINS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATAIN = DefaultField(DATAIN,'POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL',0) ;
nentities =length(DATAROM) ;
DISP3D = cell(1,nentities) ;
for ientity = 1:length(DATAROM)
    ELEMS = find(MESH2D.MaterialType == ientity) ;
    ELEMS = ELEMS(SELECTED_DOMAINS{ientity}) ;
    
    % ELEMS =
    BasisUdef = DATAROM{ientity}.BasisUdef ;
    BasisUrb = DATA_REFMESH{ientity}.BasisUrb;
    nDOM = length(ELEMS) ;
    DISP3D{ientity} =zeros(size(BasisUdef,1),nDOM) ;
    if DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 0
        DISP3D{ientity} = BasisUdef*qDEF(:,ELEMS) + BasisUrb*qRB(:,ELEMS) ;
    elseif DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 1
        DISP3D{ientity}=   BasisUrb*qRB(:,ELEMS);
    elseif DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 2
        DISP3D{ientity} =   BasisUdef*qDEF(:,ELEMS)  ;
    else
        error('OPtion not implemented')
    end
    
    %%% ROTATION OF DISPLACEMENTS   (adaption from BEAM's implementation)
    % ---------------------------
    DATA_REFMESH{ientity} = DefaultField( DATA_REFMESH{ientity},'RotationMatrixFace',[]) ;
    if ~isempty(DATA_REFMESH{ientity}.RotationMatrixFace)
        [AAAA BBBB] = cellfun(@size,DATA_REFMESH{ientity}.RotationMatrixFace) ;
        if any(AAAA)
            %%% Rotation of displacements (from domain coordinates to global coordinates )
            % -----------------------------------------------------------------------------
            % We begin by implementing the non-vectorized version
          %  ndim = size(MESH2D.ROTATIONS,1) ;
            for ielemLOC = 1:length(ELEMS)
                ielem = ELEMS(ielemLOC) ;
             %   finIND =  ndim*ielem ;
              %  iniIND = ndim*ielem - ndim+1 ;
                R = MESH2D.rotDOM{ielem} ;
                ndim = size(R,1) ; 
                dispLOC =  reshape(DISP3D{ientity}(:,ielemLOC),ndim,[]) ;
                dispLOC = R*dispLOC ;
                DISP3D{ientity}(:,ielemLOC)  =dispLOC(:) ;
            end
        end
    end
    
    
    
end


  




end


 

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function  [stressGLO_glo,stressVONMISES_glo,MAXstressVONMISES] = ...
%     FullOrderStresses(ientity,DATAIN,DATA_REFMESH,stressGLO_glo,stressVONMISES_glo,qDEF,DATAROM,ELEMS,nstrain,...
%     SELECTED_DOMAINS)
% 
% 
% 
% % ---------------------
% %% FULL-ORDER MODEL
% %% --------------------
% disp('Retrieving OFFLINE data')
% tic
% load(DATAIN.NAME_WS_MODES{ientity},'CgloDOM','Wdom','Bdom')
% toc
% disp('DONE')
% Bdom = Bdom*DATAROM{ientity}.BasisUdef ;
% CgloDOM = CgloDOM*Bdom ;
% 
% if DATAIN.PRINT_AVERAGE_STRESSES_ON_ELEMENTS == 0
%     error(['Option not implemented'])
% end
% ngaus = size(DATA_REFMESH{ientity}.posgp,2) ;
% Wst = repmat(Wdom',nstrain,1) ;
% Wst = Wst(:) ;
% nelem =  length(Wdom)/ngaus ;
% 
% if ~isempty(SELECTED_DOMAINS)
%     ELEMS = ELEMS(SELECTED_DOMAINS{ientity});
% end
% 
% 
% stressGLO_glo{ientity} = zeros(nelem*nstrain,size(length(ELEMS),2)) ;
% stressVONMISES_glo{ientity}= zeros(nelem,size(length(ELEMS),2)) ;
% MAXstressVONMISES = zeros(1,size(length(ELEMS),2))  ;
% 
% 
% 
% if DATAIN.DO_NOT_COMPUTE_STRESSES == 0
%     tic
%     disp('Computing average stress on each FE element')
%     
%     
%     for ielem = 1:length(ELEMS)
%         idom = ELEMS(ielem) ;
%         disp(['Domain =',num2str(ielem)])
%         
%         stressDOM_e = CgloDOM*qDEF(:,idom) ;  %
%         for istrain =1:nstrain
%             stressDOM_e(istrain:nstrain:end) = stressDOM_e(istrain:nstrain:end)./Wdom ;   % Cglo already includes WEIGHTS
%         end
%         [stressDOM_e] = AverageStressOnElements(stressDOM_e,Wst,nelem,nstrain,ngaus) ;
%         
%         % Von Mises
%         stressDOM_ELEM = reshape(stressDOM_e,nstrain,[]) ;
%         [ stressVONMISES_e ] =  VonMises_Stress(stressDOM_ELEM) ;
%         MAXstressVONMISES(ielem) = max(stressVONMISES_e) ;
%         stressGLO_glo{ientity}(:,ielem) = stressDOM_e ;
%         stressVONMISES_glo{ientity}(:,ielem) = stressVONMISES_e' ;
%     end
%     disp('...Done')
%     toc
% else
%     stressVONMISES=[] ;
%     stressGLO = [] ;
%     MAXstressVONMISES = [] ;
%     
% end
% 
% end



%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function  [MAXstressVONMISES] = ...
%     ReducedOrderStresses(ientity,DATAIN,DATA_REFMESH,stressGLO_glo,stressVONMISES_glo,...
%     qDEF,DATAROM,ELEMS,nstrain)
% 
% 
% % ---------------------
% %% FULL-ORDER MODEL
% %% --------------------
% % disp('Retrieving OFFLINE data')
% % tic
% % load(DATAIN.NAME_WS_MODES{ientity},'CgloDOM','Wdom','Bdom')
% % toc
% % disp('DONE')
% %%%%%%%%%%%%%%%%%%%%%%%%
% %if ~isfield(DATAROM{ientity},'HROMVAR')
% Celas_Bdom = DATAROM{ientity}.HROMVAR.Celas_Bdom;   % Product Celas times B matrix (reduced points )
% %end
% 
% if DATAIN.PRINT_AVERAGE_STRESSES_ON_ELEMENTS == 0
%     error(['Option not implemented'])
% end
% % ngaus = size(DATA_REFMESH{ientity}.posgp,2) ;
% % Wst = repmat(Wdom',nstrain,1) ;
% % Wst = Wst(:) ;
% % nelem =  length(Wdom)/ngaus ;
% 
% 
% 
% % stressGLO_glo{ientity} = [] ;
% % stressVONMISES_glo{ientity}= [] ;
% MAXstressVONMISES = zeros(1,size(length(ELEMS),2))  ;
% 
% 
% 
% if DATAIN.DO_NOT_COMPUTE_STRESSES == 0
%     tic
%     disp('Computing average stress on each FE element')
%     
%     DATAIN = DefaultField(DATAIN,'COMPUTE_STRESSES_VECTOR_FORM',1) ;
%     
%     if DATAIN.COMPUTE_STRESSES_VECTOR_FORM == 0
%         for ielem = 1:length(ELEMS)
%             idom = ELEMS(ielem) ;
%             disp(['Domain =',num2str(ielem)])
%             % Stresses at the selected Gauss points
%             stressDOM_e = Celas_Bdom*qDEF(:,idom) ;
%             stressDOM_ELEM = reshape(stressDOM_e,nstrain,[]) ;
%             [ stressVONMISES_e ] =  VonMises_Stress(stressDOM_ELEM) ;
%             [MAXstressVONMISES(ielem) INNDDDDD]= max(stressVONMISES_e) ;
%             % stressGLO_glo{ientity}(:,ielem) = stressDOM_e ;
%             % stressVONMISES_glo{ientity}(:,ielem) = stressVONMISES_e' ;
%         end
%     else
%         stressDOM_e = Celas_Bdom*qDEF ;
%         stressDOM_ELEM = reshape(stressDOM_e(:),nstrain,[]) ;
%         [ stressVONMISES_e ] =  VonMises_Stress(stressDOM_ELEM) ;
%         stressVONMISES_e = reshape(stressVONMISES_e,[],size(qDEF,2));
%         MAXstressVONMISES = max(stressVONMISES_e) ;
%     end
%     
%     
%     disp('...Done')
%     toc
% else
%     %     stressVONMISES=[] ;
%     %     stressGLO = [] ;
%     MAXstressVONMISES = [] ;
%     
% end
% 
% end





