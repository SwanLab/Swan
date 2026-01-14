function   [BasisU,SS,VV] =  DominantModesVelocityAndDisplacements(DATAoffline,DATAparamSTUDY,SNAPdisp,SNAPvel)

% Dynamic problems. Determination of basis matrix for displacements taking
% into account velocities as well
if nargin == 0
    load('tmp1.mat')
end


if DATAparamSTUDY.MethodDominantModes_DISP_AND_VELOCITIES == 1
    
    
    DATAsvd.EPSILON_GLO  =DATAoffline.errorDISP ;
    
    TOL_BLOCK = zeros(size([SNAPdisp;SNAPvel])) ;
    
    for istore = 1:length(SNAPdisp)
        normDISP = norm(SNAPdisp{istore},'fro') ;
        SNAPdisp{istore} = SNAPdisp{istore}/normDISP ;
        normVEL= norm(SNAPvel{istore},'fro') ;
        SNAPvel{istore} = SNAPvel{istore}/normVEL ;
    end
    
    
    [BasisUmixed,Svel,Vvel,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp( [SNAPdisp;SNAPvel],TOL_BLOCK,DATAsvd) ;
    disp('***********************************************************')
    disp(['Number of displacement/velocity modes =',num2str(size(BasisUmixed,2)), ' (for ERRORdisp/vel = ',num2str(DATAoffline.errorDISP),')'])
    disp('***********************************************************')
    disp(['Singular Values =',num2str(Svel')])
    
    ROWS_disp = 1:size(SNAPdisp{1},1) ;
    ROWS_vel=  size(SNAPdisp{1},1)+1:2* size(SNAPdisp{1},1) ;
    
    BasisUdisp = BasisUmixed(ROWS_disp,:) ;
    
    %   [Udi,SSSd,VVVd] = SVDT(BasisUdisp) ;
    
    BasisUvel = BasisUmixed(ROWS_vel,:) ;
    
    %  [Uvi,SSSv,VVVv] = SVDT(BasisUvel) ;
    
    
    % INTersection
    
    %  [Uint,Sint,Vint] =    SVDT(Uvi'*Udi)  ;
    
    
    DATAloc.RELATIVE_SVD =1 ;
    TOL =  1e-7;
    
    [BasisU,SS,VV] = SVDT([BasisUdisp,BasisUvel],TOL,DATAloc) ;
    
%     SNAPdisp = cell2mat(SNAPdisp) ; 
%     
%     ERROR = SNAPdisp-BasisU*(BasisU'*SNAPdisp) ; 
%     nERROR = norm(ERROR,'fro')/norm(SNAPdisp,'fro')*100
    
    
else
    
    
    TOL_BLOCK = [DATAparamSTUDY.errorDISP*ones(size([SNAPdisp])),DATAparamSTUDY.errorVELOC*ones(size([SNAPdisp]))  ] ;
    
    for istore = 1:length(SNAPdisp)
        normDISP = norm(SNAPdisp{istore},'fro') ;
        SNAPdisp{istore} = SNAPdisp{istore}/normDISP ;
        normVEL= norm(SNAPvel{istore},'fro') ;
        SNAPvel{istore} = SNAPvel{istore}/normVEL ;
    end
    
    DATAsvd = [] ;
    [BasisU,SS,VV,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp( [SNAPdisp,SNAPvel],TOL_BLOCK,DATAsvd) ;
    disp('***********************************************************')
    disp(['Number of displacement  modes =',num2str(size(BasisU,2)), ' (for ERRORdisp/vel = ',num2str(DATAoffline.errorDISP),';',num2str(DATAparamSTUDY.errorVELOC),')'])
    disp('***********************************************************')
    disp(['Singular Values =',num2str(SS')])
    
    %     ROWS_disp = 1:size(SNAPdisp{1},1) ;
    %     ROWS_vel=  size(SNAPdisp{1},1)+1:2* size(SNAPdisp{1},1) ;
    %
    %     BasisUdisp = BasisUmixed(ROWS_disp,:) ;
    %
    %     BasisUvel = BasisUmixed(ROWS_vel,:) ;
    %
    %     DATAloc.RELATIVE_SVD =1 ;
    %     TOL = 1e-6;
    %
    %     [BasisU,SS,VV] = SVDT([BasisUdisp,BasisUvel],TOL,DATAloc) ;
    
    
    
    
    
end

