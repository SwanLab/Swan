function   [DATAOUT,BASES,TEXTP ]= RVEStiffnessMatrix(BASES,DATA_REFMESH,DATAIN,TEXTP)

if nargin == 0
    load('tmp2.mat')
    %    DATAIN.INCLUDE_FLUCTUATION_MODES_FACES = 1;
end

ndim = size(DATA_REFMESH.COOR,2);
% Deformational displacement modes
% -----------------------------------
BasisUdef = BASES.DISPLACEMENTS.U ;  % Basis displacements
SinvVal_Udef = BASES.DISPLACEMENTS.S ; % aSSOCIATED SINGULAR VALUES
DATAIN = DefaultField(DATAIN,'nMODES_DISP',size(BasisUdef,2)) ;
BasisUdef = BasisUdef(:,1:DATAIN.nMODES_DISP)  ; %  Displacement modes (nMODES_DISP is specified in the input data file)
DATAOUT.BasisUdef = BasisUdef ;
% Self-equilibrated displacement modes
% -----------------------------------
BasisRdef = BASES.REACTIONS.U ;  % Basis displacements
SinvVal_Rdef = BASES.REACTIONS.S ; % aSSOCIATED SINGULAR VALUES

%DATAIN.KINEMATIC_CONSTRAINTS_MODES.ACTIVE = 1;   % New variable (21st-Apr-2019). Enable kinematic contraints
DATAIN = DefaultField(DATAIN,'KINEMATIC_CONSTRAINTS_MODES',[]) ;
DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'ACTIVE',0) ;


DATAIN = DefaultField(DATAIN,'nMODES_REACTIONS',size(BasisRdef,2)) ;
DATAIN = DefaultField(DATAIN,'EQUAL_NUMBER_REACTIONS_NUMBER_MODES_REACTION_DRIVEN',1) ;

if DATAIN.EQUAL_NUMBER_REACTIONS_NUMBER_MODES_REACTION_DRIVEN  == 1
    % Number of reaction modes cannot be greater than number of
    % displacement modes
    DATAIN.nMODES_REACTIONS = min(DATAIN.nMODES_REACTIONS,DATAIN.nMODES_DISP) ;
end

BasisRdef = BasisRdef(:,1:DATAIN.nMODES_REACTIONS)  ;
SinvVal_Rdef = SinvVal_Rdef(1:size(BasisRdef,2))  ;
BASES.REACTIONS.U = BasisRdef ;
BASES.REACTIONS.S =SinvVal_Rdef ;
% -----------------------------------------------------------------
DATAIN = DefaultField(DATAIN,'nMODES_STRESSES',size(BASES.STRESSES.U,2)) ;
DATAIN.nMODES_STRESSES = min(DATAIN.nMODES_STRESSES,size(BASES.STRESSES.U,2)) ;
BASES.STRESSES.U = BASES.STRESSES.U(:,1:DATAIN.nMODES_STRESSES) ;
BASES.STRESSES.S =  BASES.STRESSES.S(1:DATAIN.nMODES_STRESSES) ;

% PERIODIC STRUCTURE --> Number of total modes must be an even number
if DATAIN.KINEMATIC_CONSTRAINTS_MODES.ACTIVE == 1
    if ndim == 2
        nRB = 3;  % Number of rigid body modes
    else
        nRB  = 6 ;
    end
    nTOTAL =  nRB + DATAIN.nMODES_REACTIONS ;
    if mod(nTOTAL,2) ~=0
        DATAIN.nMODES_REACTIONS = DATAIN.nMODES_REACTIONS-1;
    end
    BasisRdef = BasisRdef(:,1:DATAIN.nMODES_REACTIONS)  ;
end

% Rigid body modes
% -----------------------------------
BasisUrb =   DATA_REFMESH.BasisUrb ;
% % Stiffness matrix
%K = DATA_REFMESH.K ;
disp('Retrieving Stiffness matrix')
tic
load(DATAIN.NAME_WS_MODES,'K');
toc
disp('...Done')
% Nodes faces f1 and f2, f3 , f4
NODES_faces = DATA_REFMESH.NODES_faces ;
NODES_faces_all = cell2mat(NODES_faces) ;
f = small2large(NODES_faces_all,ndim) ;
% Degrees of freedom faces f1 and f2, f3 and f4
for i=1:length(NODES_faces)
    fI{i} = small2large(NODES_faces{i},ndim) ;
    DATAOUT.fI{i} =fI{i};
end
DATAOUT.f = f;
%-------------------------------------------
DATAIN = DefaultField(DATAIN,'DeformationalInterfaceModes_AlignmentMethod',0) ;
DATAIN = DefaultField(DATAIN,'INTERFACE_MODES_REACTIONS_WORK',[] )  ; %.ACTIVE
DATAIN.INTERFACE_MODES_REACTIONS_WORK = DefaultField(DATAIN.INTERFACE_MODES_REACTIONS_WORK,...
    'ACTIVE',0)  ; %.ACTIVE
FACES_GROUPS = {[1,3],[2,4]} ;

DD = cell(size(fI)) ;
DATA_REFMESH = DefaultField(DATA_REFMESH,'RotationMatrixFace',DD) ;
%SROTrows = [] ;
%if ~isempty(DATA_REFMESH.RotationMatrixFace)
[SROTrows,SROTcols ]= cellfun(@size,DATA_REFMESH.RotationMatrixFace) ;
% end


DATAIN = DefaultField(DATAIN,'DeformationalInterfaceModes_From_Other_Project',[]) ; 


if isempty(DATAIN.DeformationalInterfaceModes_From_Other_Project) 
% Determination of reaction and interface modes (encapsulated 25-June-2019)
[DATAOUT,BASES,V,TEXTP,BasisUdef,BasisRdef,DATAIN]  = ...
    DetermineInterfaceModesRVEgen(DATAIN,SROTrows,DATA_REFMESH,DATAOUT,BasisUdef,SinvVal_Udef,...
    f,fI,SinvVal_Rdef,FACES_GROUPS,BASES,BasisRdef,TEXTP) ;
else
    % Interface  modes loaded from other project         
   load(DATAIN.DeformationalInterfaceModes_From_Other_Project,'DATAROM')  ; 
   TEXTP = {}  ;
   [~,NNN,~] = fileparts(DATAIN.DeformationalInterfaceModes_From_Other_Project) ; 
   TEXTP{end+1} = ['Interface modes from project:  ',NNN] ; 
      TEXTP{end+1} = ['Number of disp. modes:  ',num2str(size(BasisUdef,2))] ; 

   DATAOUT.BasisIntRB = DATAROM.BasisIntRB ; 
   DATAOUT.BasisInt = DATAROM.BasisInt ; 
   V  = DATAROM.BasisInt  ; 
   DATAOUT.BasisUdef =BasisUdef ;
   DATAOUT.BasisRdef  =BasisRdef ; 
   DATAOUT.SingVal_Udef  =SinvVal_Udef ; 
   
end


%%%%
%
%     [DATAOUTold] = LOAD_DATAOUT(DATAIN.DeformationalInterfaceModes_From_Other_Project)
% end



% PLOT CANDIDATES FOR BEING  INTERFACE MODES
% ---------------------------------------------
for ifgroup=1:length(FACES_GROUPS)
    refFACE = FACES_GROUPS{ifgroup}(1) ;
    
    COOR =DATA_REFMESH.COOR(NODES_faces{refFACE},:) ;
    CNref =  RenumberConnectivities( DATA_REFMESH.CONNECTb{refFACE},1:length(NODES_faces{refFACE})) ;
    TypeElementB = DATA_REFMESH.TypeElementB ;
    posgp = [] ;
    LEGENDG= ['INTF_MODES_',num2str(FACES_GROUPS{ifgroup}(1)),'_',num2str(FACES_GROUPS{ifgroup}(2))] ;
    NAME_MODES = [DATAIN.NAME_WS_MODES(1:end-4),LEGENDG ];
    
    DATAMODES = [] ;
    if isstruct(V)
        Vplot = full(V.BasisINTFall_cell{refFACE});
    else
        Vplot = V{refFACE};
    end
   TEXTP =  GidPostProcessModes_dom(COOR,CNref,TypeElementB,Vplot,posgp,NAME_MODES,DATAMODES,LEGENDG,TEXTP);
    
end





% -------------------------------------------------------------
KdomRED = BasisUdef'*(K*BasisUdef) ;  % Reduced stiffness matrix
BasisRdef_f = BasisRdef(f,:) ;   % REduced-stiffness at the boundary
BasisUdef_f = BasisUdef(f,:) ;
%  \Hqr{e} =  \BasisUdefT{e}{f}  \BasisRdef{e}{f}
Hqr = BasisUdef_f'*BasisRdef_f ;
DATAOUT.HqrRB =  BasisUrb(f,:)'*BasisRdef_f ; % 13-Jun-2019

if rank(Hqr) < size(BasisRdef_f,2)
    error('Ill-conditioned matrix. Probably the training set has rank less than 18 (if it is a 3D problem)')
end
DATAOUT.Hqr = Hqr ;
DATAOUT.KdomRED = KdomRED ;
%  \Kbeam{e}{} \defeq (\Hqr{e^T}\KdomRED{e^{-1}}{} \Hqr{e})^{-1},
Kbeam_inv = Hqr'*(KdomRED\Hqr) ;
Kbeam = inv(Kbeam_inv) ;
DATAOUT.Kbeam = Kbeam ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \BasisINTdom{1 }{e^T} (\BasisRdef{e}{\face{1}}\Kbeam{e}{}\BasisRdef{e^T}{\face{1}}) \BasisINTdom{1}{e}
BasisRrb  = DATA_REFMESH.M*BasisUrb ;  ; %

if  ~isstruct(V)
    % OLD METHOD (but good one), NOT CONSTRAINED KINEMATICALLY
    VB = cell(1,length(V)) ;
    VB_rb = VB;
    for iface = 1:length(V)
        %  ROTATION !!!!! ()
        if ~isempty(DATA_REFMESH.RotationMatrixFace{iface})
            Rotation_Intf_domain = DATA_REFMESH.RotationMatrixFace{iface} ;
            VB{iface}   = RotateMatricesProduct_interface(Rotation_Intf_domain',BasisRdef(fI{iface},:),V{iface}) ;
            VB_rb{iface}   = RotateMatricesProduct_interface(Rotation_Intf_domain',BasisRrb(fI{iface},:),V{iface})' ;
        else
            VB{iface}  = V{iface}'*BasisRdef(fI{iface},:)  ;
            VB_rb{iface} = (V{iface}'*BasisRrb(fI{iface},:))' ;
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reduced stiffness matrix
    % -------------------------
    Kskel =cell(length(VB),length(VB)) ;
    % Kskel_11 = VB1*Kbeam*VB1' ;
    % Kskel_12 = VB1*Kbeam*VB2' ; ;
    % Kskel_21 = VB2*Kbeam*VB1' ; ;
    % Kskel_22 =VB2*Kbeam*VB2' ; ;
    for iii = 1:length(VB)
        for jjj=1:length(VB)
            Kskel{iii,jjj} = VB{iii}*Kbeam*VB{jjj}';
        end
    end
    %%%% STIFNESS MATRIX % ------------------------------
    Kskel = cell2mat(Kskel);
    DATAOUT.Tcomp = VB ;
    DATAOUT.TcompRB = VB_rb ; % 13-Jun-2019,
    
    Tcomp = cell2mat(VB')' ;
    
    
    
    TEXTP{end+1} = '________________ SUMMARY __________________' ;
    TEXTP{end+1} = ['Number of interface modes = ', num2str(size(V{1},2)), ' ',num2str(size(V{2},2)), ' ',num2str(size(V{3},2)), ' ',num2str(size(V{4},2)), ' '] ;
    
else
    
    % New method
    warning('This method is not reliable !!!! ')
    Tcomp = BasisRdef(f,:)'*V.BasisINTF  ;
    Kskel =  Tcomp'*Kbeam*Tcomp ;
    
    
    
    
end
DATAOUT.Kskel  = Kskel ;



% Tcomp matrix
[UT,ST,VT] = SVDT(Tcomp) ;

TEXTP{end+1} = ['RANK Tcomp = ',num2str(rank(Tcomp))] ;
TEXTP{end+1} = ['s(end)/s(end-1) = ',num2str(ST(end)/ST(end-1))] ;






% ------------------------------------------------------
nBASES_RVE.REACTIONS = size(BasisRdef,2) ;
nBASES_RVE.DISPLACEMENTS = size(BasisUdef,2) ;
% % Stiffness constant matrix
% BeamPropertiesStiffness(BasisRdef(:,1:nBASES_RVE.REACTIONS),f1,f2,Vrb,L,AREA,DATA_REFMESH,DATAIN,...
%     BasisUdef(:,1:nBASES_RVE.DISPLACEMENTS),K) ;
% -----------------------------------------------------
% Operators relating external forces in the ROM and the full-order model
if ~isstruct(V)
    [DATAOUT]= ExtForcesROMrve_operators(BasisUrb,fI,f,BasisRdef,Kbeam,Hqr,KdomRED,...
        DATA_REFMESH,BasisUdef,V,ndim,DATAOUT,DATAIN) ;
else
    % New method, Apr-2019, kinematically constrained  (copy of .... ExtForcesROMrve_operPLATES)
    
    [DATAOUT]= ExtForcesROMrve_operKINEM(BasisUrb,fI,f,BasisRdef,Kbeam,Hqr,KdomRED,...
        DATA_REFMESH,BasisUdef,V,ndim,DATAOUT,DATAIN) ;
    %     [DATAOUT]= ExtForcesROMrve_operKINEM(BasisUrb,fI,f,BasisRdef,Kbeam,Hqr,KdomRED,...
    %     DATA_REFMESH,BasisUdef,V_fc,ndim,DATAOUT,DATAIN)
end
% -----------------------------------
% % FLUCTUATION MODES
% % -----------------
% DATAIN = DefaultField(DATAIN,'IS_A_JOINT',0) ;
% if DATAIN.IS_A_JOINT == 0
%     [DATAOUT.BasisINTfluct,DATAOUT.kFLUCT ]= ...
%         ComputeFluctuationModes(BasisUdef,f1,nBASES_RVE,DATAIN,Vrb,M,...
%         NODES_faces12,DATA_REFMESH,DATAsnap);
% end
% % -----------------------------------
% end


% FILE_BATCH  = [cd,'/INFO.txt'] ;
% fid =fopen(FILE_BATCH,'w');
% for i = 1:length(TEXTP)
%     fprintf(fid,[TEXTP{i},'\n']);
% end
% fod =fclose(fid);
%
% open(FILE_BATCH);




