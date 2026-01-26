function   [DATAOUT,MSG] = BeamStiffnessMatrix_MODESnu_g_nr(BASES,DATA_REFMESH,DATAIN,nBASES_BEAM,DATAsnap,MSG)
% Modification of BeamStiffnessMatrix.m
% GOAL: Explore what happens if nU >nR
% MAy-13th-2020 (61th QUARANTINE COVID-19)

if nargin == 0
    load('tmp1.mat')
end


ndim = size(DATA_REFMESH.COOR,2);
% Deformational displacement modes
% -----------------------------------
BasisUdef = BASES.DISPLACEMENTS.U ;  % Basis displacements
SinvVal_Udef = BASES.DISPLACEMENTS.S ; % aSSOCIATED SINGULAR VALUES
DATAIN = DefaultField(DATAIN,'nMODES_DISP',size(BasisUdef,2)) ;
BasisUdef = BasisUdef(:,1:DATAIN.nMODES_DISP)  ;   % Truncation (prescribed by the user)
DATAOUT.BasisUdef = BasisUdef ;
% Self-equilibrated displacement modes
% -----------------------------------
BasisRdef = BASES.REACTIONS.U ;  % Basis reactions
DATAIN = DefaultField(DATAIN,'nMODES_REACTIONS',size(BasisRdef,2)) ;
BasisRdef = BasisRdef(:,1:DATAIN.nMODES_REACTIONS)  ;
DATAOUT.BasisRdef = BasisRdef ;
if ~isempty(BASES.REACTIONS.S) && length(BASES.REACTIONS.S) >= DATAIN.nMODES_REACTIONS
    SinvVal_Rdef = BASES.REACTIONS.S(1:DATAIN.nMODES_REACTIONS) ; % aSSOCIATED SINGULAR VALUES
else
    SinvVal_Rdef = ones(length(DATAIN.nMODES_REACTIONS),1) ;
end

% Rigid body modes (slice)
% -----------------------------------
BasisUrb =   DATA_REFMESH.BasisUrb ;
nrigid = size(BasisUrb,2) ;

% % Stiffness matrix
%K = DATA_REFMESH.K ;
disp('Retrieving Stiffness matrix')
tic
load(DATAIN.NAME_WS_MODES,'K');
toc
disp('...Done')
% Nodes faces f1 and f2
NODES_faces12 = DATA_REFMESH.NODES_faces12 ;
% Degrees of freedom faces f1 and f2
f1 = small2large(NODES_faces12{1},ndim) ;
f2 = small2large(NODES_faces12{2},ndim) ;
f = [f1;f2] ;
DATAOUT.f1 = f1 ;
DATAOUT.f2 = f2 ;
%-------------------------------------------
% Rigid body modes for interface
Vrb = DATA_REFMESH.RigidBodyModesInterface ;
DATAOUT.BasisIntRB = Vrb ;

nfaces = 2;
M1d = DATA_REFMESH.GeometricMassMatrixInterface ;
M = cell(2,1) ;
for iface= 1:nfaces
    M{iface} = sparse(size(M1d{iface},1)*ndim,size(M1d{iface},2)*ndim) ;
    for idim = 1:ndim
        M{iface}(idim:ndim:end,idim:ndim:end) = M1d{iface} ;
    end
end



iface= 1;
BasisINTbeam = Vrb{iface} ;
DATAIN = DefaultField(DATAIN,'TOL_SINGULAR_VALUES_Hqr',0.1) ;
nBOUNDARY_INTFMODES = [] ;
DATAIN = DefaultField(DATAIN,'DoNotUseRigidBodyForInterface',0) ;
DATAIN = DefaultField(DATAIN,'INTERFACE_MODES_REACTIONS_WORK',[]) ;
DATAIN.INTERFACE_MODES_REACTIONS_WORK = DefaultField(DATAIN.INTERFACE_MODES_REACTIONS_WORK, ...
    'ACTIVE',0) ;

[BasisRdef,V,MSG,BasisUdef,SinvVal_Udef] = ReactionAndInterfaceLocalModes_nu_g_nr(BasisUdef,BasisRdef,f1,f2,DATAIN.TOL_SINGULAR_VALUES_Hqr,...
    nBASES_BEAM,DATA_REFMESH,BasisINTbeam,M{iface},DATAIN,SinvVal_Udef,SinvVal_Rdef,BasisUrb,...
    DATA_REFMESH.M,MSG)  ;


MSG{end+1} = ['**********************************************+'] ;
MSG{end+1} = ['**********************************************+'] ;
MSG{end+1} = ['Number of  displacement modes =',num2str(size(BasisUdef,2))];
MSG{end+1} = ['Number of  reaction modes =',num2str(size(BasisRdef,2))];
MSG{end+1} =   ['Number of   interface modes =',num2str(size(V,2))];


MSG{end+1} = ['**********************************************+'] ;
MSG{end+1} = ['**********************************************+'] ;
%         for iiii = 1:length(MSGR)
%             MSG{end+1} = MSGR{iiii} ;
%         end

MSG{end+1} = ['**********************************************+'] ;
MSG{end+1} = ['**********************************************+'] ;
%         for iiii = 1:length(MSG2)
%             MSG{end+1} = MSG2{iiii} ;
%         end

%  disp(['Number of additional interface modes (boundary) =',num2str(nBOUNDARY_INTFMODES)])

DATAOUT.nBOUNDARY_INTFMODES = nBOUNDARY_INTFMODES ;
DATAOUT.BasisInt = V ;
DATAOUT.BasisRdef = BasisRdef ;
DATAOUT.BasisUdef = BasisUdef ;
DATAOUT.SingVal_Udef = SinvVal_Udef ;







% -------------------------------------------------------------
KdomRED = BasisUdef'*(K*BasisUdef) ;  % Reduced stiffness matrix, local coordinates
BasisRdef_f = BasisRdef(f,:) ;   % REduced-stiffness at the boundary
BasisUdef_f = BasisUdef(f,:) ;

%  \Hqr{e} =  \BasisUdefT{e}{f}  \BasisRdef{e}{f}
Hqr = BasisUdef_f'*BasisRdef_f ;
DATAOUT.HqrRB =  BasisUrb(f,:)'*BasisRdef_f ;

if rank(Hqr) < size(BasisRdef_f,2)
    error('Ill-conditioned matrix')
end
DATAOUT.Hqr = Hqr ;
DATAOUT.KdomRED = KdomRED ;


% DIAGNOSING_ILL_POSED = 1;
% if DIAGNOSING_ILL_POSED == 1
%     [RR,ANGLES ]= IntersectionSubspaces(BasisUdef_f,BasisRdef_f,0,0) ;
% end

% OLD_METHOD =0 ;


%if OLD_METHOD == 1

%  \Kbeam{e}{} \defeq (\Hqr{e^T}\KdomRED{e^{-1}}{} \Hqr{e})^{-1},
Kbeam_inv = Hqr'*(KdomRED\Hqr) ;
Kbeam = inv(Kbeam_inv) ;
%
% else
%     % This is how is actually implemented in the paper , eq 79
%     Hinv = inv(Hqr) ;
%     Kbeam = Hinv*KdomRED*Hinv' ;
%
% end
DATAOUT.Kbeam = Kbeam ;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \BasisINTdom{1 }{e^T} (\BasisRdef{e}{\face{1}}\Kbeam{e}{}\BasisRdef{e^T}{\face{1}}) \BasisINTdom{1}{e}

% This is the calculation of the transpose of  "T" matrix
DATA_REFMESH = DefaultField(DATA_REFMESH,'RotationMatrixFace',cell(1,2)) ;

BasisRrb  = DATA_REFMESH.M*BasisUrb ;  ; %

if iscell(V)
    iface= 1;
    VB1 = V{iface}'*BasisRdef(f1,:) ;
    VB1_rb = V{iface}'*BasisRrb(f1,:) ;
    
    iface= 2;
    Rotation_Intf_domain = DATA_REFMESH.RotationMatrixFace{2}' ;
    VB2 = RotateMatricesProduct_interface(Rotation_Intf_domain,BasisRdef(f2,:),V{iface}) ;
    VB2_rb = RotateMatricesProduct_interface(Rotation_Intf_domain,BasisRrb(f2,:),V{iface}) ;
else
    % iface= 1;
    VB1 = V'*BasisRdef(f1,:) ;
    VB1_rb = V'*BasisRrb(f1,:) ;
    % iface= 2;
    %  VB2 = V'*BasisRdef(f2,:) ;
    Rotation_Intf_domain = DATA_REFMESH.RotationMatrixFace{2}' ;
    VB2 = RotateMatricesProduct_interface(Rotation_Intf_domain,BasisRdef(f2,:),V) ;
    VB2_rb = RotateMatricesProduct_interface(Rotation_Intf_domain,BasisRrb(f2,:),V) ;
    
end

Tcomp{1} = VB1 ;
Tcomp{2} = VB2 ;
DATAOUT.Tcomp= Tcomp ;

TcompRB{1} = VB1_rb' ;
TcompRB{2} = VB2_rb' ;
DATAOUT.TcompRB = TcompRB ;

Kskel_11 = VB1*Kbeam*VB1' ;
Kskel_12 = VB1*Kbeam*VB2' ; ;
Kskel_21 = VB2*Kbeam*VB1' ; ;
Kskel_22 =VB2*Kbeam*VB2' ; ;

s1 = svd(Kskel_11) ;
MSG{end+1} = ['************************'] ;
MSG{end+1} = ['ratio  SingVal_K11 =',num2str(s1(end)/s1(1))] ;
s1 = svd(Kskel_22) ;
MSG{end+1} = ['ratio  SingVal_K22 =',num2str(s1(end)/s1(1))] ;
KK = Kskel_22+Kskel_11 ;
s1 = svd(KK) ;
MSG{end+1} = ['ratio  SingVal_(K22 +K11)=',num2str(s1(end)/s1(1))] ;
MSG{end+1} = ['************************'] ;



% Kskel_11 = V'*(BasisRdef(f1,:)*Kbeam*BasisRdef(f1,:)')*V ;
% Kskel_12 = V'*(BasisRdef(f1,:)*Kbeam*BasisRdef(f2,:)')*V ;
% Kskel_21 = V'*(BasisRdef(f2,:)*Kbeam*BasisRdef(f1,:)')*V ;
% Kskel_22 = V'*(BasisRdef(f2,:)*Kbeam*BasisRdef(f2,:)')*V ;

%% CHECKING RANK
% S = svd(Kskel_22) ;
%%%% STIFNESS MATRIX % ------------------------------
Kskel = [Kskel_11, Kskel_12; Kskel_21, Kskel_22] ;
DATAOUT.Kskel  = Kskel ;
% ------------------------------------------------------
L = DATA_REFMESH.LENGTH ;
AREA = DATA_REFMESH.AREA;
if isempty(nBASES_BEAM)
    nBASES_BEAM.REACTIONS = size(BasisRdef,2) ;
    nBASES_BEAM.DISPLACEMENTS = size(BasisUdef,2) ;
end
% Stiffness constant matrix
if isempty(DATA_REFMESH.RotationMatrixFace{2})
    MSG =   BeamPropertiesStiffness(BasisRdef(:,1:nBASES_BEAM.REACTIONS),f1,f2,Vrb,L,AREA,DATA_REFMESH,DATAIN,...
        BasisUdef(:,1:nBASES_BEAM.DISPLACEMENTS),K,Tcomp,MSG) ;
end
% -----------------------------------------------------
% Operators relating external forces in the ROM and the full-order model
[DATAOUT]= ExtForcesROM_operators(BasisUrb,f1,f2,f,BasisRdef,Kbeam,Hqr,KdomRED,DATA_REFMESH,BasisUdef,...
    V,ndim,DATAOUT,DATAIN,Tcomp) ;
% -----------------------------------
% FLUCTUATION MODES
% -----------------
DATAIN = DefaultField(DATAIN,'IS_A_JOINT',0) ;
DATAIN = DefaultField(DATAIN,'CALCULATE_FLUCTUTATION_MODES',0) ;
if DATAIN.IS_A_JOINT == 0 && DATAIN.CALCULATE_FLUCTUTATION_MODES == 1
    error('Option no longer maintained....')
    [DATAOUT.BasisINTfluct,DATAOUT.kFLUCT ]= ...
        ComputeFluctuationModes(BasisUdef,f1,nBASES_BEAM,DATAIN,Vrb,M,...
        NODES_faces12,DATA_REFMESH,DATAsnap);
else
    MSG{end+1} = '------------------------------------' ;
    MSG{end+1} = 'INTERFACE MODES' ;
    MSG{end+1} = '------------------------------------' ;
    
    MSG =  PrintInterfaceModes(V, DATAIN, DATA_REFMESH,NODES_faces12,MSG);
    MSG{end+1} = '------------------------------------' ;
    
end
% -----------------------------------









end
