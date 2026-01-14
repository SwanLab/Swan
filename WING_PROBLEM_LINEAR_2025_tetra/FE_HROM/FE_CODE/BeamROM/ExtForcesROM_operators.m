function [DATAOUT]= ExtForcesROM_operators(BasisUrb,f1,f2,f,BasisRdef,Kbeam,Hqr,KdomRED,...
    DATA_REFMESH,BasisUdef,V,ndim,DATAOUT,DATAIN,Tcomp)


% .GAMMAforceVOL,DATAOUT.UPSILONforceTRAC,...
%   DATAOUT.Lgravity,DATAOUT.Q_tractions


if nargin == 0
    load('tmp.mat')
end

% ---------------------------------------------------------
% Operator relating reduced external forces with external body forces
% --------------------------------------------------------------------
%\FskelDOM{e}{i}  &=     \BasisINTdom{ }{T} (    \J_{i} +  \Z_{i})  \Par{\fextVOLloc{e}{} + \fextTRACloc{e}{}}

% ---------------------------------------------------------
% Operator Z_1 and Z_2
% ---------------------------------------------------------------
%   \Z_{i} \defeq -\BasisUrb{e}{\face{i}} \Par{(\BasisUrbT{e}{f}  \BasisUrb{e}{f})^{-1}  \BasisUrbT{e}{}  }
% -------------------------------------------------------------------------

M = DATA_REFMESH.M ; % Mass matrix (Geometric)


%Beta_0_old = (BasisUrb(f,:)'*BasisUrb(f,:))\BasisUrb'  ; % This will be used for obtaining rRB
gRB = BasisUrb(f,:)'*(M(f,f)*BasisUrb(f,:)) ;
DATAOUT.gRB = gRB ;
Beta_0 =gRB\BasisUrb'  ; % This will be used for obtaining rRB


% VZ_1_old = (-V'*BasisUrb(f1,:)*Beta_0);
% VZ_2_old = (-V'*BasisUrb(f2,:)*Beta_0);
if ~iscell(V)
    Vnew{1} = V;
    Vnew{2} = V;
    V = Vnew ;
end

fGLO = {f1,f2} ;
VZ = {} ;
Vrotated = cell(size(V)) ;
for iface = 1:length(fGLO)
    % Change, 18-Jan-2019 (to improve efficiency)
    %     Z = -M(fGLO{iface},fGLO{iface})*BasisUrb(fGLO{iface},:)*Beta_0 ;
    Rotation_Intf_domain = DATA_REFMESH.RotationMatrixFace{iface} ;
    %     VZ{iface} =  RotateMatricesProduct_interface(Rotation_Intf_domain,Z,V{iface}) ;
    Vrotated  = RotateMatrix(Rotation_Intf_domain,V{iface}) ;
    VZ{iface} = -Vrotated'*M(fGLO{iface},fGLO{iface})*BasisUrb(fGLO{iface},:)*Beta_0;
end

% VZ_1 = (-V{1}'*M(f1,f1)*BasisUrb(f1,:)*Beta_0);
% VZ_2 = (-V{2}'*M(f2,f2)*BasisUrb(f2,:)*Beta_0);


%%
% Operators J_1 and J_2
%   \J_{i}  \defeq \BasisRdef{e}{\face{i}} \Par{ \Kbeam{e}{} \Hqr{e^T}\KdomRED{e^{-1}}{} \Par{\BasisUdefT{e}{}  -\BasisUdefT{e}{f} \BasisUrb{e}{f} (\BasisUrbT{e}{f}  \BasisUrb{e}{f})^{-1}  \BasisUrbT{e}{} } }

% J_i = J_Ai*JBi
% %
J_A = Kbeam*(Hqr'*inv(KdomRED))  ;
%J_A1 = BasisRdef(f1,:)*Kbeam*(Hqr'*inv(KdomRED)) ;
%J_B_old = BasisUdef'-  BasisUdef(f,:)'*(BasisUrb(f,:)*inv(BasisUrb(f,:)'*BasisUrb(f,:))*BasisUrb' ;

J_B = BasisUdef'-  (BasisUdef(f,:)'*M(f,f))*BasisUrb(f,:)*inv(BasisUrb(f,:)'*M(f,f)*BasisUrb(f,:))*BasisUrb' ;
Jbar = J_A*J_B ;

% Rotation matrices
% -------------------
% Code before rotation matrices were introduced
% ---------------------------------
% VJ_1 = (V{1}'*BasisRdef(f1,:))*Jbar;
% VJ_2 = (V{2}'*BasisRdef(f2,:))*Jbar ;

% New code
VJ = cell(1,2);
for iface = 1:length(fGLO)
    VJ{iface} = Tcomp{iface}*Jbar;
end

% ----------------------
% Product by V
Y = {} ;
for iface = 1:length(fGLO)
    Y{iface} = sparse((VJ{iface} - VZ{iface})) ;
end


% ---------------------------------------------
%% BODY FORCES PROPORTIONAL TO DENSITY
% ----------------------------------------
%% operator relating F_ext_vol of the reduced order model with its full-order model counterpart
GAMMAforce{1} = zeros(size(V{1},2),ndim) ;
GAMMAforce{2} = zeros(size(V{2},2),ndim) ;
%
dens = DATA_REFMESH.density ;
ngaus = size(DATA_REFMESH.posgp_RHS,2) ;
densGAUSS = zeros(ngaus*length(dens),1) ;
for igaus = 1:ngaus
    densGAUSS(igaus:ngaus:end) = dens ;
end

disp('Retrieving Nst')
tic
load(DATAIN.NAME_WS_MODES,'Ndom','WdiagRHS') ;
toc
disp('Done')

NstW =  WdiagRHS*Ndom ;

Lgravity = zeros(size(Kbeam,1),ndim) ;
Dgravity = zeros(size(J_B,1),ndim) ;

BETA_gravity = zeros(size(Beta_0,1),ndim) ;



for idim = 1:ndim
    GAMMAforce{1}(:,idim) =Y{1}*(NstW(idim:ndim:end,:)'*densGAUSS);
    GAMMAforce{2}(:,idim) =Y{2}*(NstW(idim:ndim:end,:)'*densGAUSS);
    Lgravity(:,idim) = Jbar*(NstW(idim:ndim:end,:)'*densGAUSS) ;
    Dgravity(:,idim) = J_B*(NstW(idim:ndim:end,:)'*densGAUSS) ;
    
    BETA_gravity(:,idim) = Beta_0*(NstW(idim:ndim:end,:)'*densGAUSS) ;
end

DATAOUT.GAMMAforce = GAMMAforce;
DATAOUT.Lgravity = Lgravity ;
DATAOUT.Dgravity = Dgravity ;

DATAOUT.BETA_gravity = BETA_gravity ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TRACTION FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------
Nbnd = DATA_REFMESH.Nbnd;
Wbnd = DATA_REFMESH.Wbnd ;
nfaces_contact = 2;
nfaces_lateral = length(Nbnd)-nfaces_contact ;

UPSILONforceTRAC.OPER = cell(2,length(Nbnd)) ;  % Operators
UPSILONforceTRAC.COOR = cell(1,length(Nbnd)) ;  % Coordinate nodal points of the face
UPSILONforceTRAC.NORMALS = cell(1,length(Nbnd)) ;  % Coordinate nodal points of the face

Q_tractions = cell(1,length(Nbnd)) ;
D_tractions = cell(1,length(Nbnd)) ;

BETA_tractions = cell(1,length(Nbnd)) ;


VALUES_ON_GAUSS_POINTS = 1;  % If = 1, it assumes forces are given at the Gauss points of the boundary elements


for iface = nfaces_contact+1:length(Nbnd)
    NbndFACE = DATA_REFMESH.Nbnd{iface} ;  % ngausBND*nelemBND x nnodeDOM
    WbndFACE = DATA_REFMESH.Wbnd{iface} ;  % ngausBND*nelemBND x ngausBND*nelemBND
    % Initialization
    
    if VALUES_ON_GAUSS_POINTS == 0
        iREF = 2 ;
    else
        iREF = 1;
    end
    OPER_1 = sparse(size(V{1},2),size(NbndFACE,iREF)*ndim) ;
    OPER_2 = sparse(size(V{2},2),size(NbndFACE,iREF)*ndim) ;
    
    OPER_Q = sparse(size(Kbeam,1),size(NbndFACE,iREF)*ndim) ;
    OPER_D = sparse(size(J_B,1),size(NbndFACE,iREF)*ndim) ;
    
    OPER_BETA = sparse(size(Beta_0,1),size(NbndFACE,iREF)*ndim) ;
    
    
    
    % List of connectivities
    % This is computed in PropertiesReferenceDomain.m
    CNb = DATA_REFMESH.CONNECTb{iface} ;
    
    
    
    % Coordinates face nodes
    NODESface = unique(CNb(:)); % Nodes pertaining to the surface "iface"
    if VALUES_ON_GAUSS_POINTS == 0
        UPSILONforceTRAC.COOR{iface} = DATA_REFMESH.COOR(NODESface,:) ;  % Coordinates of nodes pertaining to this face
    else
        %  UPSILONforceTRAC.COOR{iface} =[] ;
        % We have to compute the centroid of each element
        % -----------------------------------------------
        UPSILONforceTRAC.COOR{iface} = CentroidBoundary(DATA_REFMESH.COOR,CNb)  ;
        
        
        UPSILONforceTRAC.NORMALS{iface} = NormalsBoundary(DATA_REFMESH.COOR,CNb)  ;
    end
    % Assembly of operator FACE 1
    
    
    for idim = 1:ndim
        
        if VALUES_ON_GAUSS_POINTS == 0
            OPER_1(:,idim:ndim:end) =   Y{1}(:,idim:ndim:end)*(WbndFACE*NbndFACE)'*NbndFACE ;
            OPER_2(:,idim:ndim:end) =   Y{2}(:,idim:ndim:end)*(WbndFACE*NbndFACE)'*NbndFACE ;
            OPER_Q(:,idim:ndim:end) =   Jbar(:,idim:ndim:end)*(WbndFACE*NbndFACE)'*NbndFACE ;
            OPER_D(:,idim:ndim:end) =   J_B(:,idim:ndim:end)*(WbndFACE*NbndFACE)'*NbndFACE ;
            OPER_BETA(:,idim:ndim:end) =   Beta_0(:,idim:ndim:end)*(WbndFACE*NbndFACE)'*NbndFACE ;
        else
            OPER_1(:,idim:ndim:end) =   Y{1}(:,idim:ndim:end)*(WbndFACE*NbndFACE)' ;
            OPER_2(:,idim:ndim:end) =   Y{2}(:,idim:ndim:end)*(WbndFACE*NbndFACE)' ;
            OPER_Q(:,idim:ndim:end) =   Jbar(:,idim:ndim:end)*(WbndFACE*NbndFACE)' ;
            OPER_D(:,idim:ndim:end) =   J_B(:,idim:ndim:end)*(WbndFACE*NbndFACE)' ;
            OPER_BETA(:,idim:ndim:end) =   Beta_0(:,idim:ndim:end)*(WbndFACE*NbndFACE)' ;
        end
    end
    
    if VALUES_ON_GAUSS_POINTS == 0
        
        DOFsFACE = small2large(NODESface,ndim) ;
        
        
    else
        DOFsFACE = 1:size(OPER_1,2) ;
        
    end
    
    UPSILONforceTRAC.OPER{1,iface} =  OPER_1(:,DOFsFACE) ;
    UPSILONforceTRAC.OPER{2,iface} =  OPER_2(:,DOFsFACE) ;
    
    
    Q_tractions{iface} =  OPER_Q(:,DOFsFACE) ;
    D_tractions{iface} =  OPER_D(:,DOFsFACE) ;
    BETA_tractions{iface} =  OPER_BETA(:,DOFsFACE) ;
    
    %     % VALIDATION
    %     VALIDATION = 0;
    %     if VALIDATION == 1
    %
    %         Tnod = [0,1000,0]';
    %         Tnod = repmat(Tnod,1,size(NODESface,1));
    %         Tnod = Tnod(:) ;
    %         F_1 = UPSILONforceTRAC.OPER{1,iface}*Tnod;
    %
    %
    % %         FORCE = zeros(size(NbndFACE,2)*ndim,1) ;
    % %         M = (WbndFACE*NbndFACE)'*NbndFACE ;
    % %         M = M(:,NODESface) ;
    % %         for idim = 1:ndim
    % %             FORCE(idim:ndim:end) = M*Tnod(idim:ndim:end) ;
    % %         end
    % %
    % %         W_1 = BasisUdef'*FORCE ;
    % %         W_2 = BasisUdef(f,:)'*BasisUrb(f,:)*inv(BasisUrb(f,:)'*BasisUrb(f,:))*(BasisUrb'*FORCE) ;
    % %
    % %         F_e_ast = W_1-W_2 ;
    % %
    % %         R_e_ast = inv(BasisUdef'*BasisRdef)*F_e_ast ;
    % %
    % %         V_Psi_1 = V'*BasisRdef(f1,:) ;
    % %          V_Psi_2 = V'*BasisRdef(f2,:) ;
    % %
    % %          F_1 = V_Psi_1*R_e_ast ;
    % %          F_2 = V_Psi_2*R_e_ast ;
    %
    %     end
    
end

% ---------------

DATAOUT.BETA_tractions = BETA_tractions;
DATAOUT.UPSILONforceTRAC = UPSILONforceTRAC ;
DATAOUT.Q_tractions = Q_tractions ;
DATAOUT.D_tractions = D_tractions ;


end
