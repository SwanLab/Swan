function [DATAOUT]= ExtForcesROMrve_operators(BasisUrb,fI,f,BasisRdef,Kbeam,Hqr,KdomRED,...
    DATA_REFMESH,BasisUdef,V,ndim,DATAOUT,DATAIN)

% Copy of function ExtForcesROM_operators (for beams)

% .GAMMAforceVOL,DATAOUT.UPSILONforceTRAC,...
%   DATAOUT.Lgravity,DATAOUT.Q_tractions


if nargin == 0
    load('tmp1.mat')
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
%Beta_0 = (BasisUrb(f,:)'*M(f,f)*BasisUrb(f,:))\BasisUrb'  ; % This will be used for obtaining rRB
gRB = BasisUrb(f,:)'*(M(f,f)*BasisUrb(f,:)) ; 
DATAOUT.gRB = gRB ; 
Beta_0 =gRB\BasisUrb'  ; % This will be used for obtaining rRB  . Non-linear cases. 13-Jun-2019

%Beta_0 = (BasisUrb(f,:)'*M(f,f)*BasisUrb(f,:))\BasisUrb'  ;  % Correction 16-sept-2018

if ~isstruct(V)
    % OLD IMPLEMENTATION (SLICES)
    % VZ_1 = (-V{1}'*M(f1,f1)*BasisUrb(f1,:)*Beta_0);
    % VZ_2 = (-V{2}'*M(f2,f2)*BasisUrb(f2,:)*Beta_0);
    % NEW  IMPLEMENTATION
    
    
   
    
    
    VZ = cell(size(V)) ;
    Vrotated = cell(size(V)) ; 
    for ifaceLOC=1:length(V)
        Rotation_Intf_domain = DATA_REFMESH.RotationMatrixFace{ifaceLOC} ;
      %  Vrotated  = RotateMatrix(Rotation_Intf_domain,V{ifaceLOC}) ;
      %  Change 3-July-2019
      Vrotated{ifaceLOC}  = RotateMatrix(Rotation_Intf_domain,V{ifaceLOC}) ;
      
     %   VZ{ifaceLOC} = -V{ifaceLOC}'*M(fI{ifaceLOC},fI{ifaceLOC})*BasisUrb(fI{ifaceLOC},:)*Beta_0 ;
    
     %     VZ{ifaceLOC} =
       %     -Vrotated'*M(fI{ifaceLOC},fI{ifaceLOC})*BasisUrb(fI{ifaceLOC},:)*Beta_0;
       %     % 3-July-2019
            VZ{ifaceLOC} = -Vrotated{ifaceLOC}'*M(fI{ifaceLOC},fI{ifaceLOC})*BasisUrb(fI{ifaceLOC},:)*Beta_0;
            % 3-July-2019
       

    end
else
    % Kinematically constrained method . From 26-Apr-2019
    % No longer  used...
    error('This method proved to be unreliable')
    VZ = -V.BasisINTF'*M(f,f)*BasisUrb(f,:)*Beta_0 ;
end


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

if ~isstruct(V)
    VJ = cell(size(V)) ;
    %VJ_1 = (V{1}'*BasisRdef(f1,:))*Jbar;
    %VJ_2 = (V{2}'*BasisRdef(f2,:))*Jbar ;
    
    % MODIFICATION 3-JULY-2019
    % ----------------------------------- ROTATION OF V
    
    for iface = 1:length(V)
        %  Rotation_Intf_domain = DATA_REFMESH.RotationMatrixFace{ifaceLOC} ;
     %   VJ{iface} = (V{iface}'*BasisRdef(fI{iface},:))*Jbar;  % Change 3-July-2019
        
        VJ{iface} = (Vrotated{iface}'*BasisRdef(fI{iface},:))*Jbar;  % Change 3-July-2019
    end
    % ----------------------
    % Product by V
    Y = cell(size(V));
    % Y_1 = sparse((VJ_1 - VZ_1)) ;
    % Y_2 = sparse(( VJ_2 - VZ_2 ))  ;
    for iface = 1:length(V)
        Y{iface} = sparse(VJ{iface}-VZ{iface}) ;
    end
else
    % New method , constrained kinematically, Apr-2019
    VJ = V.BasisINTF'*BasisRdef(f,:)*Jbar ;
    Y = VJ-VZ ;
    nfaces = length(V.BasisINTFall_cell) ; 
    
end


% % What if we ignore J
% Y_1 = sparse(V'*( - Z_1)) ;
% Y_2 = sparse(V'*(  - Z_2 ))  ;

% ---------------------------------------------
%% BODY FORCES PROPORTIONAL TO DENSITY
% ----------------------------------------
%% operator relating F_ext_vol of the reduced order model with its full-order model counterpart
% GAMMAforce{1} = zeros(size(V{1},2),ndim) ;
% GAMMAforce{2} = zeros(size(V{2},2),ndim) ;
if ~isstruct(V)
    nfaces =length(V) ;
    for iface= 1:nfaces
        GAMMAforce{iface} = zeros(size(V{iface},2),ndim);
    end
else
    % New version, apr-2019
    GAMMAforce = zeros(size(V.BasisINTF,2),ndim); 
end
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
    if ~isstruct(V)
    for iface = 1:nfaces
        GAMMAforce{iface}(:,idim) =Y{iface}*(NstW(idim:ndim:end,:)'*densGAUSS);
        %  GAMMAforce{2}(:,idim) =Y_2*(NstW(idim:ndim:end,:)'*densGAUSS);
    end
    else
        GAMMAforce(:,idim) =Y*(NstW(idim:ndim:end,:)'*densGAUSS);
    end
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
nfaces_contact = nfaces;
nfaces_lateral = length(Nbnd)-nfaces_contact ;

UPSILONforceTRAC.OPER = cell(2,length(Nbnd)) ;  % Operators
UPSILONforceTRAC.COOR = cell(1,length(Nbnd)) ;  % Coordinate nodal points of the face
UPSILONforceTRAC.NORMALS = cell(1,length(Nbnd)) ;  % Coordinate nodal points of the face

Q_tractions = cell(1,length(Nbnd)) ;
D_tractions = cell(1,length(Nbnd)) ;

BETA_tractions = cell(1,length(Nbnd)) ;


VALUES_ON_GAUSS_POINTS = 1;  % If = 1, it assumes forces are given at the Gauss points of the boundary elements


for iface = nfaces_contact+1:length(Nbnd)
    
%    error('Revise this part for the new approach ..... Apr-2019')
    
    NbndFACE = DATA_REFMESH.Nbnd{iface} ;  % ngausBND*nelemBND x nnodeDOM
    WbndFACE = DATA_REFMESH.Wbnd{iface} ;  % ngausBND*nelemBND x ngausBND*nelemBND
    % Initialization
    
    if VALUES_ON_GAUSS_POINTS == 0
        iREF = 2 ;
    else
        iREF = 1;
    end
    OPER = cell(size(V));
    for ifaceLOC=1:nfaces
        OPER{ifaceLOC} = sparse(size(V{ifaceLOC},2),size(NbndFACE,iREF)*ndim) ;
        %  OPER_2 = sparse(size(V{2},2),size(NbndFACE,iREF)*ndim) ;
    end
    
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
            error('Option not implemented')
            %             OPER_1(:,idim:ndim:end) =   Y_1(:,idim:ndim:end)*(WbndFACE*NbndFACE)'*NbndFACE ;
            %             OPER_2(:,idim:ndim:end) =   Y_2(:,idim:ndim:end)*(WbndFACE*NbndFACE)'*NbndFACE ;
            %             OPER_Q(:,idim:ndim:end) =   Jbar(:,idim:ndim:end)*(WbndFACE*NbndFACE)'*NbndFACE ;
            %             OPER_D(:,idim:ndim:end) =   J_B(:,idim:ndim:end)*(WbndFACE*NbndFACE)'*NbndFACE ;
            %             OPER_BETA(:,idim:ndim:end) =   Beta_0(:,idim:ndim:end)*(WbndFACE*NbndFACE)'*NbndFACE ;
        else
            for ifaceLOC=1:nfaces
                OPER{ifaceLOC}(:,idim:ndim:end) =   Y{ifaceLOC}(:,idim:ndim:end)*(WbndFACE*NbndFACE)' ;
                %     OPER_2(:,idim:ndim:end) =   Y_2(:,idim:ndim:end)*(WbndFACE*NbndFACE)' ;
            end
            OPER_Q(:,idim:ndim:end) =   Jbar(:,idim:ndim:end)*(WbndFACE*NbndFACE)' ;
            OPER_D(:,idim:ndim:end) =   J_B(:,idim:ndim:end)*(WbndFACE*NbndFACE)' ;
            OPER_BETA(:,idim:ndim:end) =   Beta_0(:,idim:ndim:end)*(WbndFACE*NbndFACE)' ;
        end
    end
    
    if VALUES_ON_GAUSS_POINTS == 0
        
        DOFsFACE = small2large(NODESface,ndim) ;
        
        
    else
        DOFsFACE = 1:size(OPER{1},2) ;
        
    end
    for ifaceLOC = 1:nfaces
        UPSILONforceTRAC.OPER{ifaceLOC,iface} =  OPER{ifaceLOC}(:,DOFsFACE) ;
        %    UPSILONforceTRAC.OPER{2,iface} =  OPER_2(:,DOFsFACE) ;
    end
    
    
    
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
