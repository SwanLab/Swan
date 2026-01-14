function   DATAOUT = RVEStiffnessMatrix_PLATES(BASES,DATA_REFMESH,DATAIN)
% Copy of RVEStiffnessMatrix.m
% Version for plate structures
% JAHO, 5-Sept-2018
% ------------------------------
if nargin == 0
    load('tmp1.mat')
end

ndim = 3;
% Deformational displacement modes
% -----------------------------------
BasisUdef = BASES.DISPLACEMENTS.U ;  % Basis displacements
SinvVal_Udef = BASES.DISPLACEMENTS.S ; % aSSOCIATED SINGULAR VALUES
DATAIN = DefaultField(DATAIN,'nMODES_DISP',size(BasisUdef,2)) ;
BasisUdef = BasisUdef(:,1:DATAIN.nMODES_DISP)  ;
DATAOUT.BasisUdef = BasisUdef ;
% Self-equilibrated displacement modes
% -----------------------------------
BasisRdef = BASES.REACTIONS.U ;  % Basis displacements
DATAIN = DefaultField(DATAIN,'nMODES_REACTIONS',size(BasisRdef,2)) ;
BasisRdef = BasisRdef(:,1:DATAIN.nMODES_REACTIONS)  ;



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
NODES_faces_all = cell2mat(DATA_REFMESH.NODES_faces) ;
% Degrees of freedom faces f1 and f2, f3 and f4
for i=1:length(DATA_REFMESH.NODES_faces)
    fI{i} = small2large(DATA_REFMESH.NODES_faces{i},ndim) ;
    DATAOUT.fI{i} =fI{i};
end
% _______________________________________________________-
% Boundary DOFs: c1,c2,c3,c4,s1,s2,s3,s4
% ----------------------------------------
NODES_corners_all = cell2mat(DATA_REFMESH.NODES_CORNERS) ;
% Degrees of freedom corners c1 and c2, c3 and f4
for i=1:length(DATA_REFMESH.NODES_CORNERS)
    cornerI{i} = small2large(DATA_REFMESH.NODES_CORNERS{i},ndim) ;
    DATAOUT.cornerI{i} =cornerI{i};
end
% Boundary DOFs: c1,c2,c3,c4,s1,s2,s3,s4
% ----------------------------------------
NODES_sides_all = cell2mat(DATA_REFMESH.NODES_SIDES) ;
% Degrees of freedom sides s1 and s2, s3 and s4
for i=1:length(DATA_REFMESH.NODES_SIDES)
    sideI{i} = small2large(DATA_REFMESH.NODES_SIDES{i},ndim) ;
    DATAOUT.sideI{i} =sideI{i};
end

% Boundary DOFS
NODES_BND =  [NODES_corners_all;NODES_sides_all];
f = small2large(NODES_BND,ndim);
DATAOUT.f = f;
%-------------------------------------------
% Rigid body modes for FACES
Vrb = DATA_REFMESH.RigidBodyModesInterface ;
DATAOUT.BasisIntRB = Vrb ;
DATAIN = DefaultField(DATAIN,'RIGID_BODY_MODES_INTERFACE',1) ;
%%% Displacements interface boundaries
%DISP_BOUND = BASES.DISPLACEMENTS_INTERFACES_BOUNDARY ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selection of reaction modes to meet stability conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATAIN = DefaultField(DATAIN,'TOL_SINGULAR_VALUES_Hqr',0.1) ;
BasisRdef = SelectReactionModesStability(BasisUdef,BasisRdef,DATAOUT.f,DATAIN.TOL_SINGULAR_VALUES_Hqr) ;
DATAOUT.BasisRdef = BasisRdef ;

%%%%%



% if  DATAIN.RIGID_BODY_MODES_INTERFACE ==1
%     BasisIntBeam = Vrb ;
% else
%    % DATAOUT.BasisInt = BasisUdef(f1,:) ;
%     BasisIntBeam =  BasisUdef(f1,:) ;
% end

% ------------------------
% Geometric mass matrices
% ------------------------
nfaces = length(Vrb);
M1d = DATA_REFMESH.GeometricMassMatrixInterface ;
M = cell(nfaces,1) ;
for iface= 1:nfaces
    M{iface} = sparse(size(M1d{iface},1)*ndim,size(M1d{iface},2)*ndim) ;
    for idim = 1:ndim
        M{iface}(idim:ndim:end,idim:ndim:end) = M1d{iface} ;
    end
end
nBASES_RVE =[] ;

if isempty(nBASES_RVE)
    % if  DATAIN.RIGID_BODY_MODES_INTERFACE ==1
    % DATAOUT.BasisInt = DATA_REFMESH.BasisINTfc ;
    BasisINTfc = DATA_REFMESH.BasisINTfc ;  % Basis matrix
    %else
    % V = []
    %   error('This option does not work')
    %  V =  [BasisUdef(f2,:) BasisUdef(f1,:)]; %- Vrb*((Vrb'*M*Vrb)\(Vrb'*M*BasisUdef(f1,:)))  ;
    %   [V S]= svd(V) ;
    %    DATAOUT.BasisInt = V ;
    % end
else
    error('Option not implemented yet')
    % if  DATAIN.RIGID_BODY_MODES_INTERFACE ==1
    %
    % SO FAR, this module is only prepared for slices. Accordingly,
    % BasisINTbeam{1} = BasisINTbeam{2}
    %     iface= 1;
    %     BasisINTbeam = Vrb{iface} ;
    %     %else
    %     %    BasisINTbeam = BasisUdef(f2,1:nBASES_RVE.DISPLACEMENTS) ;
    %     %   BasisINTbeam=  BasisINTbeam - Vrb*((Vrb'*M*Vrb)\(Vrb'*M*BasisINTbeam))  ;
    %     %   BasisINTbeam = orth(BasisINTbeam) ;
    %     %  end
    %     % Choose reaction modes and interface modes
    %
    %     %         [BasisRdef,V,nBOUNDARY_INTFMODES] = ReactionAndInterfaceLocalModes(BasisUdef,BasisRdef,f1,f2,DATAIN.TOL_SINGULAR_VALUES_Hqr,...
    %     %          nBASES_RVE,DATA_REFMESH,BasisINTbeam,M,DISP_BOUND,DATAIN)  ;
    %
    %     [BasisRdef,V,nBOUNDARY_INTFMODES] = ReactionAndInterfaceLocalModes_new(BasisUdef,BasisRdef,f1,f2,DATAIN.TOL_SINGULAR_VALUES_Hqr,...
    %         nBASES_RVE,DATA_REFMESH,BasisINTbeam,M{iface},DISP_BOUND,DATAIN,SinvVal_Udef)  ;
    %     disp(['Number of additional displacement modes =',num2str(size(BasisUdef,2)-nBASES_RVE.DISPLACEMENTS)])
    %     disp(['Number of additional reaction modes =',num2str(size(BasisRdef,2)-nBASES_RVE.REACTIONS)])
    %     disp(['Number of additional interface modes =',num2str(size(V,2)-6)])
    %     disp(['Number of additional interface modes (boundary) =',num2str(nBOUNDARY_INTFMODES)])
    %
    %     DATAOUT.nBOUNDARY_INTFMODES = nBOUNDARY_INTFMODES ;
    %     DATAOUT.BasisInt = V ;
    %     DATAOUT.BasisRdef = BasisRdef ;
end
% -------------------------------------------------------------
KdomRED = BasisUdef'*(K*BasisUdef) ;  % Reduced stiffness matrix
BasisRdef_f = BasisRdef(f,:) ;   % REduced-stiffness at the boundary
BasisUdef_f = BasisUdef(f,:) ;
%  \Hqr{e} =  \BasisUdefT{e}{f}  \BasisRdef{e}{f}
Hqr = BasisUdef_f'*BasisRdef_f ;
if rank(Hqr) < size(BasisRdef_f,2)
    error('Ill-conditioned matrix. Probably the training set has rank less than 18')
end
DATAOUT.Hqr = Hqr ;
DATAOUT.KdomRED = KdomRED ;
%  \Kbeam{e}{} \defeq (\Hqr{e^T}\KdomRED{e^{-1}}{} \Hqr{e})^{-1},
Kbeam_inv = Hqr'*(KdomRED\Hqr) ;
Kbeam = inv(Kbeam_inv) ;
DATAOUT.Kbeam = Kbeam ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Expression for Tcomp
% --------------------
%  \Tcomp{e^T} \defeq       \BasisRdefT{e}{f}    \BasisINTdom{fc}{e}
Tcomp = BasisINTfc'*BasisRdef_f ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduced stiffness matrix
% \KskelDOM{e}{}  = \Tcomp{e} \Kbeam{e}{} \Tcomp{e^T}
Kskel = Tcomp*Kbeam*Tcomp' ;

% CHECK_IMPLE = 0;
% if CHECK_IMPLE == 1
%     % Checl 1 
%     BENDING_INDICES{1} = [3,4,5]';
%     for iii= 2:4
%         BENDING_INDICES{iii} = BENDING_INDICES{iii-1} + 5 ;
%     end
%     BENDING_INDICES = cell2mat(BENDING_INDICES') ;
%     Kbending = Kskel(BENDING_INDICES,BENDING_INDICES) ;
%     
%     Kbending = mat2cell(Kbending,[3,3,3,3],[3,3,3,3]) ;
%     
%     %%   check 2
%     a = zeros(20,1) ;
%     IND = [11,16] ;
%     u = 1;
%     a(IND) = u ;
%     
%     F = Kskel*a ;
%     
%     % Check 3 
%     % -------
%     DOFr = [1:5 6] ; 
%     IND = [11,16] ;
%     DOFl = setdiff(1:20,DOFr) ; 
%     F = zeros(20,1) ; 
%     F(IND) = 1e3; 
%     a = zeros(20,1) ; 
%     a(DOFl) = Kskel(DOFl,DOFl)\F(DOFl) ; 
%     
% end



DATAOUT.Kskel  = Kskel ;
% ------------------------------------------------------
nBASES_RVE.REACTIONS = size(BasisRdef,2) ;
nBASES_RVE.DISPLACEMENTS = size(BasisUdef,2) ;
% % Stiffness constant matrix
% BeamPropertiesStiffness(BasisRdef(:,1:nBASES_RVE.REACTIONS),f1,f2,Vrb,L,AREA,DATA_REFMESH,DATAIN,...
%     BasisUdef(:,1:nBASES_RVE.DISPLACEMENTS),K) ;
% -----------------------------------------------------
% Operators relating external forces in the ROM and the full-order model
[DATAOUT]= ExtForcesROMrve_operPLATES(BasisUrb,fI,f,BasisRdef,Kbeam,Hqr,KdomRED,DATA_REFMESH,BasisUdef,BasisINTfc,ndim,DATAOUT,DATAIN) ;
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

end


function BasisRdef = SelectReactionModesStability(BasisUdef,BasisRdef,f,TOL_SINGULAR_VALUES_Hqr)



% BasisUdef_Lf1 = BasisUdef_L(f_reference,:) ;
% BasisRdef_Lf1 = BasisRdef_L(f_reference,:) ;

BasisUdef_f = BasisUdef(f,:) ;
BasisRdef_f = BasisRdef(f,:) ;

% Determining number of reaction modes
%b -----------------------------------
imode = 1;

MODES_INCLUDE =[] ;
IMODES_INCLUDE = [] ;
while imode <=size(BasisRdef_f,2) &&  length(IMODES_INCLUDE) <size(BasisUdef_f,2)   % It cannot be greater than nU
    NEW_MODES = [MODES_INCLUDE,BasisRdef_f(:,imode) ] ;
    HqrT = NEW_MODES'*BasisUdef_f;
    SSVAL = svd(HqrT) ;
    if imode == 1
        ratioSV =1 ; % SSVAL(1) ;
    else
        ratioSV = SSVAL(end)/SSVAL(end-1) ;
    end
    if ratioSV >= TOL_SINGULAR_VALUES_Hqr
        MODES_INCLUDE = NEW_MODES ;
        IMODES_INCLUDE(end+1) = imode ;
    end
    imode = imode + 1;
    
    
end

nmodesR = imode-1;
BasisRdef = BasisRdef(:,IMODES_INCLUDE);

% DATAOUT.BasisRdef = BasisRdef ;

end
