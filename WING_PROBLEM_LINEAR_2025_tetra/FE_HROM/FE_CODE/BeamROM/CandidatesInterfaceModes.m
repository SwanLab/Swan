function [nBOUNDARY_INTFMODES,BasisREFERENCE,Vrb_Morth] =...
    CandidatesInterfaceModes(BasisUdef_L,f1,f2,DISP_BOUND,SingVal_Udef_L,...
    DATAIN,BasisUdef_B,M,Vrb,BasisUrb)

if nargin ==0
    load('tmp2.mat')
end

% %
%  WEIGH_BY_SINGULAR_VALUES = 1 ;
% 
%  if WEIGH_BY_SINGULAR_VALUES ==1
%      SingVal_Udef_L = SingVal_Udef_L/SingVal_Udef_L(1) ;
%      BasisUdef_L =    bsxfun(@times,BasisUdef_L',SingVal_Udef_L)' ;
%  end

% nRB = size(Vrb,2) ; 
% if nRB == 3
%     ndimSP = 2; 
% else
%     ndimSP = 3; 
% end

nfaces = 2;  

BasisINTFdefCAND = cell(nfaces,1) ;
%BasisRdom = cell(nfaces,1) ;
BasisUdom = cell(nfaces,1) ;
%    BasisUdom{1} = [BasisUrb(f1,:),SingVal_Udef_L(f1,:)] ; 
%     BasisUdom{2} = [BasisUrb(f2,:),SingVal_Udef_L(f2,:)] ; 

DATAIN = DefaultField(DATAIN,'ROTATION_LOC',[]) ; 
R = DATAIN.ROTATION_LOC ;
if isempty(R)
    BasisREFERENCE = [BasisUdef_L(f1,:), BasisUdef_L(f2,:)];
else
    % Displacements of face 2 are given in domain coordinates.
    % Hence, they should be transformed back to local coordinates of the
    % interfae ---by multiplying by the transpose of the rotation matrix
    BasisUdef_Lf2= RotateMatrix(R',BasisUdef_L(f2,:))  ;
    BasisREFERENCE = [BasisUdef_L(f1,:), BasisUdef_Lf2];
     BasisUdom{2} =  RotateMatrix(R', BasisUdom{2})  ;
    
end

%%% PURGE ROTATION MODES as WELL AS FLUCTUATION MODES
% ----------------------------------------------------
% Matrix of fluctuation modes
% .---------------------------
%PurgeFluctuations = 0;
% if  PurgeFluctuations == 1
%     BasisBeamInterface = [BasisUdef_B(f1,:), BasisUdef_B(f2,:)];
%     PG = (Vrb'*M*Vrb)  ;
%     BasisBeamInterface = BasisBeamInterface  - Vrb*(PG\(Vrb'*M*BasisBeamInterface)) ;
%     TOL = 1e-6 ;
%     [BeamFluctuationModes,S,Vbar] = SVDT(BasisBeamInterface,TOL) ;
%     BeamModesInterface = [Vrb,BeamFluctuationModes]  ;
% else
    BeamModesInterface = Vrb;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we make the candidates for interface modes M-orthogonal to
% matrix BeamModesInterface 
Mchol = chol(M) ;
nmodesINTF = size(BasisUdef_L,2) ;
BasisREFERENCE = MorthogonalMatrix(M,Mchol,BasisREFERENCE,BeamModesInterface)  ;
%BasisREFERENCE = BasisREFERENCE(:,1:nmodesINTF) ;  It is not advisable to 
%truncate at this level (9-Dec-2018)
%%%%%%%%%%%%%%%
% Now we turn BeamModesInterface (rigid body) M-orthogonal
Xbar = Mchol*BeamModesInterface ; 
TOL = 0; 
[Ubar,S,Vbar] = SVDT(Xbar,TOL) ;
Vrb_Morth = Mchol\Ubar ;




% INTERFACE MODES for boundaries
nBOUNDARY_INTFMODES = [] ;
if ~isempty(DISP_BOUND.MODES)
    % The first modes to try are those corresponding to boundary
    % displacements
    MODES_DISP = [] ;
    nBOUNDARY_INTFMODES = size(DISP_BOUND.MODES,2) ;
    for imodes = 1:nBOUNDARY_INTFMODES
        if DISP_BOUND.FACES == 1
            fREF = f1;
        else
            fREF  = f2;
        end
        MODES_DISP = [MODES_DISP,DISP_BOUND.MODES(fREF,imodes)] ;
    end
    % Now we make MODES_DISP M-orthogonal with respect to rigid body
    % motions ..................
    MODES_BOUNDARY = MorthogonalMatrix(M,Mchol,MODES_DISP,Vrb)  ;
    
    % Finally, we make BasisREFERENCE M-orthogonal to  MODES_BOUNDARY
    MODES_INTERIOR = MorthogonalMatrix(M,Mchol,BasisREFERENCE,MODES_BOUNDARY)  ;
    
    %     %
    %          [MODES_BOUNDARY,SS,VbarS] = SVDT(MODES_DISP,0) ;
    %          MODES_INTERIOR = BasisREFERENCE - MODES_BOUNDARY*(MODES_BOUNDARY'*BasisREFERENCE);
    %
    %
    %     [MODES_INTERIOR,S,Vbar] = SVDT(MODES_INTERIOR,0) ;
    
    
    BasisREFERENCE = [MODES_BOUNDARY,MODES_INTERIOR] ;
end
