function [nBOUNDARY_INTFMODES,BasisREFERENCE] = CandidatesInterfaceModes_New(BasisUdef_L,f1,f2,DISP_BOUND,SingVal_Udef_L,...
    DATAIN,BasisUdef_B,M,Vrb)

if nargin ==0
    load('tmp2.mat')
end

%
%  WEIGH_BY_SINGULAR_VALUES = 1 ;
%
%  if WEIGH_BY_SINGULAR_VALUES ==1
%      SingVal_Udef_L = SingVal_Udef_L/SingVal_Udef_L(1) ;
%      BasisUdef_L =    bsxfun(@times,BasisUdef_L',SingVal_Udef_L)' ;
%  end

BasisREFERENCE = [BasisUdef_L(f1,:), BasisUdef_L(f2,:)];

%%% PURGE ROTATION MODES as WELL AS FLUCTUATION MODES
% ----------------------------------------------------
% Matrix of fluctuation modes
% .---------------------------
PurgeFluctuations = 0 ;
if  PurgeFluctuations == 1
    BasisBeamInterface = [BasisUdef_B(f1,:), BasisUdef_B(f2,:)];
    PG = (Vrb'*M*Vrb)  ;
    BasisBeamInterface = BasisBeamInterface  - Vrb*(PG\(Vrb'*M*BasisBeamInterface)) ;
    TOL = 1e-6 ;
    [BeamFluctuationModes,S,Vbar] = SVDT(BasisBeamInterface,TOL) ;
    BeamModesInterface = [Vrb,BeamFluctuationModes]  ;
else
    BeamModesInterface = Vrb;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we make the candidates for interface modes M-orthogonal to
% matrix BeamModesInterface
Mchol = chol(M) ;
nmodesINTF = size(BasisUdef_L,2) ;
BasisREFERENCE = MorthogonalMatrix(M,Mchol,BasisREFERENCE,BeamModesInterface)  ;
BasisREFERENCE = BasisREFERENCE(:,1:nmodesINTF) ;
%%%%%%%%%%%%%%%


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
    % Now we make MODES_DISP M-orthogonal to BeamModesInterface and BasisREFERENCE
    BeamInterior_rb = [BeamModesInterface,BasisREFERENCE] ;
    MODES_BOUNDARY = MorthogonalMatrix(M,Mchol,MODES_DISP,BeamInterior_rb)  ;
    
    %
    %     [MODES_BOUNDARY,SS,VbarS] = SVDT(MODES_DISP,0) ;
    %     MODES_INTERIOR = BasisREFERENCE - MODES_BOUNDARY*(MODES_BOUNDARY'*BasisREFERENCE);
    %     [MODES_INTERIOR,S,Vbar] = SVDT(MODES_INTERIOR,0) ;
    
    
    BasisREFERENCE = [MODES_BOUNDARY,BasisREFERENCE] ;
end
