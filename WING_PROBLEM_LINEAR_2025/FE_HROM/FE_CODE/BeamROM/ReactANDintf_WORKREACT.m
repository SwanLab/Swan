function [BasisRdef,BasisINT,TEXTR,BasisUdef,SingVal_Udef] = ...
    ReactANDintf_WORKREACT(BasisUdef,BasisRdef,f1,f2,TOL_SINGULAR_VALUES_Hqr,...
    nBASES_BEAM,DATA_REFMESH,Vrb,M,DATAIN,SingVal_Udef,SingVal_Rdef,BasisUrb,Mdom)
if nargin == 0
    load('tmp1.mat')
end

DATALOC = DATAIN.INTERFACE_MODES_REACTIONS_WORK ;
DATALOC = DefaultField(DATALOC,'TOLERANCE_SVD_DEF_MODES',1e-6) ; % TOLERANCE FOR SVD DEFORMATIONAL MODES
DATALOC = DefaultField(DATALOC,'TOLERANCE_ANGLE_INTERSECTION_REACTIONS',0.01) ; % TOLERANCE FOR SVD DEFORMATIONAL MODES
DATALOC = DefaultField(DATALOC,'TOLERANCE_SVD_REACTION_MODES_INTERSECTION',1e-6) ; % TOLERANCE FOR SVD reaction MODES

f = [f1;f2] ;

% ----------------------------------------------
% Rotation
% ---------
if ~isempty(DATA_REFMESH.RotationMatrixFace)
    iface = 2;
    DATAIN.ROTATION_LOC = DATA_REFMESH.RotationMatrixFace{iface} ;
    R = DATAIN.ROTATION_LOC ;
else
    R = [] ;
end


% STEP1 ) Selection of displacement modes   --- by maximizing the alignment with reaction modes  
% ---------------------------------------------------------------------------------------
[BasisUdef,SingVal_Udef, TEXTR]= ReactionDispModesEqual(BasisUdef,BasisRdef,f,...
    SingVal_Udef) ;


% STEP2 ) Candidate for being interface modes
% -------------------------------------
% Vdef_Morth --> Rigid body modes (always included)
% Vrv_Morth --> Additional modes (DEFORMATIONAL MODES)
DATAIN.TOLERANCE_TRUNCATE_CANDIDATE_INTERFACE_MODES  = DATALOC.TOLERANCE_SVD_DEF_MODES ;
[Vdef_Morth,Vrb_Morth,TEXTR,BasisRdefROT] = CandidatesInterfaceModesW(BasisUdef,f1,f2,SingVal_Udef,...
    DATAIN,M,Vrb,TEXTR,BasisRdef ) ;


% STEP3) DECIDING HOW MANY MODES TO INCLUDE

METHOD_DECIDE_NUMBER_OF_MODES = 0 ;

if METHOD_DECIDE_NUMBER_OF_MODES == 1    
    % METHOD BASED ON THE TRUNCATED SINGULAR VALUE (ABANDONED IN FAVOR OF THE ONE SHOWN BELOW)    
    DATALOC = DefaultField(DATALOC,'TOL_MIN_SINGUL_VAL_REACTIONS_FORCES',0.1);
    TOL_SVD_INCLUDE = DATALOC.TOL_MIN_SINGUL_VAL_REACTIONS_FORCES ;
    S1 = svd(BasisRdefROT{1}) ;
    S1 = S1/S1(1) ;
    S2 = svd(BasisRdefROT{2}) ;
    S2 = S2/S2(1) ;
    nmodes1 = length(find(S1>TOL_SVD_INCLUDE)) ;
    nmodes2 = length(find(S2>TOL_SVD_INCLUDE)) ;
    nmodesR = min(nmodes1,nmodes2) ;    
    TEXTR{end+1} = ['Number of singular values of BasisRdef_1 over TOL = ',num2str(TOL_SVD_INCLUDE),' --> ',num2str(nmodes1)] ;
    TEXTR{end+1} = ['Number of singular values of BasisRdef_2 over TOL = ',num2str(TOL_SVD_INCLUDE),' --> ',num2str(nmodes2)] ;
       
else    
    % INTERSECTION OF SUBSPACES 
    % ------------------------------
    TOL_SVD =   DATALOC.TOLERANCE_SVD_REACTION_MODES_INTERSECTION; 
    TOL_ANGLE = DATALOC.TOLERANCE_ANGLE_INTERSECTION_REACTIONS ;
    [RR,ANGLES ]= IntersectionSubspaces(BasisRdefROT{1},BasisRdefROT{2},TOL_SVD,TOL_ANGLE) ;
    nmodesR = size(RR,2) ;
    
    figure(5623)
    hold on
    plot(ANGLES)
    xlabel('Modes')
    ylabel('Angles (degrees)')
    title('INTERSECTION BETWEEN SUBSPACES BasisRdef_{f1} and BasisRdef_{f2}')
 
    TEXTR{end+1} = ['Dimension of the intersection of reaction subspaces (for TOL ANGLE = ',num2str(TOL_ANGLE),') --> ',num2str(nmodesR)] ;
end


% STEP 4) CHOOSING INTERFACE MODES 
% --------------------------------
% Method based on picking up nmodesR linearly independent columns from
% BasisRdef'*U
[BasisINT,TEXTR] =  ChooseIntModesWORK1(Vrb_Morth,Vdef_Morth,SingVal_Rdef,BasisRdefROT,TEXTR,nmodesR,Vrb) ; 



%
% if DATAIN.KINEMATIC_CONSTRAINTS_MODES.DEIM_BASED_METHOD == 1
%     % May-14th- New method. Similar to the one implemented for RVEs. See
%     % InterfaceRVEKinematicConstraint.m
%     BasisUdom = [BasisUrb,BasisUdef] ;
%
%     [BasisINT,TEXTR ]= InterfaceRVEKinematicConstraint_BEAMS(Vrb,M,BasisRdef,Vdef_Morth,...
%         TOL_SINGULAR_VALUES_Hqr,f1,f2,R,DATAIN,Vrb_Morth,Mdom,BasisUdom,BasisUrb,TEXTR) ;
%
%
% else
%
%     if DATAIN.DeformationalInterfaceModes_AlignmentMethod == 1
%         warning('This implementation is not reliable')
%         DATAIN.RotationMatrixLocal{1} = [] ;
%         DATAIN.RotationMatrixLocal{2} = R' ;
%         iface = 1;
%         [BasisINT,~,TEXTR ]=  DeformModesInterface_AlignmentMethod(BasisRdef,f1,f2,...
%             Vrb,M,DATAIN,BasisUdef,SingVal_Udef,SingVal_Rdef,iface,TEXTR) ;
%        % TEXTB ={} ;
%     else
%
%         % --------
%         %% Final selection of interface modes
%         [BasisINT,TEXTR ]= SelectCandidatesInterfaceModes_BEAMSnew(Vrb,M,BasisRdef,Vdef_Morth,...
%             TOL_SINGULAR_VALUES_Hqr,f1,f2,R,DATAIN,Vrb_Morth) ;
%
%     end
%
% end
