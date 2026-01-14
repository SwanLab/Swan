function [BasisINT,TEXTB ]= SelectCandidatesInterfaceModes_BEAMSnew(Vrb,M,BasisRdef,Vdef_Morth,...
    TOL_SINGULAR_VALUES_Hqr,f1,f2,R,DATAIN,Vrb_Morth)

if nargin == 0
    load('tmp.mat')
end


% We construct the set of candidates, including the rigid body modes
% (which has been orthogonalized with respect to M)
BasisREFERENCE = [Vrb_Morth,Vdef_Morth] ;

%%%%%
MODES_INCLUDE =Vrb_Morth ; % rigid body vectors are always included
IMODES_INCLUDE = [1:size(Vrb,2)] ;
imode = size(Vrb,2)+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOL_REL_MAX_SV = 1e-4;

if DATAIN.UnifiedApproachForModes == 1
    % We have to check that axial, shear and bending modes are present in
    % the basis matrices
    COV_f1 = Vrb_Morth'*BasisRdef(f1,:);
    if isempty(R)
        COV_f2 = Vrb_Morth'*BasisRdef(f2,:);
    else
        NEW_MODES_R= RotateMatrix(R,Vrb_Morth)  ;
        COV_f2 = NEW_MODES_R'*BasisRdef(f2,:);
    end
    SSVAL_f1 = svd(COV_f1) ;
    SSVAL_f2 = svd(COV_f2) ;
    
    if  length(SSVAL_f1) >= size(Vrb_Morth,2) &&  length(SSVAL_f2) >= size(Vrb_Morth,2)
        ratioSV_ABS_1 = SSVAL_f1(end)/SSVAL_f1(1) ;
        ratioSV_ABS_2 = SSVAL_f2(end)/SSVAL_f2(1) ;
        if ratioSV_ABS_1 <= TOL_REL_MAX_SV || ratioSV_ABS_2 <= TOL_REL_MAX_SV
            
            error('The set of FE training tests may  not be sufficient to generate beam modes')
        end
    else
        COV_f1
        COV_f2
        error('The set of FE training tests is not sufficient to generate beam modes')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATAIN  = DefaultField(DATAIN,'TOL_REL_MAX_SVD_FIRST_AND_END',TOL_SINGULAR_VALUES_Hqr*0.1) ;
TOL_REL_MAX_SV = DATAIN.TOL_REL_MAX_SVD_FIRST_AND_END; %DATAIN.TOL_REL_MAX_SVD_FIRST_AND_END ;
% The number of interface modes cannot be greater than the number of
% reaction modes
%DATAIN = DefaultField(DATAIN,'FILTERING_INF_MODES_MINIMUM_CORRELATION_ON_JUST_ONE_FACE',0);



TEXTB  = {} ;
RATIOS_END{1} = [] ;
RATIOS_END{2} = [] ;

DATAIN.CriterioJumpSVD_end_firstVAlues = 2 ;

DATAIN = DefaultField(DATAIN,'FILTER_OUT_REACTION_MODES_RATIO_SINGULAR_VALUES',[]) ;


% Filtering reaction modes 
[BasisRdef1,BasisRdef2] = FilterBeforeReactionModesInterface(DATAIN,BasisRdef,f1,f2) ; 
% ----------------------


while  length(IMODES_INCLUDE) <size(BasisRdef,2) && imode <=size(BasisREFERENCE,2)
    
    newINTFmode = BasisREFERENCE(:,imode)  ;
    NEW_MODES = [MODES_INCLUDE,newINTFmode ] ;
    COV_f1 = NEW_MODES'*BasisRdef1;
    if isempty(R)
        COV_f2 = NEW_MODES'*BasisRdef2;
    else
        NEW_MODES_R= RotateMatrix(R,NEW_MODES)  ;
        COV_f2 = NEW_MODES_R'*BasisRdef2;
    end
    SSVAL_f1 = svd(COV_f1) ;
    SSVAL_f2 = svd(COV_f2) ;
    sNEWMODES = min(size(NEW_MODES)) ;
    if imode == 1
        %  if SSVAL_f1(1) > 1e-10 &
        ratioSV_f1 = 1  ;
        ratioSV_f2 = 1  ;
        % end
    else
        ratioSV_f1 = SSVAL_f1(end)/SSVAL_f1(end-1) ;
        ratioSV_f2 = SSVAL_f2(end)/SSVAL_f2(end-1) ;
    end
    [TEXTB,MODES_INCLUDE,IMODES_INCLUDE,rr1,rr2] = FilteringLocalModes(ratioSV_f1,TOL_SINGULAR_VALUES_Hqr,ratioSV_f2,SSVAL_f1,SSVAL_f2,sNEWMODES,TOL_REL_MAX_SV,...
        NEW_MODES,imode,IMODES_INCLUDE,TEXTB,MODES_INCLUDE,DATAIN)  ;
    
    if rr1 ~= 0
        RATIOS_END{1}(end+1) =  rr1 ;
    end
    if rr1 ~= 0
        RATIOS_END{2}(end+1) =  rr2 ;
    end
    
    imode = imode + 1;
    
    
end
TEXTB{end+1} =['Selected interface modes = ',num2str(IMODES_INCLUDE)] ;

TEXTB{end+1} ='*********************************************************************' ;
TEXTB{end+1} ='*********************************************************************' ;

TEXTB{end+1} =['RATIOS svd FACE 1:   ',num2str(RATIOS_END{1})]
TEXTB{end+1} =['RATIOS svd FACE 2:  ',num2str(RATIOS_END{2})]


% figure(1970)
% hold on
% xlabel('IMODE')
% ylabel('RATIO')
% hh=plot(RATIOS_END{1})
% hhh= plot(RATIOS_END{2},'r')
% DATAIN = DefaultField(DATAIN,'nINTFadditional',size(MODES_INCLUDE,2))  ;  % Additional modes interfaces
%
% if ~isempty(DATAIN.nINTFadditional)
%     BasisINTdef = MODES_INCLUDE(:,1:DATAIN.nINTFadditional) ;
% else
BasisINTdef = MODES_INCLUDE  ;
%end

BasisINT = [Vrb,BasisINTdef(:,size(Vrb,2)+1:end)] ;