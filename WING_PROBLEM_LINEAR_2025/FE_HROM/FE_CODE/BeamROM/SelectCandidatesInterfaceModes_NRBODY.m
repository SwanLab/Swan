function BasisINT = SelectCandidatesInterfaceModes_NRBODY(Vrb,M,BasisRdef,Vdef_Morth,...
    TOL_SINGULAR_VALUES_Hqr,f1,f2,R,DATAIN,Vrb_Morth)

if nargin == 0
    load('tmp.mat')
end


% We construct the set of candidates, including the rigid body modes
% (which has been orthogonalized with respect to M)
BasisREFERENCE = [Vdef_Morth] ;

%%%%%
MODES_INCLUDE =[] ; % rigid body vectors are always included
IMODES_INCLUDE =[];
imode = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATAIN  = DefaultField(DATAIN,'TOL_REL_MAX_SVD_FIRST_AND_END',TOL_SINGULAR_VALUES_Hqr) ;
TOL_REL_MAX_SV = DATAIN.TOL_REL_MAX_SVD_FIRST_AND_END ;
TOL_REL_MAX_SV = 1e-4
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
        
        error('The set of FE training tests is not sufficient to generate beam modes')
    end
else
    COV_f1
    COV_f2
    error('The set of FE training tests is not sufficient to generate beam modes')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The number of interface modes cannot be greater than the number of
% reaction modes
while  length(IMODES_INCLUDE) <size(BasisRdef,2) && imode <=size(BasisREFERENCE,2)
    
    newINTFmode = BasisREFERENCE(:,imode)  ;
    NEW_MODES = [MODES_INCLUDE,newINTFmode ] ;
    COV_f1 = NEW_MODES'*BasisRdef(f1,:);
    if isempty(R)
        COV_f2 = NEW_MODES'*BasisRdef(f2,:);
    else
        NEW_MODES_R= RotateMatrix(R,NEW_MODES)  ;
        COV_f2 = NEW_MODES_R'*BasisRdef(f2,:);
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
    if (ratioSV_f1 >= TOL_SINGULAR_VALUES_Hqr && ratioSV_f2 >= TOL_SINGULAR_VALUES_Hqr) ...
            && length(SSVAL_f1) == sNEWMODES && length(SSVAL_f2) == sNEWMODES
        if isempty(TOL_REL_MAX_SV)
            MODES_INCLUDE = NEW_MODES ;
            IMODES_INCLUDE(end+1) = imode ;
        else
            ratioSV_ABS_1 = SSVAL_f1(end)/SSVAL_f1(1) ;
            ratioSV_ABS_2 = SSVAL_f2(end)/SSVAL_f2(1) ;
            if ratioSV_ABS_1 >= TOL_REL_MAX_SV && ratioSV_ABS_2 >= TOL_REL_MAX_SV
                MODES_INCLUDE = NEW_MODES ;
                IMODES_INCLUDE(end+1) = imode ;
            end
        end
    end
    imode = imode + 1;
    
    
end


% DATAIN = DefaultField(DATAIN,'nINTFadditional',size(MODES_INCLUDE,2))  ;  % Additional modes interfaces
%
% if ~isempty(DATAIN.nINTFadditional)
%     BasisINTdef = MODES_INCLUDE(:,1:DATAIN.nINTFadditional) ;
% else
BasisINT = MODES_INCLUDE  ;
%end

%BasisINT = [Vrb,BasisINTdef(:,size(Vrb,2)+1:end)] ;