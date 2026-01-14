function BasisINT = SelectCandidatesInterfaceModes_BEAMS(Vrb,M,BasisRdef_L,BasisREFERENCE,...
    TOL_SINGULAR_VALUES_Hqr,f1,f2,R,DATAIN)




%
PG = (Vrb'*M*Vrb)  ;
%


MODES_INCLUDE =[] ;
IMODES_INCLUDE = [] ;
imode = 1;

% The number of interface modes cannot be greater than the number of
% reaction modes
while  length(IMODES_INCLUDE) <size(BasisRdef_L,2)
 %   while  imode <size(BasisRdef_L,2)

    %    if USE_REACTIONS == 1
    newINTFmode = BasisREFERENCE(:,imode)  - Vrb*(PG\(Vrb'*M*BasisREFERENCE(:,imode))) ;
    %   else
    %     newINTFmode = BasisRdef_L(f1,imode)  - Vrb*(PG\(Vrb'*M*BasisRdef_L(f1,imode))) ;
    %  end
    NEW_MODES = [MODES_INCLUDE,newINTFmode ] ;
    COV_f1 = NEW_MODES'*BasisRdef_L(f1,:);
    if isempty(R)
        COV_f2 = NEW_MODES'*BasisRdef_L(f2,:);
    else
        NEW_MODES_R= RotateMatrix(R,NEW_MODES)  ;
        COV_f2 = NEW_MODES_R'*BasisRdef_L(f2,:);
    end
    SSVAL_f1 = svd(COV_f1) ;
    SSVAL_f2 = svd(COV_f2) ;
    if imode == 1
        %  if SSVAL_f1(1) > 1e-10 &
        ratioSV_f1 = 1  ;
        ratioSV_f2 = 1  ;
        % end
    else
        ratioSV_f1 = SSVAL_f1(end)/SSVAL_f1(end-1) ;
        ratioSV_f2 = SSVAL_f2(end)/SSVAL_f2(end-1) ;
    end
    if ratioSV_f1 >= TOL_SINGULAR_VALUES_Hqr && ratioSV_f2 >= TOL_SINGULAR_VALUES_Hqr
        MODES_INCLUDE = NEW_MODES ;
        IMODES_INCLUDE(end+1) = imode ;
    end
    imode = imode + 1;
    
    
    end


DATAIN = DefaultField(DATAIN,'nINTFadditional',size(MODES_INCLUDE,2))  ;  % Additional modes interfaces

if ~isempty(DATAIN.nINTFadditional)
    BasisINTdef = MODES_INCLUDE(:,1:DATAIN.nINTFadditional) ;
else
    BasisINTdef = MODES_INCLUDE  ;
end

BasisINT = [Vrb,BasisINTdef] ;