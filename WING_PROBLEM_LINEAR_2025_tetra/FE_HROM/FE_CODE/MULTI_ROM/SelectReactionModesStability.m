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
