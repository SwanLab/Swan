function BasisRdef_L = ReactionModesLocalEffects(BasisUdef_L,BasisRdef_L,f,TOL_SINGULAR_VALUES_Hqr)



 
% BasisUdef_Lf1 = BasisUdef_L(f_reference,:) ;
% BasisRdef_Lf1 = BasisRdef_L(f_reference,:) ;

BasisUdef_Lf = BasisUdef_L(f,:) ;
BasisRdef_Lf = BasisRdef_L(f,:) ;

% Determining number of reaction modes
%b -----------------------------------
imode = 1;

MODES_INCLUDE =[] ;
IMODES_INCLUDE = [] ;
while  imode <=size(BasisRdef_Lf,2)
    NEW_MODES = [MODES_INCLUDE,BasisRdef_Lf(:,imode) ] ;
    HqrT = NEW_MODES'*BasisUdef_Lf;
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
BasisRdef_L = BasisRdef_L(:,IMODES_INCLUDE);
