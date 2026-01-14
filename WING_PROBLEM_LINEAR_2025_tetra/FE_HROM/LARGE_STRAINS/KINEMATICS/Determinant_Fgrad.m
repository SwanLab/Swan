function detF = Determinant_Fgrad(FgradST,ndim) ; 
% Determinat of FgradST 
% See symbolic function  MATLAB_CODES/FE_HROM/LARGE_STRAINS/SYMBOLIC/Determinant_F.m
% JAHO, 10-dec-2020


if ndim == 2
    nF = 4 ;
    nelem_ngaus = length(FgradST)/nF  ;
    
    FROWS = cell(1,nF) ;
    for icols  =1:nF
        FROWS{icols} = icols:nF:size(FgradST,1) ;
    end
    
    % Generated automatically by  MATLAB_CODES/FE_HROM/LARGE_STRAINS/SYMBOLIC/Determinant_F.m
detF =FgradST(FROWS{1}).*FgradST(FROWS{2}) - FgradST(FROWS{3}).*FgradST(FROWS{4});

    
else
    
    
     nF = 9 ;
    nelem_ngaus = length(FgradST)/nF  ;
    
    FROWS = cell(1,nF) ;
    for icols  =1:nF
        FROWS{icols} = icols:nF:size(FgradST,1) ;
    end
        % Generated automatically by  MATLAB_CODES/FE_HROM/LARGE_STRAINS/SYMBOLIC/Determinant_F.m

        detF =FgradST(FROWS{1}).*FgradST(FROWS{2}).*FgradST(FROWS{3}) - FgradST(FROWS{1}).*FgradST(FROWS{4}).*FgradST(FROWS{7}) - FgradST(FROWS{2}).*FgradST(FROWS{5}).*FgradST(FROWS{8}) - FgradST(FROWS{3}).*FgradST(FROWS{6}).*FgradST(FROWS{9}) + FgradST(FROWS{4}).*FgradST(FROWS{6}).*FgradST(FROWS{8}) + FgradST(FROWS{5}).*FgradST(FROWS{7}).*FgradST(FROWS{9});

end