function [COEFFSpol,LelemSCALING,coorREF,POLYINFO] = CoeffsPolyShapeFunctions(elemCONTAINER,NODESloc,VAR_SMOOTH_FE,POLYINFO)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/ContinuousECM/EvaluateBasisFunctionDIRECTFIT.m
if nargin == 0
    load('tmp1.mat')
end

% INDEXES ASSOCIATED GAUSS POINTS OF THE ELEMENT
if ~isempty(POLYINFO.COEFFSpolynomial{elemCONTAINER})
    % THE COEFFICIENTS HAVE BEEN ALREADY CALCULATED
    COEFFSpol = POLYINFO.COEFFSpolynomial{elemCONTAINER} ;
    LelemSCALING = POLYINFO.SCALING_VARIABLES.LENGTH(elemCONTAINER,:) ;
    coorREF = POLYINFO.SCALING_VARIABLES.coorREF(elemCONTAINER,:) ;
else
    % cOORDINATES OF THE GAUSS POINTS OF THE ELEMENT CONTAINING THE POINT
    % UNDER CONSIDERATION
    COORelement = VAR_SMOOTH_FE.COORg(NODESloc,:) ; %   coordinates of the Gauss points
    if length(COORelement) == 1
        error(['The number of Gauss points per element is to be higher than one'])
    end 
    % Scaling (required to avoid that P is badly scaled)
    ireferenceCOOR = 1;
    coorREF  = COORelement(ireferenceCOOR,:) ;  % Coordinates reference point
    % Now let us calculate the difference between all the
    % coordinates and the reference coordinates
    COORelement_REF = bsxfun(@minus,COORelement',coorREF')';
    LelemSCALING = max(abs(COORelement_REF)) ;
    
    POLYINFO.SCALING_VARIABLES.LENGTH(elemCONTAINER,:) = LelemSCALING ;
    POLYINFO.SCALING_VARIABLES.coorREF(elemCONTAINER,:) = coorREF ;
    % Transformed coordinates (reltive to coorREF), and scaled by   LelemSCALING
    for idim = 1:length(LelemSCALING)
        COORelement(:,idim) = COORelement_REF(:,idim)/LelemSCALING(idim);
    end
    % Coefficientes of  the polynomial
    P = CoordinateMatrixPolynomial(COORelement,VAR_SMOOTH_FE.ORDER_POLYNOMIALS)  ;
    
    COEFFSpol =  inv(P) ; % P\eye(nnodeE);  % Coefficients of the polynomials
    POLYINFO.COEFFSpolynomial{elemCONTAINER} = COEFFSpol ;
end