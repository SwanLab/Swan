function DATAOUT = ShapeFunCoefficients(COORnodes,ORDER_POLYNOMIALS)

 % Scaling (required to avoid that P is badly scaled)
% ireferenceCOOR = 1;
% DATAOUT.coorREF  = COORnodes(ireferenceCOOR,:) ;  % Coordinates reference point
DATAOUT.coorREF  = sum(COORnodes,1)/size(COORnodes,1) ;  % Introduced 14-Nov-2022. It amends
% issues when the order of the polynomial is high


% Now let us calculate the difference between all the
% coordinates and the reference coordinates
COORnodes_REF = bsxfun(@minus,COORnodes',DATAOUT.coorREF')';
DATAOUT.LelemSCALING = max(abs(COORnodes_REF)) ;
 
% Transformed coordinates (reltive to coorREF), and scaled by   LelemSCALING
for idim = 1:length(DATAOUT.LelemSCALING)
    COORnodes(:,idim) = COORnodes_REF(:,idim)/DATAOUT.LelemSCALING(idim);
end
% Coefficientes of  the polynomial
P = CoordinateMatrixPolynomial(COORnodes,ORDER_POLYNOMIALS)  ;
DATAOUT.COEFFSpol =  inv(P) ; % P\eye(nnodeE);  % Coefficients of the polynomials