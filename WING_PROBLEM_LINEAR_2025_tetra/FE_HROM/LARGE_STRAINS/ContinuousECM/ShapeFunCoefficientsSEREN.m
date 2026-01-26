function DATAOUT = ShapeFunCoefficientsSEREN(COORnodes,ORDER_POLYNOMIALS)
% Serendipity elements
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/03_20nodeHEX.mlx
if nargin == 0
    load('tmp.mat')
end

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

if size(P,1) == 20
    % Columns to be eliminated 
%     000  001  002     010  011  012 -  020  021  022     ** 1:9   Eliminate  -->  9
% 100  101  002     110  011  012 -  120  121  122   ** 10:18: Eliminate--> 18 
% 200  201  202     210  211  212 -  220  221  222    :19:27:  Eliminate -->  21, 24, 25,26,27
colelim = [9,18,21,24,25,26,27] ; 
colrema = setdiff(1:27,colelim) ; 
P = P(:,colrema) ; 

 
elseif size(P,1) == 8 
    % COLUMNS TO BE ELIMINATED 
    % --------------------------
       % Columns to be eliminated 
% 00  01  02     10  11 12    20  21  22  --- > 9 
colelim = 9 ; 
colrema = setdiff(1:9,colelim) ; 
P = P(:,colrema) ; 
    
    
else
    error('Option not implemented')
end

DATAOUT.COEFFSpol =  inv(P) ; % P\eye(nnodeE);  % Coefficients of the polynomials
DATAOUT.IS_SEREND = 1; 
DATAOUT.COLUMNS_TO_REMAIN = colrema ; 