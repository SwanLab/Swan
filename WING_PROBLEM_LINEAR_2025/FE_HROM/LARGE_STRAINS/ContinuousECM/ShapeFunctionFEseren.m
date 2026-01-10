function [N, dN,DATAOUT ]=    ShapeFunctionFEseren(xLIM,x,DATA)

if nargin ==0
    load('tmp.mat')
end

DATAOUT = [] ;


%COEFFSpol = ShapeFunCoefficients(COORnodes,ORDER_POLYNOMIALS) ;
%DATA.Integrand.COEFFSpol = COEFFSpol;
%DATA.Integrand.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
%DATA.Integrand.NameFunctionGenerate  = 'ShapeFunctionFE';
DATAshape = DATA.DATAshape;  
   x  = (x-DATAshape.coorREF)./DATAshape.LelemSCALING   ;  % Transformed coordinate of the point under study
        factorDERIVATIVE = 1./DATAshape.LelemSCALING ;    % This for computing the derivatives


[Pevaluate, PevaluateDER]= CoordinateMatrixPolynomial(x,DATA.ORDER_POLYNOMIALS)  ;

COLUMNS_REMAIN = DATA.DATAshape.COLUMNS_TO_REMAIN; 
Pevaluate = Pevaluate(:,COLUMNS_REMAIN) ; 
for i = 1:length(PevaluateDER)
    PevaluateDER{i} = PevaluateDER{i}(:,COLUMNS_REMAIN) ; 
end


N = Pevaluate*DATAshape.COEFFSpol ;  % Shape functions at the given points COORevaluate
%    PHIk_y(inew,:)  =   N*VAR_SMOOTH_FE.BasisIntegrand(NODESloc,:) ;
dN = cell(size(PevaluateDER)) ;
for idim = 1:length(PevaluateDER)
    dN{idim} =  factorDERIVATIVE(idim)*PevaluateDER{idim}*DATAshape.COEFFSpol ;
    %       dPHIk_y{idim}(inew,:) = factorDERIVATIVE(idim)*dN*VAR_SMOOTH_FE.BasisIntegrand(NODESloc,:) ;
end