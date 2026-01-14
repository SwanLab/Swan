function [NtN, dNtN,DATAOUT ]=    ShapeFunctionFE_NtN(xLIM,x,DATA)

if nargin ==0
    load('tmp.mat')
end

DATAOUT = [] ;

[N, dN,DATAOUT ]=    ShapeFunctionFE(xLIM,x,DATA) ; 

ndim = size(x,2) ; 
nfun = size(N,2) ; 
NtN = zeros(size(N,1),nfun^2) ; 
dNtN = cell(1,ndim) ; 
ifunacum = 0; 
for ifun = 1:nfun
    for jfun = 1:nfun
        ifunacum = ifunacum +1 ; 
        NtN(:,ifunacum) = N(:,ifun).*N(:,jfun) ; 
        for idim = 1:size(x,2)
            dNtN{idim}(:,ifunacum) = dN{idim}(:,ifun).*N(:,jfun) + N(:,ifun).*dN{idim}(:,jfun)  ; 
        end
    end
end

% 
% %COEFFSpol = ShapeFunCoefficients(COORnodes,ORDER_POLYNOMIALS) ;
% %DATA.Integrand.COEFFSpol = COEFFSpol;
% %DATA.Integrand.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
% %DATA.Integrand.NameFunctionGenerate  = 'ShapeFunctionFE';
% DATAshape = DATA.DATAshape;  
%    x  = (x-DATAshape.coorREF)./DATAshape.LelemSCALING   ;  % Transformed coordinate of the point under study
%         factorDERIVATIVE = 1./DATAshape.LelemSCALING ;    % This for computing the derivatives
% 
% 
% [Pevaluate, PevaluateDER]= CoordinateMatrixPolynomial(x,DATA.ORDER_POLYNOMIALS)  ;
% N = Pevaluate*DATAshape.COEFFSpol ;  % Shape functions at the given points COORevaluate
% %    PHIk_y(inew,:)  =   N*VAR_SMOOTH_FE.BasisIntegrand(NODESloc,:) ;
% dN = cell(size(PevaluateDER)) ;
% for idim = 1:length(PevaluateDER)
%     dN{idim} =  factorDERIVATIVE(idim)*PevaluateDER{idim}*DATAshape.COEFFSpol ;
%     %       dPHIk_y{idim}(inew,:) = factorDERIVATIVE(idim)*dN*VAR_SMOOTH_FE.BasisIntegrand(NODESloc,:) ;
% end