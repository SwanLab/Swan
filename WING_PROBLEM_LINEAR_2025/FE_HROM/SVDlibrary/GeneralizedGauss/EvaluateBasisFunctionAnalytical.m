function  [FUNk dFUNm,DATAOUT_irreg]= ...
    EvaluateBasisFunctionAnalytical(xNEW,VSinv,DATALOC,EVALUATE_GRADIENT,EVALUATE_FUN)

if nargin == 0
    load('tmp2.mat')
end
% NAME_FUNCTION_TO_INTEGRATE
%---------------------------------
% 1) [FUNk,FUNk_x,FUNk_y]   =    LagrangePolynomial2D(DATALOC.xLIM,DATALOC.PORDER,xNEW) ;

% 2) [FUNk, FUNk_x,FUNk_y,DATAOUT_irreg ]=    LagrangePolynomial2D_B_B(DATALOC.xLIM,DATALOC.PORDER,xNEW,DATALOC)

% 3) [FUNk,FUNk_x,FUNk_y,FUNk_z]   =    LagrangePolynomial3D(DATALOC.xLIM,DATALOC.PORDER,xNEW,DATALOC) ;

% 4)  [FUNk, dFUN,DATAOUT_irreg ]=    LagrangePolynomial_ALLD_B_B(DATALOC.xLIM,DATALOC.PORDER,xNEW,DATALOC)  ;

%       if iscell(FUNk)
%               FUNk = cell2mat(FUNk) ;
%               FUNk_x = cell2mat(dFUN(1,:)) ;   FUNk_y = cell2mat(dFUN(2,:)) ;    FUNk_z = cell2mat(dFUN(3,:));
%           else
%               FUNk_x = dFUN{1} ;  FUNk_y = dFUN{2} ;   FUNk_z = dFUN{3} ;
%           end
%

NAME_FUNCTION_TO_INTEGRATE = DATALOC.NAME_FUNCTION_TO_INTEGRATE ;
PHIk_y = [] ;
dPHIk_y ={} ;
DATAOUT_irreg = [] ;
DATALOC.EVALUATE_GRADIENT = EVALUATE_GRADIENT;
DATALOC.EVALUATE_FUN = EVALUATE_FUN;

[FUNk,dFUN,DATAOUTfun]= feval(NAME_FUNCTION_TO_INTEGRATE,DATALOC.xLIM,DATALOC.PORDER,xNEW,DATALOC) ;

ndim = size(xNEW,2) ;

if iscell(FUNk)
    FUNk = cell2mat(FUNk) ;
    dFUNm = cell(1,ndim) ;
    for idim =1:ndim
        dFUNm{idim} = cell2mat(dFUN(idim,:)) ;
    end
    
else
    dFUNm = dFUN ;
    
    
end





if EVALUATE_FUN == 1
    
    if ~isempty(VSinv)
        FUNk =FUNk*VSinv ;
    end
    
end

if EVALUATE_GRADIENT == 1
    if ~isempty(VSinv)
        for idim = 1:ndim
            dFUNm{idim} = dFUNm{idim}*VSinv ;
        end
    end
    
    
end


