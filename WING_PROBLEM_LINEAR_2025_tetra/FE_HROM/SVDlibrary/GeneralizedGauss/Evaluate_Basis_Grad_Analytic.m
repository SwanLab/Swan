function  [g,Grad_g]=  Evaluate_Basis_Grad_Analytic(xNEW,VSinv,DATALOC,EVALUATE_GRADIENT,DATAFITTING)

if nargin == 0
    load('tmp.mat')
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


% IF DATALOC.Evaluate_Fun_Gradients_via_FITTING ==1, then both the function and its gradient are calculated
% numerically (via splines)


DATALOC = DefaultField(DATALOC,'Evaluate_Fun_Gradients_via_FITTING',0);

if DATALOC.Evaluate_Fun_Gradients_via_FITTING == 0
    NAME_FUNCTION_TO_INTEGRATE = DATALOC.TYPEFUN ;
    DATALOC.EVALUATE_GRADIENT = EVALUATE_GRADIENT;
    [f,Grad_f,~]= feval(NAME_FUNCTION_TO_INTEGRATE,DATALOC.xLIM,xNEW,DATALOC) ;
    
    
    ndim = size(xNEW,2) ;
    
    if iscell(f)
        f = cell2mat(f) ;
        Grad_fm = cell(1,ndim) ;
        for idim =1:ndim
            Grad_fm{idim} = cell2mat(Grad_f(idim,:)) ;
        end
    else
        Grad_fm = Grad_f ;
    end
    
    
    
    
    
    if ~isempty(VSinv)
        g =(f*VSinv)' ;
    else
        g = f' ;
    end
    
    Grad_g = cell(size(Grad_fm)) ;
    if EVALUATE_GRADIENT == 1
        if ~isempty(VSinv)
            for idim = 1:ndim
                Grad_g{idim} = (Grad_fm{idim}*VSinv)' ;
            end
        else
            for idim = 1:ndim
                Grad_g{idim} = dFUNm{idim}' ;
            end
        end
        
        
    end
elseif DATALOC.Evaluate_Fun_Gradients_via_FITTING == 1
    [g,Grad_g] = FUN_GRAD_approximate(xNEW,DATAFITTING);
else
    error('Option not implemented')
    
end

