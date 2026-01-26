function  [PHIk_y dPHIk_y,DATAOUT_irreg]= ...
    EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,EVALUATE_GRADIENT,EVALUATE_FUN)

if nargin == 0
    load('tmp2.mat')
end

DATAOUT_irreg = [] ; 
DATALOC = DefaultField(DATALOC,'FUN_GRADIENT_DIRECTLY',[]) ; 
DATALOC.FUN_GRADIENT_DIRECTLY = DefaultField(DATALOC.FUN_GRADIENT_DIRECTLY,'NAME_FUN',[]) ; 

NAMEFUN = DATALOC.FUN_GRADIENT_DIRECTLY.NAME_FUN  ;
NAME_GRAD = DATALOC.FUN_GRADIENT_DIRECTLY.NAME_FUN ;
PHIk_y = [] ;
dPHIk_y ={} ;
  DATALOC.EVALUATE_GRADIENT = EVALUATE_GRADIENT; 
switch  NAMEFUN
    case 'Poly2Dfun'
        if EVALUATE_FUN == 1
            FUNk = Poly2Dfun(DATALOC.PORDER,xNEW(:,1),xNEW(:,2)) ;
            
            %  PHIk_y =FUNk_y*VSinv ;
        end
        if EVALUATE_GRADIENT == 1
            [FUNk_x  FUNk_y]= dPoly2Dfun(DATALOC.PORDER,xNEW(:,1),xNEW(:,2)) ;
            %
            %             dPHIk_y{1}  =FUNk_x*VSinv ;
            %             dPHIk_y{2}  =FUNk_y*VSinv ;
        end
        
    case 'LagrangePOLY2D'   
        
  %      if EVALUATE_FUN == 1
   %         FUNk   =    LagrangePolynomial2D(DATALOC.xLIM,DATALOC.PORDER,xNEW) ; % Poly2Dfun(DATALOC.PORDER,xNEW(:,1),xNEW(:,2)) ;
            
            %  PHIk_y =FUNk_y*VSinv ;
    %    end
     %   if EVALUATE_GRADIENT == 1
    %       [FUNk,FUNk_x,FUNk_y]   =    LagrangePolynomial2D(DATALOC.xLIM,DATALOC.PORDER,xNEW) ;
           
           if ~isempty(DATALOC.IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA)
               % See LagrangePolynomial2D_irregular_aux.mlx
                [FUNk,FUNk_x,FUNk_y,DATAOUT_irreg] =    LagrangePolynomial2D_irregular(DATALOC.xLIM,DATALOC.PORDER,xNEW,DATALOC) ;
           else
                [FUNk,FUNk_x,FUNk_y]   =    LagrangePolynomial2D(DATALOC.xLIM,DATALOC.PORDER,xNEW) ;
           end
           
      %  end
      
    case 'LagrangePOLY2D_B_B'
        
        % Terms int(B^T*B*detJ)
        %------------------------
      
          [FUNk, FUNk_x,FUNk_y,DATAOUT_irreg ]=    LagrangePolynomial2D_B_B(DATALOC.xLIM,DATALOC.PORDER,xNEW,DATALOC)  ; 
          
     case 'LagrangePOLY2D_N_N'
        
        % Terms int(B^T*B*detJ)
        %------------------------
       
          [FUNk, FUNk_x,FUNk_y,DATAOUT_irreg ]=    LagrangePolynomial2D_N_N(DATALOC.xLIM,DATALOC.PORDER,xNEW,DATALOC)  ;      
      
       case 'LagrangePOLY3D' 
          
         [FUNk,dFUN]   =    LagrangePolynomial3D(DATALOC.xLIM,DATALOC.PORDER,xNEW,DATALOC) ;
         
         FUNk_x = dFUN{1} ;  FUNk_y = dFUN{2} ;    FUNk_z = dFUN{3} ; 
         
           case 'LagrangePOLY3D_B_B'
        
        % Terms int(B^T*B*detJ)
        %------------------------
      
          [FUNk, dFUN,DATAOUT_irreg ]=    LagrangePolynomial_ALLD_B_B(DATALOC.xLIM,DATALOC.PORDER,xNEW,DATALOC)  ;
          
          if iscell(FUNk)
              FUNk = cell2mat(FUNk) ; 
              FUNk_x = cell2mat(dFUN(1,:)) ;   FUNk_y = cell2mat(dFUN(2,:)) ;    FUNk_z = cell2mat(dFUN(3,:));
          else
              FUNk_x = dFUN{1} ;  FUNk_y = dFUN{2} ;   FUNk_z = dFUN{3} ;
          end
         
    case 'Poly3Dfun'
        if EVALUATE_FUN == 1
            FUNk = Poly3Dfun(DATALOC.PORDER,xNEW(:,1),xNEW(:,2),xNEW(:,3)) ;
            
            %  PHIk_y =FUNk_y*VSinv ;
        end
        if EVALUATE_GRADIENT == 1
            [FUNk_x  FUNk_y FUNk_z]= dPoly3Dfun(DATALOC.PORDER,xNEW(:,1),xNEW(:,2),xNEW(:,3)) ;
            
            %             dPHIk_y{1}  =FUNk_x*VSinv ;
            %             dPHIk_y{2}  =FUNk_y*VSinv ;
            %             dPHIk_y{3}  =FUNk_z*VSinv ;
        end
    case 'Poly3Dfun_const'
        if EVALUATE_FUN == 1
            FUNk = Poly3Dfun_const(DATALOC.PORDER,xNEW(:,1),xNEW(:,2),xNEW(:,3)) ;
            
            %   PHIk_y =FUNk_y*VSinv ;
        end
        if EVALUATE_GRADIENT == 1
            [FUNk_x  FUNk_y FUNk_z]= dPoly3Dfun(DATALOC.PORDER,xNEW(:,1),xNEW(:,2),xNEW(:,3)) ;
            
            %             dPHIk_y{1}  =FUNk_x*VSinv ;
            %             dPHIk_y{2}  =FUNk_y*VSinv ;
            %             dPHIk_y{3}  =FUNk_z*VSinv ;
        end
    case 'SinCosExpFun'
        %  dbstop('52')
        if EVALUATE_FUN == 1
            FUNk = SinCosExpFun(DATALOC.mu,xNEW(:,1),xNEW(:,2)) ;
            
            % PHIk_y =FUNk_y*VSinv ;
        end
        if EVALUATE_GRADIENT == 1
            [FUNk_x  FUNk_y ]= dSinCosExpFun(DATALOC.mu,xNEW(:,1),xNEW(:,2)) ;
            
            %             dPHIk_y{1}  =FUNk_x*VSinv ;
            %             dPHIk_y{2}  =FUNk_y*VSinv ;
            
        end
        
    case 'SinCosExpFun3D'
        %  dbstop('52')
        if EVALUATE_FUN == 1
            FUNk = SinCosExpFun3D(DATALOC.mu,xNEW(:,1),xNEW(:,2),xNEW(:,3)) ;
            
            
        end
        if EVALUATE_GRADIENT == 1
            [FUNk_x  FUNk_y  FUNk_z]= dSinCosExpFun3D(DATALOC.mu,xNEW(:,1),xNEW(:,2),xNEW(:,3)) ;
            
            
        end
        
    case 'SinCosExpFun3D_2param'
        %  dbstop('52')
        if EVALUATE_FUN == 1
            FUNk = SinCosExpFun3D_2param(DATALOC.mu,xNEW(:,1),xNEW(:,2),xNEW(:,3)) ;
            
            
        end
        if EVALUATE_GRADIENT == 1
            [FUNk_x  FUNk_y  FUNk_z]= dSinCosExpFun3D_2param(DATALOC.mu,xNEW(:,1),xNEW(:,2),xNEW(:,3)) ;
            
            
        end
        
    otherwise        
        error('Option not implemented')       
        
end

%dbstop('82')
DATALOC = DefaultField(DATALOC,'ImposeVolumeConstraint_directly',0);
if DATALOC.ImposeVolumeConstraint_directly == 1
    if EVALUATE_FUN == 1
        newcol = ones(size(FUNk,1),1) ;
        FUNk = [newcol FUNk] ;
    end
    if EVALUATE_GRADIENT == 1
        newcol = zeros(size(FUNk_x,1),1) ;
        
        FUNk_x = [newcol FUNk_x] ;
        FUNk_y = [newcol FUNk_y] ;
        if size(xNEW,2) ==3
            FUNk_z = [newcol FUNk_z] ;
        end
    end
end

if EVALUATE_FUN == 1
    if ~isempty(VSinv)
    PHIk_y =FUNk*VSinv ;
    else
        PHIk_y = FUNk ; 
    end
end

if EVALUATE_GRADIENT == 1
    if ~isempty(VSinv)
        if size(xNEW,2) ==2
            dPHIk_y{1}  =FUNk_x*VSinv ;
            dPHIk_y{2}  =FUNk_y*VSinv ;
        elseif size(xNEW,2) ==3
            dPHIk_y{1}  =FUNk_x*VSinv ;
            dPHIk_y{2}  =FUNk_y*VSinv ;
            dPHIk_y{3}  =FUNk_z*VSinv ;
        end        
    else
         if size(xNEW,2) ==2
            dPHIk_y{1}  =FUNk_x ;
            dPHIk_y{2}  =FUNk_y ;
        elseif size(xNEW,2) ==3
            dPHIk_y{1}  =FUNk_x ;
            dPHIk_y{2}  =FUNk_y ;
            dPHIk_y{3}  =FUNk_z ;
        end
        
    end
end


