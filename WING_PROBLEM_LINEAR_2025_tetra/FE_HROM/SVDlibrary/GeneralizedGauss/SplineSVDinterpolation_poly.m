function [Uypoly,Vxpoly,Uypoly_der,Vxpoly_der,S,nERRORx,nERRORy] = SplineSVDinterpolation_poly(x,y,Fmat,DATA)
% This function apply the SVD to Fmat, and interpolate, via splines or polynomial fitting, the
% left and right singular values. SMOOTHING can be applied to such curves
% to eliminate noise
% JAHO, 3-April-2020, 21th -Quarantine, COVID-19
% -----------------------------------------------------------------------------------
if nargin == 0
    load('tmp2.mat')
    DATA.TOL_SVD = 1e-2; 
    DATA.PLOT_APPROXIMATED = 1; 
        DATA.ORDER_GLOBAL = 8; 
        DATA.METHOD = 'SPLINE' ; 
        DATA.SMOOTH_ACTIVE  = 1; 
end

DATA = DefaultField(DATA,'TOL_SVD',1e-4) ;
DATA = DefaultField(DATA,'SMOOTH_ACTIVE',0) ;
DATA = DefaultField(DATA,'METHOD','SPLINE') ;

DATA = DefaultField(DATA,'ORDER_GLOBAL',7) ;
DATA = DefaultField(DATA,'PLOT_FUNCTIONS',0) ;
DATA = DefaultField(DATA,'NFUN',1) ;
DATA = DefaultField(DATA,'PLOT_APPROXIMATED',1) ;

% DATA.SVD_INTERPOLATION.METHOD = 'POLYNOMIAL' ; 'SPLINE' ;
% DATA.SVD_INTERPOLATION.ORDER_GLOBAL = 7;

% ---------------------------------------
TOL = DATA.TOL_SVD;
DATASVD.RELATIVE_SVD = 1;
[Uy,S,Vx] = RSVDT(Fmat,TOL,[],0,DATASVD) ;






Uypoly = cell(length(S),1) ;
Vxpoly = cell(length(S),1) ;
ERROR_x = zeros(length(S),1) ; 
ERROR_y = zeros(length(S),1) ; 


Uypoly_der = cell(length(S),1) ;
Vxpoly_der = cell(length(S),1) ;
for i = 1:length(S)
    
    if  DATA.SMOOTH_ACTIVE == 1
        Vloc  = smooth(Vx(:,i)) ;
        Uloc = smooth(Uy(:,i)) ;
    else
        Vloc = Vx(:,i) ;
        Uloc =  (Uy(:,i)) ;
    end
    
    
    
    
    switch  DATA.METHOD
        case 'SPLINE'
            Vxpoly{i} =  spline(x,Vloc) ;
            Uypoly{i} = spline(y,Uloc) ;
            Vxpoly_der{i} = ppDer(Vxpoly{i}) ;
            Uypoly_der{i} = ppDer(  Uypoly{i});
            FUN_EVAL ='ppval' ;
            
        case   'POLYNOMIAL'
            % SALIR = 0 ;
            ORDER = DATA.ORDER_GLOBAL ;
            
            [  Vxpoly{i},EEEx]=  polyfit(x(:),Vloc,ORDER) ;
            [ Uypoly{i},EEy] = polyfit(y(:),Uloc,ORDER) ;
            
            ERROR_x(i) = EEEx.normr ; 
               ERROR_y(i) = EEy.normr ; 
            [ Vxpoly_der{i}]= polyder(Vxpoly{i}) ;
            [  Uypoly_der{i} ]= polyder(  Uypoly{i});
            FUN_EVAL ='polyval' ;
            
            
        otherwise
            error('Option not implemented')
    end
    
end


if all(ERROR_x ==0)
    ERROR_x = []; 
    ERROR_y = [] ; 
    LGx ='';
    LGy = ''; 
    nERRORx =[] ; 
    nERRORy = [] ; 
else
    nERRORx =  norm(ERROR_x) ; 
    nERRORy =  norm(ERROR_y) ; 
    LGx = ['; ERROR intp. =  ',num2str(nERRORx)]; 
    LGy = ['; ERROR intp. =  ',num2str(nERRORy)]; 
end


if DATA.PLOT_FUNCTIONS == 1
    figure(100+DATA.NFUN)
    1e-3
    subplot(2,1,1)
     
    hold on
    h = [] ;
    xlabel('x')
    ylabel('V(x)')
    title(['Right singular vectors, ', 'FUN =',num2str(DATA.NFUN),LGx])
    LLL = {} ;
    for i = 1:length(S)
        h(end+1) = plot(x,Vx(:,i),'color',rand(1,3)) ;
        LLL{end+1} = ['V_',num2str(i),'; s = ',num2str(S(i))] ;
        
        if DATA.PLOT_APPROXIMATED == 1
            % Approximated
            Vx_approx = feval(FUN_EVAL,Vxpoly{i},x) ;
            h(end+1) = plot(x,Vx_approx,'color',rand(1,3)) ;
            LLL{end+1} = ['V_',num2str(i),' (approx)'  ] ;
        end
        % ------------------------------
        
        
    end
    legend(h,LLL)
    
  subplot(2,1,2) 
    hold on
    h = [] ;
    xlabel('y')
    ylabel('U(y)')
    title(['Left singular vectors',LGy])
    LLL = {} ;
    for i = 1:length(S)
        h(end+1) = plot(y,Uy(:,i),'color',rand(1,3)) ;
        LLL{end+1} = ['U_',num2str(i),'; s = ',num2str(S(i))] ;
        
         if DATA.PLOT_APPROXIMATED == 1
        % Approximated
        Uy_approx = feval(FUN_EVAL,Uypoly{i},y) ;
        h(end+1) = plot(y,Uy_approx,'color',rand(1,3)) ;
        LLL{end+1} = ['U_',num2str(i),' (approx)'  ] ;
         end
        % ------------------------------
    end
    legend(h,LLL)
    
    
    
    
    
end

