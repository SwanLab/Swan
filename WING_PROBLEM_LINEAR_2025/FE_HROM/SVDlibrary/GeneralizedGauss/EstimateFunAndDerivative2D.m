function [f_approx,d_approx,INDfit,POLYCOEFF2D] = EstimateFunAndDerivative2D(COOR,F,xnew,DATA)
% 2D-Data fitting.

% INPUTS
%
% 1.COOR: Coordinates
% 2.xnew: Point at which the function and its derivative  are to be
% estimated
% 3. F: Value of the function at each point COOOR
% 4. DATA.SELECTION_METHOD = 'RECTANGULAR',  'CIRCLE'
% 5. DATA.ORDER_polynomial =[m,n]   (x^m*y^n)
% 6. DATA.N_POINTS_FITTING = [npoints_x,npoints_y]  ; % Points used for the
% interpolation
% 7. DATA.POLYCOEFF2D.  (Coefficients Polynom., empty by default)
%
% OUTPUTs
% f_approx -->Approximated function
% d_approx --> Derivative
% POLYCOEFF2D --> Coefficients Polynom.
% Written by JAHO, 2th April 2020, 20th day QUARANTINE, COVID-19
% ---------------------------------------------------------------

DATA = DefaultField(DATA,'ORDER_polynomial',[3,3]) ;
NMON = (DATA.ORDER_polynomial(1)+1)^2;
COEFF_over = 3;

NPOINTS_FITTING_x = ceil(sqrt(NMON*COEFF_over)) ;
DATA = DefaultField(DATA,'N_POINTS_FITTING',NPOINTS_FITTING_x*[1,1]);
DATA = DefaultField(DATA,'SELECTION_METHOD','CIRCLE');

DATA =  DefaultField(DATA,'POLYCOEFF2D',[]);
DATA =  DefaultField(DATA,'INDfit',[]);


if isempty(DATA.POLYCOEFF2D)
    if isempty(DATA.INDfit)
        [Indcls,dummy] = knnsearch(COOR,xnew) ;
        xCLOS1 = COOR(Indcls,:) ;
        % Which is the nearest point to xCLOS
        [Indcls2,dummy] = knnsearch(COOR,xCLOS1,'k',2) ;
        Indcls2 = Indcls2(2) ;
        xCLOS2 = COOR(Indcls2,:) ;
        h = norm(xCLOS2-xCLOS1) ;  % Typical distance
        
        switch DATA.SELECTION_METHOD
            case 'RECTANGULAR'
                %  error('option not ready yet')
                for idim = 1:ndim
                    DISTp = (abs(COOR(:,idim) - xnew(idim)));
                    radius = h*DATA.N_POINTS_FITTING(idim)/2 ;
                    [IND{idim}] =find(DISTp <=radius )  ;
                end
                INDfit = intersect(IND{1},IND{2}) ;
                
            case 'CIRCLE'
                INDfit = knnsearch(COOR,xnew,'k',prod(DATA.N_POINTS_FITTING)) ;
        end
        
        
    else
        
        INDfit = DATA.INDfit ;
        
    end
    COORfit = COOR(INDfit,:) ;
    
    %  plot3(COORfit(:,1),COORfit(:,2),zeros(length(INDfit),1),'ro');
    % ----------------------------------------
    % E) POLYNOMIAL FITTING
    % ----------------------------------------
    % /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/SVDlibrary/GeneralizedGauss/POLY2D
    %
    Ffit = F(INDfit) ;
    POLYCOEFF2D = polyFit2D(Ffit,COORfit(:,1),COORfit(:,2),DATA.ORDER_polynomial(1),DATA.ORDER_polynomial(2)) ;
else
    POLYCOEFF2D = DATA.POLYCOEFF2D ;
end

%
f_approx =  polyVal2D(POLYCOEFF2D,xnew(1),xnew(2),DATA.ORDER_polynomial(1),DATA.ORDER_polynomial(2)) ;
%
%
d_approx = zeros(1,2);

[d_approx(1),d_approx(2)] =  polyDer2D(POLYCOEFF2D,xnew(1),xnew(2),DATA.ORDER_polynomial(1),DATA.ORDER_polynomial(2)) ;

