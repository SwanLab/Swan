function  [xGAUSS,wGAUSS,Xf,W,DATAapprox]  =  TestingGeneralizedGauss2DFUN(NAMEDATA,M,P,DATA)
% Testing the Empirical Cubature Method
% Joaquín A. Hernández, June 8-th - 2016 (jhortega@cimne.upc.edu)
% please cite paper: Dimensional hyperreduction of nonlinear parameterized finite element models via
% empirical cubaturre
%  hernandez2016dimensional.pdf
%----------------------------------------------------------------------------------------------------
% TEsting function: page 12 of "A ‘best points’ interpolation method for efficient approximation
% of parametrized functions", by N.C. Nguyen et al, 2007
%---------------------------------------------------
format long g
% -----------------------------------------------------------
% INPUTS
% ---------------------------

%NAMEDATA ='DATA_GG2D_poly' ;'DATA_GG2D_HOLE';  'DATA_GenGauss2D';'DATA_GenGauss' ; 'DATA_orthopol';'DATA_piece';   'DATA_recursiveSEVD' ; % 'DATA_2'; %  'DATA_complex2'% ; 'DATA_1'

% -------------------------------------------------
% END INPUTS
% -------------------------------------------------
DATA = DefaultField(DATA,'READ_AGAIN_INPUTS',1) ;
if DATA.READ_AGAIN_INPUTS == 1
    eval(NAMEDATA) ;
else
    DATAFUN = DATA.DATAFUN ;
end

% Generating snapshot matrix and weight vector (step 2 in Box 5.1)
[Xf,W,mu,x,MPOINTS,XY] = GenerateFunExampleNDIM2d(M,P,DATA,DATAFUN);
DATA = DefaultField(DATA,'PartitionedMethodECM',0 ) ;
DATA = DefaultField(DATA,'SVD_BEFORE_zeroaverage',0 ) ;
DATA = DefaultField(DATA,'MULTIPLY_BY_INTEGRAL',0 ) ;
DATA = DefaultField(DATA,'AN2009method',0 ) ;
DATA = DefaultField(DATA,'EMPLOY_XIAO_FORSELECTINGPOINTS',0 ) ;


DATA = DefaultField(DATA,'PERCENTAGE_POINTS_EXCLUDED_BOUNDARIES_DECM',[]) ;
DATA.ECM_POINTS_INCLUDE =[] ;
if ~isempty(DATA.PERCENTAGE_POINTS_EXCLUDED_BOUNDARIES_DECM)
    ndim = size(x,2) ;
    %   nperc_direct = ((DATA.PERCENTAGE_POINTS_EXCLUDED_BOUNDARIES_DECM)^(1/ndim));
    nperc_include = 100-DATA.PERCENTAGE_POINTS_EXCLUDED_BOUNDARIES_DECM ;
    xLIM = DATAFUN.xLIM ;
    xCENT =sum(xLIM,2)/2 ;
    xSPAN =  (xLIM(:,2)-xLIM(:,1))/2;
    
    xSPANincludeHALF = xSPAN*nperc_include/100;
    
    xINCLUDE_min = xCENT -xSPANincludeHALF ;
    xINCLUDE_max =  xCENT +xSPANincludeHALF ;
    xINC = [xINCLUDE_min,xINCLUDE_max] ;
    % Define polygon
    POLI = [xINC(1,1),xINC(2,1)
        xINC(1,2),xINC(2,1)
        xINC(1,2),xINC(2,2)
        xINC(1,1),xINC(2,2)
        xINC(1,1),xINC(2,1)] ;
    
    
    IINN =      inpolygon(x(:,1),x(:,2),POLI(:,1),POLI(:,2)) ;
    
    IndInclude = find(IINN~=0) ;
    DATA.ECM_POINTS_INCLUDE= IndInclude ;
    
    
end



DATA.x = x ;
if DATA.EMPLOY_XIAO_FORSELECTINGPOINTS==1
    Xf = bsxfun(@times,Xf,sqrt(W)) ;
    if DATA.RANDOMIZED_SVD == 0
        [PHI,S,V,ERRORsvd,TIMETOTAL] = SVD_BLOCKWISE(Xf,DATA);
    else
        [PHI,S,V,ERRORsvd,TIMETOTAL,DATA] = SVD_BLOCKWISE_random(Xf,DATA);
    end
    DATAOUT.SingV_Jmin = S ;
    DATAOUT.V_Jmin = V ;
    z = DEIMfun(PHI) ;
    b  = PHI'*sqrt(W) ;
    alpha = PHI(z,:)\b ;
    w  = alpha.*sqrt(W(z)) ;
    DATAOUT.J = PHI' ;
else
    if DATA.AN2009method == 1
        [z,w,errorGLO,DATAOUT] = OptimizedCubature(Xf,W,DATA);
    else
        if DATA.PartitionedMethodECM ==1
            [z,w,errorGLO] = EmpiricalCubatureMethodPart(Xf,W,DATA);
        elseif  DATA.PartitionedMethodECM ==0
            [z,w,errorGLO,DATAOUT] = EmpiricalCubatureMethod1dom(Xf,W,DATA);
        elseif DATA.PartitionedMethodECM ==2
            % Hybrid approach
            % ---------------
            [z,w,errorGLO,DATAOUT] = EmpiricalCubatureMethodHybr(Xf,W,DATA);
        end
    end
end

disp('---------------------------------------------')
%disp('TOTAL ERROR (SVD of Xf and integration)')
ExactIntegral = Xf'*W ;
ApproximatedInt = Xf(z,:)'*w ;


%clear Xf

disp('ERROR (%)')
EE = norm(ExactIntegral-ApproximatedInt)/norm(ExactIntegral)*100;
disp([num2str(EE)])
%EEE = norm(Exactintegra)
%%%%%%%%%%%% -ADDITIONAL STEP: NNLSQ WITH THE SELECTED SET OF POINTS





%%%% SECOND STEP (solving moment matching equations )
% ----------------------------------------------------
DATALOC = DATA.DATAREMOVEPOINTS  ;
TIMEela = tic;
DATA = DefaultField(DATA,'AN2009method',0) ;
if DATA.AN2009method ==1  | DATA.IncludeWeights_Jmin==0
    DATALOC.PHIisTRANSPOSE = 1 ;
else
    DATALOC.PHIisTRANSPOSE = 0 ;
end
DATALOC.xLIM =DATAFUN.xLIM ;
DATALOC.INITIAL_DISCRETIZATION_TENSORPRODUCT =  DATA.INITIAL_DISCRETIZATION_TENSORPRODUCT ;
DATALOC.APPROX_FUN__DERI =DATA.APPROX_FUN__DERI ;  % 2-apr-2020
DATALOC.APPROX_FUN__DERI.COOR = x ;
%dbstop('76')
DATALOC.SAVE_XGAUSS = DATA.SAVE_XGAUSS ;
DATA = DefaultField(DATA,'USE_SVD_SOLVE_EQUATIONS',0) ;


DATALOC.USE_SVD_SOLVE_EQUATIONS=  DATA.USE_SVD_SOLVE_EQUATIONS ;
DATALOC.METHOD_SOLVE_UNDERDETERMINED_SYSTEM_EQUATIONS=  DATA.METHOD_SOLVE_UNDERDETERMINED_SYSTEM_EQUATIONS ;
DATA = DefaultField(DATA,'USE_VERSION_REMOVE_POINTS_AT_END_FROM_FE_IMPLEMENTATION',0) ; % = 1;
DATALOC.IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA=  DATA.IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA ;

%if DATA.USE_VERSION_REMOVE_POINTS_AT_END_FROM_FE_IMPLEMENTATION == 0
[xGAUSS,wGAUSS,DATAapprox] = GeneralizedGauss2D(DATAOUT.J,...
    W,x,z,w,DATALOC,MPOINTS,DATAOUT.SingV_Jmin,DATAOUT.V_Jmin,XY) ;
% else
%     DATALOC.ExactIntegral = ExactIntegral ;
%
%     % Implemented 15-Oct-2021. More in line with what is written in
%     %   /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/OptimalECM
%     VAR_SMOOTH_FE = [] ;
%     [xGAUSS,wGAUSS,DATALOC,ELEMENTS_xGAUSS,HISTORY,VAR_SMOOTH_FE] ...
%     = GenGauss2D3Dfebased(x,w,DATALOC,x(z,:),VAR_SMOOTH_FE)  ;
% end


TIMEela = toc(TIMEela);
disp(['TIME gaussian =',num2str(TIMEela)])




