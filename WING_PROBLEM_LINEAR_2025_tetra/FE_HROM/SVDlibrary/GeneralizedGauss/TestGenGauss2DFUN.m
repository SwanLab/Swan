function  [xGAUSS,wGAUSS,Xf,W,DATAapprox]  =  TestGenGauss2DFUN(NAMEDATA,M,P,DATA)
% Testing the Empirical Cubature Method
% Joaquín A. Hernández, June 8-th - 2016 (jhortega@cimne.upc.edu)
% UPDATED.  APRIL-4-2020 (covid 19 quarantine)
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

eval(NAMEDATA) ;

% Generating snapshot matrix and weight vector (step 2 in Box 5.1)
[Xf,W,mu,x,MPOINTS,XY] = GenerateFunExampleNDIM2d(M,P,DATA,DATAFUN);
DATA = DefaultField(DATA,'PartitionedMethodECM',0 ) ;
DATA = DefaultField(DATA,'SVD_BEFORE_zeroaverage',0 ) ;
DATA = DefaultField(DATA,'MULTIPLY_BY_INTEGRAL',0 ) ;
DATA = DefaultField(DATA,'AN2009method',0 ) ;
DATA = DefaultField(DATA,'EMPLOY_XIAO_FORSELECTINGPOINTS',0 ) ;



DATA.x = x ;
% if DATA.EMPLOY_XIAO_FORSELECTINGPOINTS==1
%     Xf = bsxfun(@times,Xf,sqrt(W)) ;
%     if DATA.RANDOMIZED_SVD == 0
%         [PHI,S,V,ERRORsvd,TIMETOTAL] = SVD_BLOCKWISE(Xf,DATA);
%     else
%         [PHI,S,V,ERRORsvd,TIMETOTAL,DATA] = SVD_BLOCKWISE_random(Xf,DATA);
%     end
%     DATAOUT.SingV_Jmin = S ; 
%     DATAOUT.V_Jmin = V ; 
%     z = DEIMfun(PHI) ; 
%     b  = PHI'*sqrt(W) ;
%     alpha = PHI(z,:)\b ; 
%     w  = alpha.*sqrt(W(z)) ; 
%     DATAOUT.J = PHI' ; 
% else
%     if DATA.AN2009method == 1
%         [z,w,errorGLO,DATAOUT] = OptimizedCubature(Xf,W,DATA);
%     else
%         if DATA.PartitionedMethodECM ==1
%             [z,w,errorGLO] = EmpiricalCubatureMethodPart(Xf,W,DATA);
%         elseif  DATA.PartitionedMethodECM ==0
            [z,w,errorGLO,DATAOUT] = ECM1dom(Xf,W,DATA);
%         elseif DATA.PartitionedMethodECM ==2
%             % Hybrid approach
%             % ---------------
%             [z,w,errorGLO,DATAOUT] = EmpiricalCubatureMethodHybr(Xf,W,DATA);
%         end
%     end
% end

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
%DATA = DefaultField(DATA,'AN2009method',0) ;
% if DATA.AN2009method ==1  | DATA.IncludeWeights_Jmin==0
%     DATALOC.PHIisTRANSPOSE = 1 ; 
% else
      DATALOC.PHIisTRANSPOSE = 0 ; 
%end
DATALOC.xLIM =DATAFUN.xLIM ; 
DATALOC.INITIAL_DISCRETIZATION_TENSORPRODUCT =  DATA.INITIAL_DISCRETIZATION_TENSORPRODUCT ; 
DATALOC.APPROX_FUN__DERI =DATA.APPROX_FUN__DERI ;  % 2-apr-2020    
DATALOC.APPROX_FUN__DERI.COOR = x ; 
%dbstop('76')
DATALOC.SAVE_XGAUSS = DATA.SAVE_XGAUSS ;
 [xGAUSS,wGAUSS,DATAapprox] = GenGauss2DapproxFUN(DATAOUT.J,...
    W,x,z,w,DATALOC,MPOINTS,DATAOUT.SingV_Jmin,DATAOUT.V_Jmin,XY) ; 
TIMEela = toc(TIMEela);
disp(['TIME gaussian =',num2str(TIMEela)])




