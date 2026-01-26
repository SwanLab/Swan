function   DATA_regress_weightsECM = MAW_ECM_regression(qLATENT,wALL,DATA_interp ) ; 

 %--------------------------------------------------------------------------
 % JAHO, 13-Sept-2025, Thriller CAFFE BAR, Sarajevo, Bosnia. 
%--------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat') 
end
      


DATA_interp = DefaultField(DATA_interp,'LocalToleranceSVD_ECMpoints',1e-5);

DATA_regress_weightsECM= BsplinesLeastSquares_MAW_ECM(DATA_interp, qLATENT, wALL);


 