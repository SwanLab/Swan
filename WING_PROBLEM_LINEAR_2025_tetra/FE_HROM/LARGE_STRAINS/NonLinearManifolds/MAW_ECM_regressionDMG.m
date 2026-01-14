function   DATA_regress_weightsECM = MAW_ECM_regressionDMG(qLATENT,wALL,DATA_interp,qLATENTall ) ; 
 %--------------------------------------------------------------------------
 % JAHO,  28-Oct-2025, Tuesday HGs Pedralbes, Barcleona
%--------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat') 
end
      


DATA_interp = DefaultField(DATA_interp,'LocalToleranceSVD_ECMpoints',1e-5);

DATA_regress_weightsECM= BsplinesLeastSquares_MAW_ECMdmg(DATA_interp, qLATENT,qLATENTall(1,:), qLATENTall(2,:), wALL);


 