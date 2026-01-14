function [tau,tauDER1,tauDER2] = tauFUN_identity(qL,DUMMY)
if nargin == 0
    load('tmp1.mat')
end
 
tau = qL ; 
tauDER1 = speye(size(qL,1)) ; 
tauDER2 =[] ; 
