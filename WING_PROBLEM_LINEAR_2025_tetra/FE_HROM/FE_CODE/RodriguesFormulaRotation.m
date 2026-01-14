function R = RodriguesFormulaRotation(theta,n)
% Matrix representing a rotation theta around an axis n 
% https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
if nargin == 0
    n = [0 0 1]' ;
    theta = pi/2; 
end

K = [0 -n(3) n(2)
     n(3) 0  -n(1)
    -n(2) n(1) 0] ; 
R = eye(3) +K+ (K*K)/(1+cos(theta)) ; 