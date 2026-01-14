function C = DamageMatStressStress_2(S,nstrain)
% Compute matrix C = S*S^T, damage model, vectorized fashion
% JAHO, 22-Oct-2025, Balmes 185, Barcelona
if nargin == 0
    load('tmp2.mat')
    S = [S;S] ; 
end 
n = size(S,1) ; 
C = zeros(n,nstrain); 
for istrain = 1:nstrain
    istrainGLO = istrain:nstrain:n ; 
    for jstrain = 1:nstrain 
        jstrainGLO = jstrain:nstrain:n; 
        C(istrainGLO,jstrainGLO)  = S(istrainGLO)*S(jstrainGLO)' ; 
    end
end