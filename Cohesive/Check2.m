load('damageHistory.mat') 
close all

lambda = lambda(:,1);
d = d(:,1);

% Experimental
scatter(lambda,d,'.')
hold on


% Theoretical
jumpCrit  = 1.25e-7; 
jumpFinal = 0.025e-3;
isJump = lambda > 1e-12; isJump = 1;
dTheoretical = min(isJump .* (jumpFinal .* (lambda - jumpCrit)./ (lambda .* (jumpFinal - jumpCrit))),1);
dTheoretical = max(dTheoretical,0);



plot(lambda,dTheoretical);

