function ProposingAnalyticSuperEllipseExponent
rhoMin = 0.01;
rhoMax = 0.99;
rhoDif = (rhoMax-rhoMin);
q0 = 32;
q1 = 25;
beta = 300;
qRhoMax = @(rho) q0 + ((rho - rhoMin)/(rhoDif)).^(beta)*(q1-q0);

rhoS = rhoMin + 0.2*(rhoMax-rhoMin);
A = [rhoMax*(rhoMax - 2*rhoS) -rhoMin*(rhoMin - 2*rhoS);...
     2*rhoS  -2*rhoS;...
     -1 1];


vect =  1/((rhoMax - rhoMin)*((rhoMax - rhoS) - (rhoS - rhoMin)))*A*[q0;q1];
qVect = @(rho) [ones(size(rho,1),1) rho rho.^2];
qRhoMin = @(rho) sum(bsxfun(@times,qVect(rho)',vect),1);


rhoV(:,1) = linspace(rhoMin,rhoMax,10000);
figure(1)
plot(rhoV,[qRhoMax(rhoV),qRhoMin(rhoV)'])

end