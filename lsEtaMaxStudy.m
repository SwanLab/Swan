% Descifrando etaMax levelset

clear;
close all;

fig = open("lsMagicRatio.fig");

ch         = fig.Children;
iters      = 2:24;
beta       = ch(1).Children.YData(3:end);
alpha      = ch(2).Children.YData(3:end);
theta      = ch(3).Children.YData(3:end);
magicRatio = ch(4).Children.YData(3:end);
eta        = ch(8).Children.YData(3:end);
k          = ch(10).Children.YData(3:end);
close all;

etaMax = 1./magicRatio;

bVec = -20:0.1:20;
for j=1:length(bVec)
[aVec(j),error(j)]  = comptueKappaExponentFromTheta(etaMax,k,theta,bVec(j));
end

figure
plot(bVec,aVec)
grid on
grid minor
xlabel('Theta sinus exponent')
ylabel('Kappa exponent')

figure
plot(bVec,error)
grid on
grid minor
xlabel('Theta sinus exponent')
ylabel('Error')

b=10;
a=0.5;
combo = (k.^a).*(sin(theta)).^b;

figure
semilogy(iters,etaMax,iters,theta,iters,k,iters,alpha,iters,beta,iters,combo)
grid on
grid minor
legend('Real etaMax','theta','kappa','alpha','beta','Proposed combination')

function [aRes,error]=comptueKappaExponentFromTheta(x,k,theta,b)
fun   = @(a) norm((k.^a).*(sin(theta)).^b-x);
aRes  = fminbnd(fun,-20,20);
error = fun(aRes);
end