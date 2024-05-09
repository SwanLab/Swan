% Descifrando etaMax levelset

clear;
close all;

%fig = open("lsMagicRatio.fig");
fig = open("lsMagicRatio2.fig");

ch         = fig.Children;
iters      = ch(1).Children.XData(27:end);
beta       = ch(1).Children.YData(27:end);
alpha      = ch(2).Children.YData(27:end);
theta      = ch(3).Children.YData(27:end);
magicRatio = ch(4).Children.YData(27:end);
eta        = ch(8).Children.YData(27:end);
k          = ch(10).Children.YData(27:end);
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

b=0;
a=0.87;
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