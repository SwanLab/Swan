% Descifrando etaMax levelset

clear;
close all;

fig = open("lsMagicRatio5SuperSmoothLarge.fig");
% fig = open("lsMagicRatio3.fig");

ch         = fig.Children;
iters      = ch(1).Children.XData(27:end); % ex 2: 27-end; ex 3.1: 3-end; ex3.2: 3-20
beta       = ch(1).Children.YData(27:end);
alpha      = ch(2).Children.YData(27:end);
theta      = ch(3).Children.YData(27:end);
magicRatio = abs(ch(4).Children.YData(27:end));
eta        = ch(8).Children.YData(27:end);
k          = ch(10).Children.YData(27:end);
gV         = abs(ch(13).Children.YData(27:end));
close all;

etaMax = 1./magicRatio;

combo=beta./eta;
figure
semilogy(iters,magicRatio,iters,theta,iters,k,iters,alpha,iters,beta,iters,combo,iters,gV)
grid on
grid minor
legend('Predicted tau','theta','kappa','alpha','beta','Proposed combination','volume constraint')