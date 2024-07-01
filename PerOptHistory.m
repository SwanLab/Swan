
clear;
close all;

fig = open("IsoRelPer.fig");

ch       = fig.Children;
itersIso = ch(end).Children.XData(2:end);
PIso     = ch(end).Children.YData(2:end);
close all

fig = open("AniRelPerA0B38.fig");

ch        = fig.Children;
itersAni1 = ch(end).Children.XData(2:101);
PAni1     = ch(end).Children.YData(2:101);
close all

fig = open("AniRelPerA0B1940.fig");

ch        = fig.Children;
itersAni2 = ch(end).Children.XData(2:101);
PAni2     = ch(end).Children.YData(2:101);
close all

fig = open("AniRelPerA8B1940.fig");

ch        = fig.Children;
itersAni3 = ch(end).Children.XData(2:101);
PAni3     = ch(end).Children.YData(2:101);
close all

figure
plot(itersIso,PIso,itersAni1,PAni1,itersAni2,PAni2,itersAni3,PAni3)
grid on
grid minor
xlim([-5 105])
xlabel('Iteration', 'Interpreter','latex')
ylabel('$Per_\epsilon(\chi)$','Interpreter','latex')
legend('Isotropic relative perimeter','Anisotropic relative perimeter with $\alpha=0$ and $\beta=3\pi/8$',...
    'Anisotropic relative perimeter with $\alpha=0$ and $\beta=19\pi/40$',...
    'Anisotropic relative perimeter with $\alpha=\pi/8$ and $\beta=19\pi/40$','interpreter','latex')