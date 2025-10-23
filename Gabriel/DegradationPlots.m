load('HomogenizationResults.mat');
simp = @(x) x.^3;

tiledlayout(1,3)
nexttile
hold on
plot(volFrac,squeeze(Chomog(1,1,1,1,:)),'LineStyle','none','Marker','o')
fplot(Interpolation.fun(1,1,1,1),[0 1])
fplot(simp,[0 1])

nexttile
hold on
plot(volFrac,squeeze(Chomog(1,1,2,2,:)),'LineStyle','none','Marker','o')
fplot(Interpolation.fun(1,1,2,2),[0 1])
fplot(simp,[0 1])

nexttile
hold on
plot(volFrac,squeeze(Chomog(1,2,1,2,:)),'LineStyle','none','Marker','o')
fplot(Interpolation.fun(1,2,1,2),[0 1])
fplot(simp,[0 1])