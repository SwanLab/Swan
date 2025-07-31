clc,clear,close all
load('CircleAreaDerivative2.mat')

figure()
tiledlayout(1,3)
nexttile
fplot(degradation.fun{1,1,1,1},[0 1])
hold on
plot(phi,squeeze(mat(1,1,1,1,:)),'.')
title('C11')

nexttile
fplot(degradation.fun{1,1,2,2},[0 1])
hold on
plot(phi,squeeze(mat(1,1,2,2,:)),'.')
title('C12')

nexttile
fplot(degradation.fun{1,2,1,2},[0 1])
hold on
plot(phi,squeeze(mat(1,2,1,2,:)),'.')
title('C33')

figure()
tiledlayout(1,3)
nexttile
fplot(degradation.dfun{1,1,1,1},[0 1])

title('dC11')

nexttile
fplot(degradation.dfun{1,1,2,2},[0 1])
title('dC12')

nexttile
fplot(degradation.dfun{1,2,1,2},[0 1])
title('dC33')

figure()
tiledlayout(1,3)
nexttile
fplot(degradation.ddfun{1,1,1,1},[0 1])

title('ddC11')

nexttile
fplot(degradation.ddfun{1,1,2,2},[0 1])
title('ddC12')

nexttile
fplot(degradation.ddfun{1,2,1,2},[0 1])
title('ddC33')