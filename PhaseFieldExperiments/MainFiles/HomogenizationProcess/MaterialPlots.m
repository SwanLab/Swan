clc,clear,close all
load('HorizontalCrack.mat')

figure(1)
tiledlayout(3,3)
nexttile
plot(phi,squeeze(mat(1,1,1,1,:)),'.')
title('C11')
nexttile
plot(phi,squeeze(mat(1,1,2,2,:)),'.')
title('C12')
nexttile
plot(phi,squeeze(mat(1,1,1,2,:)),'.')
title('C13')
nexttile
plot(phi,squeeze(mat(2,2,1,1,:)),'.')
title('C21')
nexttile
plot(phi,squeeze(mat(2,2,2,2,:)),'.')
title('C22')
nexttile
plot(phi,squeeze(mat(2,2,1,2,:)),'.')
title('C23')
nexttile
plot(phi,squeeze(mat(1,2,1,1,:)),'.')
title('C31')
nexttile
plot(phi,squeeze(mat(1,2,2,2,:)),'.')
title('C32')
nexttile
plot(phi,squeeze(mat(1,2,1,2,:)),'.')
title('C33')

figure(2)
tiledlayout(3,3)
nexttile
hold on
plot(phi,squeeze(mat(1,1,1,1,:)),'.')
fplot(degradation.fun(1,1,1,1,:),[0,1])
title('C11')
nexttile
hold on
plot(phi,squeeze(mat(1,1,2,2,:)),'.')
fplot(degradation.fun(1,1,2,2,:),[0,1])
title('C12')
nexttile
hold on
plot(phi,squeeze(mat(1,1,1,2,:)),'.')
fplot(degradation.fun(1,1,1,2,:),[0,1])
title('C13')
nexttile
hold on
plot(phi,squeeze(mat(2,2,1,1,:)),'.')
fplot(degradation.fun(2,2,1,1,:),[0,1])
title('C21')
nexttile
hold on
plot(phi,squeeze(mat(2,2,2,2,:)),'.')
fplot(degradation.fun(2,2,2,2,:),[0,1])
title('C22')
nexttile
hold on
plot(phi,squeeze(mat(2,2,1,2,:)),'.')
fplot(degradation.fun(2,2,1,2,:),[0,1])
title('C23')
nexttile
hold on
plot(phi,squeeze(mat(1,2,1,1,:)),'.')
fplot(degradation.fun(1,2,1,1,:),[0,1])
title('C31')
nexttile
hold on
plot(phi,squeeze(mat(1,2,2,2,:)),'.')
fplot(degradation.fun(1,2,2,2,:),[0,1])
title('C32')
nexttile
hold on
plot(phi,squeeze(mat(1,2,1,2,:)),'.')
fplot(degradation.fun(1,2,1,2,:),[0,1])
title('C33')

figure(2)
tiledlayout(3,3)
nexttile
fplot(degradation.fun(1,1,1,1,:),[0,1])
title('C11')
nexttile
fplot(degradation.fun(1,1,2,2,:),[0,1])
title('C12')
nexttile
fplot(degradation.fun(1,1,1,2,:),[0,1])
title('C13')
nexttile
fplot(degradation.fun(2,2,1,1,:),[0,1])
title('C21')
nexttile
fplot(degradation.fun(2,2,2,2,:),[0,1])
title('C22')
nexttile
fplot(degradation.fun(2,2,1,2,:),[0,1])
title('C23')
nexttile
fplot(degradation.fun(1,2,1,1,:),[0,1])
title('C31')
nexttile
fplot(degradation.fun(1,2,2,2,:),[0,1])
title('C32')
nexttile
fplot(degradation.fun(1,2,1,2,:),[0,1])
title('C33')










