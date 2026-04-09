clc,clear,close all
load('HorizontalCrack.mat')

figure()
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