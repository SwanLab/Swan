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