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

figure(3)
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

%% Comparative with Rao et al. (2022)
clc,clear,close all
load('HorizontalCrack.mat')
% prescription = slmset('plot','on','decreasing','on','interiorknots','free','knots',20,...
%                        'leftvalue',2*C22,'rightvalue',0,'concavedown','on');
% 
% prescription = slmset('plot','on','decreasing','on','knots',0:0.025:1,...
%                        'leftvalue',2*C22,'rightvalue',0);

C = cell(2,2,2,2);
Cmat = zeros(2,2,2,2,length(phi));
CmatCOV = zeros(2,2,2,2,length(phi));
dCmat = zeros(2,2,2,2,length(phi));
d2Cmat = zeros(2,2,2,2,length(phi));

idx = 2;
for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                prescription = slmset('decreasing','on','knots',0:0.020:1,...
                    'leftvalue',mat(i,j,k,l,1));
                if i~=1 || j~=1 || k~=1 || l~=1
                    C{i,j,k,l} = slmengine(phi.^2,squeeze(mat(i,j,k,l,:)),prescription,'rightvalue',0);
                else
                    C{i,j,k,l} = slmengine(phi.^2,squeeze(mat(i,j,k,l,:)),prescription);
                end
                A = slmeval(phi,C{i,j,k,l});
                Cmat(i,j,k,l,:) = A;
                CmatCOV(i,j,k,l,:) = A.^(1/idx);
                dCmat(i,j,k,l,:) = slmeval(phi,C{i,j,k,l},1);
                d2Cmat(i,j,k,l,:) = slmeval(phi,C{i,j,k,l},2);
            end
        end
    end
end

save('HCsmlCOV.mat','Cmat','dCmat','d2Cmat');

figure(6)
hold on
plot(phi.^(idx),2.*squeeze(Cmat(1,1,1,1,:)),'Color','R','LineWidth',2)
plot(phi.^(idx),2.*squeeze(Cmat(2,2,2,2,:)),'Color','B','LineWidth',2)
plot(phi.^(idx),2.*squeeze(Cmat(2,2,1,1,:)),'Color','G','LineWidth',2)
plot(phi.^(idx),2.*squeeze(Cmat(1,2,1,2,:)),'Color','K','LineWidth',2)
lgd = legend('C1111','C2222','C2211','C1212')
title('Change of variable in phi values')
fontsize(lgd,14,'points')
box on
grid on
xlabel('$\varphi = \phi^2$ [-]','Interpreter','latex')
% exportgraphics(gcf,'transparent.eps',...   % since R2020a
%     'ContentType','vector',...
%     'BackgroundColor','none')

figure(7)
hold on
plot(phi,squeeze(CmatCOV(1,1,1,1,:)),'Color','R','LineWidth',2)
plot(phi,squeeze(CmatCOV(2,2,2,2,:)),'Color','B','LineWidth',2)
plot(phi,squeeze(CmatCOV(2,2,1,1,:)),'Color','G','LineWidth',2)
plot(phi,squeeze(CmatCOV(1,2,1,2,:)),'Color','K','LineWidth',2)
legend('C1111','C2222','C2211','C1212')
title('Change of variable in C values')

figure(8)
hold on
plot(phi,2.*squeeze(Cmat(1,1,1,1,:)),'Color','R','LineWidth',2)
plot(phi,2.*squeeze(Cmat(2,2,2,2,:)),'Color','B','LineWidth',2)
plot(phi,2.*squeeze(Cmat(2,2,1,1,:)),'Color','G','LineWidth',2)
plot(phi,2.*squeeze(Cmat(1,2,1,2,:)),'Color','K','LineWidth',2)
lgd = legend('C1111','C2222','C2211','C1212')
title('Original')


%%

figure(4)
hold on
plot(phi,2.*squeeze(mat(1,1,1,1,:)),'.','Color','R')
plot(phi,2.*squeeze(mat(2,2,2,2,:)),'.','Color','B')
plot(phi,2.*squeeze(mat(2,2,1,1,:)),'.','Color','G')
plot(phi,2.*squeeze(mat(1,2,1,2,:)),'.','Color','K')
legend('C1111','C2222','C2211','C1212')

figure(5)
hold on
fplot(@(x) 2.*degradation.fun{1,1,1,1}(x),[0,1],'Color','R')
fplot(@(x) 2.*degradation.fun{2,2,2,2}(x),[0,1],'Color','B')
fplot(@(x) 2.*degradation.fun{2,2,1,1}(x),[0,1],'Color','G')
fplot(@(x) 2.*degradation.fun{1,2,1,2}(x),[0,1],'Color','K')
legend('C1111','C2222','C2211','C1212')

figure(6)
hold on
plot(phi,2.*squeeze(Cmat(1,1,1,1,:)),'Color','R','LineWidth',2)
plot(phi,2.*squeeze(Cmat(2,2,2,2,:)),'Color','B','LineWidth',2)
plot(phi,2.*squeeze(Cmat(2,2,1,1,:)),'Color','G','LineWidth',2)
plot(phi,2.*squeeze(Cmat(1,2,1,2,:)),'Color','K','LineWidth',2)
lgd = legend('C1111','C2222','C2211','C1212')
fontsize(lgd,14,'points')
box on
grid on
xlabel('$\phi$ [-]','Interpreter','latex')
exportgraphics(gcf,'transparent.eps',...   % since R2020a
    'ContentType','vector',...
    'BackgroundColor','none')

figure(100)
hold on
fplot(@(x) 2.*degradation.fun{1,1,1,1}(x),[0,1],'--','Color','R')
fplot(@(x) 2.*degradation.fun{2,2,2,2}(x),[0,1],'--','Color','B')
fplot(@(x) 2.*degradation.fun{2,2,1,1}(x),[0,1],'--','Color','G')
fplot(@(x) 2.*degradation.fun{1,2,1,2}(x),[0,1],'--','Color','K')
plot(phi,2.*squeeze(Cmat(1,1,1,1,:)),'Color','R')
plot(phi,2.*squeeze(Cmat(2,2,2,2,:)),'Color','B')
plot(phi,2.*squeeze(Cmat(2,2,1,1,:)),'Color','G')
plot(phi,2.*squeeze(Cmat(1,2,1,2,:)),'Color','K')
legend('C1111 (fitting)','C2222','C2211','C1212','C1111 (splines)','C2222','C2211','C1212')

figure(7)
hold on
plot(phi,squeeze(dCmat(1,1,1,1,:)),'Color','R')
plot(phi,squeeze(dCmat(2,2,2,2,:)),'Color','B')
plot(phi,squeeze(dCmat(2,2,1,1,:)),'Color','G')
plot(phi,squeeze(dCmat(1,2,1,2,:)),'Color','K')
legend('dC1111','dC2222','dC2211','dC1212')

figure(8)
hold on
plot(phi,squeeze(d2Cmat(1,1,1,1,:)),'Color','R')
plot(phi,squeeze(d2Cmat(2,2,2,2,:)),'Color','B')
plot(phi,squeeze(d2Cmat(2,2,1,1,:)),'Color','G')
plot(phi,squeeze(d2Cmat(1,2,1,2,:)),'Color','K')
legend('dC1111','dC2222','dC2211','dC1212')




