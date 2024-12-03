clear all;  clc; close all

mesh20=1;
mesh50=1;
mesh100=1;

if mesh20==1
load('./PGDchol_20x20.mat')
residual{1}=PGDchol.residual;
err{1}=PGDchol.err;
errA{1}=PGDchol.errA;

load('./PGDjacobi_20x20.mat')
residual{2}=PGDjacobi.residual;
err{2}=PGDjacobi.err;
errA{2}=PGDjacobi.errA;

load('./PGDssor_20x20.mat')
residual{3}=PGDssor.residual;
err{3}=PGDssor.err;
errA{3}=PGDssor.errA;

load('./PGDmodes_20x20.mat')
residual{4}=PGDbasis.residual;
err{4}=PGDbasis.err;
errA{4}=PGDbasis.errA;

load('./PGDmodes20_20x20.mat')
residual{5}=PGDbasis20.residual;
err{5}=PGDbasis20.err;
errA{5}=PGDbasis20.errA;

load('./GD_20x20.mat')
residual{6}=GD.residual;
err{6}=GD.err;
errA{6}=GD.errA;

tit='Mesh 20x20';

elseif mesh50==1
load('./PGCchol_50x50.mat')
residual{1}=PCGchol.residual;
err{1}=PCGchol.err;
errA{1}=PCGchol.errA;

load('./PGCjacobi_50x50.mat')
residual{2}=PCGjacobi.residual;
err{2}=PCGjacobi.err;
errA{2}=PCGjacobi.errA;

load('./PGCssor_50x50.mat')
residual{3}=PCGssor.residual;
err{3}=PCGssor.err;
errA{3}=PCGssor.errA;

load('./PGCmodes_50x50.mat')
residual{4}=PCGmodes.residual;
err{4}=PCGmodes.err;
errA{4}=PCGmodes.errA;

load('./PGCmodes20_50x50.mat')
residual{5}=PCGmodes20.residual;
err{5}=PCGmodes20.err;
errA{5}=PCGmodes20.errA;

load('./GC_50x50.mat')
residual{6}=CG.residual;
err{6}=CG.err;
errA{6}=CG.errA;

tit='Mesh 50x50';

elseif mesh100==1
load('./PGCchol_100x100.mat')
residual{1}=PCGchol.residual;
err{1}=PCGchol.err;
errA{1}=PCGchol.errA;

load('./PGCjacobi_100x100.mat')
residual{2}=PCGjacobi.residual;
err{2}=PCGjacobi.err;
errA{2}=PCGjacobi.errA;

load('./PGCssor_100x100.mat')
residual{3}=PCGssor.residual;
err{3}=PCGssor.err;
errA{3}=PCGssor.errA;

load('./PGCmodes_100x100.mat')
residual{4}=PCGmodes.residual;
err{4}=PCGmodes.err;
errA{4}=PCGmodes.errA;

load('./PGCmodes20_100x100.mat')
residual{5}=PCGmodes20.residual;
err{5}=PCGmodes20.err;
errA{5}=PCGmodes20.errA;

load('./GC_100x100.mat')
residual{6}=CG.residual;
err{6}=CG.err;
errA{6}=CG.errA;

tit='Mesh 100x100';
end

figure
t=tiledlayout(3,1);
title(t,tit)
nexttile
for i=1:length(residual)
    plot(residual{i},'LineWidth',1.3)
    hold on
end
set(gca, 'YScale', 'log')
xlabel('iteration')
ylabel('residual')
legend('Cholesky','Jacobi','SSOR','Eig 10 Basis','Eig 20 Basis','CG')

nexttile
for i=1:length(err)
    plot(err{i},'LineWidth',1.3)
    hold on
end
set(gca, 'YScale', 'log')
xlabel('iteration')
ylabel('||x-x*||')
legend('Cholesky','Jacobi','SSOR','Eig 10 Basis','Eig 20 Basis','CG')

nexttile
for i=1:length(errA)
    plot(errA{i},'LineWidth',1.3)
    hold on
end
set(gca, 'YScale', 'log')
xlabel('iteration')
ylabel('||x-x*||_{A}')
legend('Cholesky','Jacobi','SSOR','Eig 10 Basis','Eig 20 Basis','CG')

