function PlottinMaxWithPNorm

p = 64;

%nIteration = 498;
%fCase = 'ExperimentingPlot';
%folder = ['/home/alex/git-repos/Swan/Output/',fCase];

%nIteration = 545;
%fCase = 'LatticeExperimentInputCantileverSymmetricMeshSuperEllipsePDE';
%folder = '/media/alex/My Passport/LatticeResults/CantileverSymmetricMeshSuperEllipsePDEEpsilonEhGradient2000';

nIteration = 235;
fCase = 'LatticeExperimentInputCantileverSymmetricMeshSuperEllipsePDE';
folder = '/media/alex/My Passport/LatticeResults/CantileverSymmetricMeshSuperEllipsePDEGradientEpsilonH';

%nIteration = 2362;
%fCase = 'LatticeExperimentInputCantileverSymmetricMeshSuperEllipsePDE';
%folder = '/media/alex/My Passport/LatticeResults/CantileverSymmetricMeshSuperEllipsePDE10000iteration';
for iter = 1:nIteration
    
   
    s.fileName = [fCase,num2str(iter)];
    s.folderPath = fullfile(folder);
    
    wM = WrapperMshResFiles(s);
    wM.compute();
    
    mesh = wM.mesh;
    
    quad = Quadrature.set(mesh.type);
    quad.computeQuadrature('CONSTANT');
    
    dvolum = mesh.computeDvolume(quad);
    sl2Norm = wM.dataRes.StressNormGauss;
    sl2NormP = (sl2Norm).^(p);
    int = sl2NormP.*dvolum';
    intOpt = int(:);
    stressNorml2LpNorm = sum(intOpt);
    stresLpNorm(iter)  = stressNorml2LpNorm.^(1/p);
    stressMax(iter)    = max(abs(sl2Norm));
    iter
    
    if iter == 1 || iter/10 == floor(iter/10)
        print(folder,fCase,stresLpNorm,stressMax,iter)
    end
    
    
end
print(folder,fCase,stresLpNorm,stressMax,iter)
end


function print(folder,fCase,stresLpNorm,stressMax,nIteration)

fid = fopen(fullfile(folder,[fCase,'.txt']));
tLine = textscan(fid,'%s','delimiter','\n', 'headerlines',0);

a = tLine{1};
costA = a(end-9,1);
costAs = split(costA);
costAt = str2double(costAs(2));

a = tLine{1};
costR = a(end,1);
costRs = split(costR);
costRt = str2double(costRs(2));

const = costRt/costAt;

fNameMon = fullfile(folder,'Monitoring.fig');
h = openfig(fNameMon);
handles = findobj(h,'Type','line');
iterC = get(handles(5),'Xdata');
cost = get(handles(5),'Ydata');
close all
f = figure();
pN{1} = plot(iterC,cost*const);
hold on
pN{2} = plot(1:nIteration,stresLpNorm);
hold on
pN{3} = plot(1:nIteration,stressMax);

leg = legend({'Cost function','L^{64}(\Omega) norm','max norm'});
set(leg);

savefig(f,fullfile(folder,[fCase,'CostVsLinf.fig']))

p = plotPrinter(f,pN);
p.print(fullfile(folder,[fCase,'Cost']));

end