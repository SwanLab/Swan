close all;
clear;clc;

load('Gripping.mat');
clc;

mesh           = createMesh(mesh);
designVariable = createDesignVariable(mesh,fun);
clear('fun');
designVariable.fun.plot();

e    = 6:-1:1;
pTar = 0.2:0.1:1;
pp   = 1:16;

for ii = 1:size(pTar,2)
    PpValues = zeros(size(e,2),size(pp,2));
    gxValues = zeros(size(e,2),size(pp,2));
    for jj = 1:size(e,2)
        for kk = 1:size(pp,2)
            perimeter = createPerimeterConstraint(mesh,pp(kk),pTar(ii),e(jj));
            perimeter.computeFunctionAndGradient(designVariable);
            PpValues(jj,kk) = perimeter.Pp;
            gxValues(jj,kk) = perimeter.gx;
        end
    end
    createPlot(PpValues,gxValues,pp,pTar(ii),e);
end


function mesh = createMesh(m)
    s.connec = m.connec;
    s.coord  = m.coord;
    mesh     = Mesh.create(s);
end


function designVariable = createDesignVariable(mesh,f)
    sF.fValues     = f.fValues;
    sF.mesh        = mesh;
    sF.order       = 'P1';
    s.fun          = LagrangianFunction(sF);
    s.mesh         = mesh;
    s.type         = 'Density';
    s.plotting     = false;
    dens           = DesignVariable.create(s);
    designVariable = dens;
end

function perimeter = createPerimeterConstraint(mesh,p,pTarget,eps)
    s.mesh            = mesh;
    s.perimeterTarget = pTarget;
    s.p               = p;
    s.eps             = eps;
    s.gradientTest    = LagrangianFunction.create(mesh,1,'P1');
    perimeter         = PerimeterNormPFunctionalTest(s);
end

function createPlot(PpValues,gxValues,pp,pTar,e)
    fileLocation = 'C:\Users\Biel\Desktop\UNI\TFG\ResultatsNormP_Density\01. Sentit físic Ptarget';

    figure

    hold on
    for i = 1:size(e,2)
        plot(pp,PpValues(i,:),'DisplayName',sprintf('e = %d',e(i)));
    end
    grid on
    xlabel('p')
    ylabel('Perimeter p-norm')
    title(sprintf('p vs Perimeter p-norm (pTarget = %.2f)',pTar));
    legend();
    
    screenSize = get(0, 'Screensize');
    set(gcf, 'Position', [screenSize(1), screenSize(2), screenSize(3)/2, screenSize(4)/2]);
    fileName1 = fullfile(fileLocation,sprintf('Monitoring_Perimeter_Gripping_pTar%.2f_epsilon%.2f-%.2f.png',pTar,min(e),max(e)));
    fileName2 = fullfile(fileLocation,sprintf('Monitoring_Perimeter_Gripping_pTar%.2f_epsilon%.2f-%.2f.fig',pTar,min(e),max(e)));
    print(fileName1,'-dpng','-r300');
    savefig(fileName2);
    close

    figure
    hold on
    for j = 1:size(e,2)
        plot(pp,gxValues(j,:),'DisplayName',sprintf('e = %d',e(j)));
    end
    grid on
    xlabel('p')
    ylabel('Constraint p-norm')
    title(sprintf('p vs Constraint p-norm (pTarget = %.2f)',pTar));
    legend();
    
    screenSize = get(0, 'Screensize');
    set(gcf, 'Position', [screenSize(1), screenSize(2), screenSize(3)/2, screenSize(4)/2]);
    fileName1 = fullfile(fileLocation,sprintf('Monitoring_Constraint_Gripping_pTar%.2f_epsilon%.2f-%.2f.png',pTar,min(e),max(e)));
    fileName2 = fullfile(fileLocation,sprintf('Monitoring_Constraint_Gripping_pTar%.2f_epsilon%.2f-%.2f.fig',pTar,min(e),max(e)));
    print(fileName1,'-dpng','-r300');
    savefig(fileName2);
    close
end