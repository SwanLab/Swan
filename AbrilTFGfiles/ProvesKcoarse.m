%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS
% r=1e-6:0.05:0.999; 
%r=1e-6:0.1:0.999; 
r=0:0.05:0.999;
% r=0.4;

p.Inclusion  = 'Material';        %'Material'/'Hole'/'HoleRaul'
p.nelem      = 20;

K_mult        = zeros(8,8,size(r,2));
K_eife        = zeros(8,8,size(r,2));
K_eife_dirac  = zeros(8,8,size(r,2));


%% DATA GENERATION

for j = 1:size(r,2)
    radius = r(j);
    mR              = createReferenceMesh(p,radius);
    g=computeLevelSet(radius);

    % MULTISCALE
    p.Training   = 'Multiscale';
    p.Sampling = 'Isolated';
    materialIso   = createMaterial(mR,[1 1],p.Inclusion,g);
    mesh       = mR;
    bMesh      = mesh.createSingleBoundaryMesh();
    s.mesh=bMesh;
    s.type='continuous';
    cf=CoarseFunctions(s);
    cf_a = cf.getAnalytical();
    s.mesh          = mesh;
    s.uFun          = LagrangianFunction.create(mesh, mesh.ndim, 'P1');
    s.lambdaFun     = LagrangianFunction.create(bMesh,mesh.ndim, 'P1');
    s.material      = materialIso;
    s.dirichletFun  = cf_a;
    e  = ElasticHarmonicExtension(s);
    [~,~,~,Kcoarse] = e.solve();
    
    K_mult(:,:,j)=Kcoarse;

    % EIFEM
    p.Training   = 'EIFEM';
    p.Sampling   = 'Oversampling';
    [nS,dI]      = defineNumberOfSubdomains(p.Sampling);
    [material,mT]     = createMaterial(mR,nS,p.Inclusion,g);
    s.mesh           = mR;
    s.r              = radius;
    s.material       = material;
    s.domainIndices  = dI;
    s.nSubdomains    = nS;
    s.unfittedMesh = mT.unfittedMesh;
    m= EIFEMTraining(s);
    data          = m.train();
    data.material = materialIso;
    data.dirac    = false;
    z = OfflineDataProcessor(data);
    EIFEoper   = z.computeROMbasis();
    K_eife(:,:,j)  = EIFEoper.Kcoarse;

    % EIFEM DIRAC
    data.dirac    = true;
    z = OfflineDataProcessor(data);
    EIFEoper   = z.computeROMbasis();
    K_eife_dirac(:,:,j) = EIFEoper.Kcoarse;
   
end

%% IMPORT FROM DIRECTORY

files = dir(fullfile("AbrilTFGfiles/Data/LatticeDirac_3D/", 't1_0_10_*.mat'));

for k = 1:1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    load(filePath,'Kcoarse',"tCross");
    K_Over(:,:,k)=Kcoarse;
    disp(['Loaded: ', files(k).name]);  % Display the file being loaded
    t2(k)=tCross;
end

files2 = dir(fullfile("AbrilTFGfiles/Data/LatticeMult_3D/", 't1_0_1*.mat'));

for k = 1:1:length(files2)
    filePath = fullfile(files2(k).folder, files2(k).name);
    load(filePath,'Kcoarse',"tCross");
    K_Mult(:,:,k)=Kcoarse;
    disp(['Loaded: ', files2(k).name]);  % Display the file being loaded
    t2(k)=tCross;
end


pairs = [];


for i = 1:24
    for j = i:24
        pairs = [pairs; i j];
    end
end

figure
t = tiledlayout(5,6,'TileSpacing','compact','Padding','compact');
for k = 1:30
    i = pairs(k,1);
    j = pairs(k,2);

    nexttile
    hold on
    y_over  = squeeze(K_Over(i,j,:));
    y_mult  = squeeze(K_Mult(i,j,:));
    plot(t2, y_over, 'LineWidth', 1.5);
    plot(t2, y_mult, 'LineWidth', 1.5);
    legend('EIFEM over','Multiscale')
    xlabel('t2')

    title(sprintf('$K_{%d,%d}$', i, j), ...
        'Interpreter','latex','FontSize',12)
end

figure
t = tiledlayout(5,6,'TileSpacing','compact','Padding','compact');

for k = 31:60
    i = pairs(k,1);
    j = pairs(k,2);
    
    nexttile
    hold on
    y_over  = squeeze(K_Over(i,j,:));
    y_mult  = squeeze(K_Mult(i,j,:));
    plot(t2, y_over, 'LineWidth', 1.5);
    plot(t2, y_mult, 'LineWidth', 1.5);
    legend('EIFEM over','Multiscale')
    xlabel('t2')
    
    title(sprintf('$K_{%d,%d}$', i, j), ...
        'Interpreter','latex','FontSize',12)
end



%%  PLOTS
% 
% figure
% plot(r,K_mult,r,K_eife,r,K_eife_dirac,LineWidth=1.5)
% legend('Multiscale','EIFEM Oversampling', 'EIFEM oversampling +Dirac')
% xlabel('r')
% ylabel('K(2,2)')
% title('K(2,2) vs r')
idx=1;

pairs = [];

for i = 1:8
    for j = i:8
        pairs = [pairs; i j];
    end
end

figure
t = tiledlayout(3,6,'TileSpacing','compact','Padding','compact');

for k = 1:18
    i = pairs(k,1);
    j = pairs(k,2);

    nexttile
    hold on
    y_mult  = squeeze(K_mult(i,j,:));
    y_eife  = squeeze(K_eife(i,j,:));
    y_eifeD = squeeze(K_eife_dirac(i,j,:));

    plot(r, y_mult, 'LineWidth', 1.5);
    plot(r, y_eife, 'LineWidth', 1.5);
    plot(r, y_eifeD, 'LineWidth', 1.5);
    xlabel('r')
    legend('Multiscale','EIFEM', 'EIFEM+Dirac')

    title(sprintf('$K_{%d,%d}$', i, j), ...
        'Interpreter','latex','FontSize',12)
end


figure
t = tiledlayout(3,6,'TileSpacing','compact','Padding','compact');

for k = 19:36
    i = pairs(k,1);
    j = pairs(k,2);
    
    nexttile
    hold on
    y_mult  = squeeze(K_mult(i,j,:));
    y_eife  = squeeze(K_eife(i,j,:));
    y_eifeD = squeeze(K_eife_dirac(i,j,:));

    plot(r, y_mult, 'LineWidth', 1.5);
    plot(r, y_eife, 'LineWidth', 1.5);
    plot(r, y_eifeD, 'LineWidth', 1.5);
    xlabel('r')
    legend('Multiscale','EIFEM', 'EIFEM+Dirac')
    
    title(sprintf('$K_{%d,%d}$', i, j), ...
        'Interpreter','latex','FontSize',12)
end


% tiledlayout(6,6, 'TileSpacing', 'compact', 'Padding', 'compact')
% for i = 1:8
%     for j = i:8 
%         nexttile
%         hold on
%         y_mult  = squeeze(K_mult(i,j,:));
%         y_eife  = squeeze(K_eife(i,j,:));
%         y_eifeD = squeeze(K_eife_dirac(i,j,:));
% 
%         plot(r, y_mult, 'LineWidth', 1.5);
%         plot(r, y_eife, 'LineWidth', 1.5);
%         plot(r, y_eifeD, 'LineWidth', 1.5);
%         xlabel('r')
%         title(sprintf('K(%d,%d)', i, j))
%         legend('Multiscale','EIFEM Oversamp', 'EIFEM Oversamp +Dirac')
% 
%         idx = idx + 1;
%     end
% end
%% FUNCTIONS

function mS = createReferenceMesh(p,r)
    
    switch p.Inclusion
        case 'Material'
            mS = createStructuredMesh(p);
        case 'Hole'
            mS = createStructuredMesh(p);
            lvSet    = createLevelSetFunction(mS,r);
            uMesh    = computeUnfittedMesh(mS,lvSet);
            mS       = uMesh.createInnerMesh();
        case 'HoleRaul'
            switch p.nelem
                case 10
                    mS=mesh_rectangle_via_triangles(r,1,-1,1,-1,7,6,0,0);   % 10x10
                case 20
                    % mS=mesh_rectangle_via_triangles(r,1,-1,1,-1,15,12,0,0);
                    mS=mesh_rectangle_via_triangles(r,1,-1,1,-1,7,14,0,0); % 20x20
                case 50
                    mS=mesh_rectangle_via_triangles(r,1,-1,1,-1,34,35,0,0);  % 50x50
            end
    end
end


function mS = createStructuredMesh(p)
    n =p.nelem;
    x1      = linspace(-1,1,n);
    x2      = linspace(-1,1,n);
    [xv,yv] = meshgrid(x1,x2);
    [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
    s.coord  = V(:,1:2);
    s.connec = F;
    % 
    % mesh= QuadMesh(1,1,n,n);
    % s.coord= mesh.coord;
    % s.connec=mesh.connec;

    obj.xmin = min(x1);            
    obj.xmax = max(x1);
    obj.ymin = min(x2);
    obj.ymax = max(x2);

    delta = 1e-9;
    s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
        s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-0*delta];
    s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
        s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,0*delta];
    s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
        s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[delta,-0*delta];
    s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
        s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[delta,0*delta];

    mS = Mesh.create(s); 
   
end

function levelSet = createLevelSetFunction(bgMesh,r)
    sLS.type        = 'CircleInclusion';
    sLS.xCoorCenter = 0;
    sLS.yCoorCenter = 0;
    sLS.radius      = r;
    g               = GeometricalFunction(sLS);
    lsFun           = g.computeLevelSetFunction(bgMesh);
    levelSet        = lsFun.fValues;
end

function uMesh = computeUnfittedMesh(bgMesh,levelSet)
    sUm.backgroundMesh = bgMesh;
    sUm.boundaryMesh   = bgMesh.createBoundaryMesh();
    uMesh              = UnfittedMesh(sUm);
    uMesh.compute(levelSet);
end

function [nS,dI] = defineNumberOfSubdomains(type)
    switch type
        case 'Isolated'
            nS = [1 1]; %nx ny
            dI = [1 1];
        case 'Oversampling'
            nS = [5 5]; %nx ny
            dI = [3 3];
    end
end

function g=computeLevelSet(r)
    gPar.type         = 'Circle';
    gPar.radius       = r;
    gPar.xCoorCenter  = 0;
    gPar.yCoorCenter  = 0;
    g                 = GeometricalFunction(gPar);
end

function [material,m] = createMaterial(mesh,nSubdomains,inclusionType,g)
    s.mesh           = mesh;
    s.inclusionType  = inclusionType;
    s.nSubdomains    = nSubdomains;
    s.geomFun        = g;
    m = MaterialTraining(s);
    material = m.create();
end