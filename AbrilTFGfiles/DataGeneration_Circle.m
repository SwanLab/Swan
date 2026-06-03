%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS
% r=1e-6:0.05:0.999; 
% r=1e-6:0.1:0.999; 
% r=0:0.05:0.999;
% r=0.2:0.05:0.6;
r=0.8;

p.Training   = 'Multiscale';      % 'EIFEM'/'Multiscale'
p.Inclusion  = 'Hole';        %'Material'/'Hole'/'HoleRaul'
p.Sampling   = 'Isolated';  %'Isolated'/'Oversampling'
p.nelem      = 20;
meshName     = p.nelem+"x"+p.nelem;


%% DATA GENERATION

for j = 1:size(r,2)
    radius = r(j);
    mR              = createReferenceMesh(p,radius);
    g=computeLevelSet(radius);
    switch p.Training
        case 'Multiscale'
            p.Sampling = 'Isolated';
            material   = createMaterial(mR,[1 1],p.Inclusion,g);
            mesh       = mR;
            bMesh      = mesh.createSingleBoundaryMesh();
            s.mesh          = mesh;
            s.uFun          = LagrangianFunction.create(mesh, mesh.ndim, 'P1');
            s.lambdaFun     = LagrangianFunction.create(bMesh,mesh.ndim, 'P1');
            s.material      = material;
            s.dirichletFun  = createDirichletFunction(bMesh);
            e  = ElasticHarmonicExtension(s);
            [T,lambda,K,Kcoarse] = e.solve();

        case 'EIFEM'
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
            data.material = createMaterial(mR,[1 1],p.Inclusion,g);
            data.dirac=true;
            z = OfflineDataProcessor(data);

            EIFEoper = z.computeROMbasis();
            T        = EIFEoper.U;
            mesh     = data.mesh;
            Kcoarse  = EIFEoper.Kcoarse;
            
    end
    R        = r(j);

    % Initialization for K_all and T_all
    if j==1
        K_all=zeros([size(Kcoarse), length(r)]);
        T_all=zeros(mesh.nnodes,size(T,2)*mesh.ndim+mesh.ndim+1,length(r));
    end

    K_all(:,:,j)=Kcoarse; 

    % Reshapes U data and adds coordinates  % Adds the radius and coordinates column
    t_all=[];
    for k = 1:size(T,2)
        t_k = reshape(T(:,k), mesh.ndim, []).';
        t_all = [t_all, t_k];
    end

    t_aux = [r(j)*ones(size(mesh.coord,1),1), mesh.coord, t_all];
    T_all(:,:,j)=t_aux;   % Saves the result for each radius


    %Designa un nom per cada linea corresponent a un radi
    meshName=p.nelem+"x"+p.nelem;
    string = strrep("r"+num2str(r(j), '%.4f'), ".", "_")+"-"+meshName+".mat";

    % Guarda el .mat per cert radi
    FileName=fullfile('AbrilTFGfiles','Data',"Circle",p.Training,p.Inclusion,p.Sampling,meshName,string);

    % FileName=fullfile('AbrilTFGfiles','Data','CircleDirac',string);


    switch p.Training
        case 'Multiscale'
            save(FileName,"T","Kcoarse","mesh","R"); 
        case 'EIFEM'
            save(FileName, "EIFEoper","T","Kcoarse","mesh","R"); 
    end

    
end



%% Reshapes the T data and saves it in a csv file
% 
% Redimensioning the U_all1
TData=[];
for n=1:size(T_all,3)
    TData=[TData;T_all(:,:,n)];
end

T=array2table(TData,"VariableNames",{'r','x','y','Tx1','Ty1','Tx2','Ty2','Tx3','Ty3','Tx4','Ty4' ...
    'Tx5','Ty5','Tx6','Ty6','Tx7','Ty7','Tx8','Ty8'});

uFileName = fullfile('AbrilTFGfiles','Data','Circle',p.Training,p.Inclusion,p.Sampling,'DataT.csv');
writematrix(TData,uFileName);


%% Reshapes the K data and saves it in a csv file

kdata=zeros(size(r,2),36);
for n=1:size(r,2)
    triangSup=triu(K_all(:,:,n));  %gets the triangular superior matrix
    clear row;
    row=[];
    for i=1:8
        for j=i:8
            row(end+1)=triangSup(i,j);
        end
    end
    kdata(n,:)=row;
end

kdata=[r.',kdata];
kFileName = fullfile('AbrilTFGfiles','Data',"Circle",p.Training,p.Inclusion,p.Sampling,'dataK.csv');
writematrix(kdata,kFileName);


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
        s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-delta];
    s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
        s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,+delta];
    s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
        s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[+delta,-delta];
    s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
        s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[+delta,+delta];
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

function dF = createDirichletFunction(bMesh)
s.mesh = bMesh;
s.type = 'continuous';
s.order = 2;
cf = CoarseFunctions(s);
dF = cf.getAnalytical();
end	
