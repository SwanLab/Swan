%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS

r=1e-6:0.05:0.96;
% r=1e-6;

p.Training   = 'EIFEM';      % 'EIFEM'/'Multiscale'
p.Sampling   = 'Oversampling';  %'Isolated'/'Oversampling'
p.nelem      = 5;
meshName     = p.nelem+"x"+p.nelem;

%% DATA GENERATION
for j = 1:size(r,2)
    radius = r(j);
        mR              = createReferenceMesh(p);
        g=computeLevelSet(radius);
        switch p.Training
            case 'Multiscale'
                p.Sampling = 'Isolated';
                material   = createMaterial(mR,[1 1],'Material',g);
                mesh       = mR;
                bMesh      = mesh.createSingleBoundaryMesh();
                s.mesh          = mesh;
                s.uFun          = LagrangianFunction.create(mesh, mesh.ndim, 'P1');
                s.lambdaFun     = LagrangianFunction.create(bMesh,mesh.ndim, 'P1');
                s.material      = material;
                s.dirichletFun  = createDirichletFunction(bMesh);
                e  = ElasticHarmonicExtension(s);
                [T,lambda,K,Kcoarse] = e.solve();
                V  = Integrator.compute(mTr.designVariable.fun,mR,2);
    
            case 'EIFEM'
                [nS,dI]      = defineNumberOfSubdomains(p.Sampling);
                [material,mT]    = createMaterial(mR,nS,'Material',g);
                s.mesh           = mR;
                s.material       = material;
                s.domainIndices  = dI;
                s.nSubdomains    = nS;  
                s.levelSet = mT.levelSet;
                s.unfittedMesh = mT.unfittedMesh;
                m= EIFEMTraining(s);
                data          = m.train();
                [data.material,mTr] = createMaterial(mR,[1 1],'Material',g);
                z = OfflineDataProcessor(data);
    
                EIFEoper = z.computeROMbasis();
                T        = EIFEoper.U;
                mesh     = data.mesh;
                Kcoarse  = EIFEoper.Kcoarse;
                vT       = mR.computeVolume();
                V        = Integrator.compute(mTr.designVariable.fun,mR,2);
                Vfr = V/vT;
        end
    
        R        = r(j);

        % Initialization for K_all and T_all
        if j==1
            zeros([size(Kcoarse), length(r)]);
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

        %Designa un nom per cada valor diferent del parametre
        meshName=p.nelem+"x"+p.nelem;
        string = strrep("r"+num2str(r(j), '%.4f'), ".", "_")+"-"+meshName+".mat";
    
        % Guarda el .mat per cert radi
        FileName=fullfile('AbrilTFGfiles','Data',p.Training,"Sphere",string);
    
        switch p.Training
            case 'Multiscale'
                save(FileName,"T","Kcoarse","mesh","R"); 
            case 'EIFEM'
                save(FileName, "EIFEoper","T","Kcoarse","R","V","Vfr"); 
        end
    
end

%% Reshapes the T data and saves it in a csv file
% 
% % Redimensioning the U_all1
TData=[];
for n=1:size(T_all,3)
    TData=[TData;T_all(:,:,n)];
end

% T=array2table(TData,"VariableNames",{'r','x','y','Tx1','Ty1','Tx2','Ty2','Tx3','Ty3','Tx4','Ty4' ...
%     'Tx5','Ty5','Tx6','Ty6','Tx7','Ty7','Tx8','Ty8'});

uFileName = fullfile('AbrilTFGfiles','Data',p.Training,'Sphere','DataT.csv');
writematrix(TData,uFileName);


%% Reshapes the K data and saves it in a csv file

kdata=zeros(size(r,2),300);
for n=1:size(r,2)
    triangSup=triu(K_all(:,:,n));  %gets the triangular superior matrix
    clear row;
    row=[];
    for i=1:24
        for j=i:24
            row(end+1)=triangSup(i,j);
        end
    end
    kdata(n,:)=row;
end

kdata=[r.',kdata];
kFileName = fullfile('AbrilTFGfiles','Data',p.Training,'Sphere','dataK.csv');
writematrix(kdata,kFileName);



%% FUNCTIONS

function mS = createReferenceMesh(p)
    n =p.nelem;
    m = TetraMesh(1,1,1,n,n,n);

    s.coord=m.coord;
    s.connec=m.connec;
    
    maxC= max(s.coord);
    minC = min(s.coord);
    
    obj.xmin = minC(1);
    obj.xmax = maxC(1);
    obj.ymin = minC(2);
    obj.ymax = maxC(2);
    obj.zmin = minC(3);
    obj.zmax = maxC(3);

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[0,0,1E-9];

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)+[0,0,1E-9];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)-[0,0,1E-9];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[0,0,1E-9];

    s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,1E-9];

    s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,1E-9];

    s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,1E-9];

    s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,1E-9];

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:)-[0,1E-9,0];

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:)+[0,1E-9,0];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:)-[0,1E-9,0];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:)+[0,1E-9,0];

    mS = Mesh.create(s);

end

function dF = createDirichletFunction(bMesh)
s.mesh = bMesh;
s.type = 'continuous';
cf = CoarseFunctions(s);
dF = cf.getAnalytical();
end	


function [nS,dI] = defineNumberOfSubdomains(type)
    switch type
        case 'Isolated'
            nS = [1 1 1]; %nx ny
            dI = [1 1 1];
        case 'Oversampling'
            nS = [3 3 3]; %nx ny
            dI = [2 2 2];
    end
end

function g=computeLevelSet(r)
    gPar.type         = 'Sphere';
    gPar.radius       = r;
    gPar.xCoorCenter  = 0;
    gPar.yCoorCenter  = 0;
    gPar.zCoorCenter  = 0;
    g                 = GeometricalFunction(gPar);
end

function [material,m] = createMaterial(mesh,nSubdomains,inclusionType,g)
    s.mesh       = mesh;
    s.inclusionType  = inclusionType;
    s.nSubdomains    = nSubdomains;
    s.geomFun       = g;
    m = MaterialTraining(s);
    material = m.create();
end