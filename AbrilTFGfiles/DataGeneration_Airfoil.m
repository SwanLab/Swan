%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS

p.Training   = 'EIFEM';      % 'EIFEM'/'Multiscale'
p.Sampling   = 'Oversampling';  %'Isolated'/'Oversampling'

%% DATA GENERATION
mR              = createReferenceMesh();
switch p.Training
    case 'Multiscale'
        material   = createMaterial(mR);
        mesh       = mR;
        bMesh      = mesh.createSingleBoundaryMesh();
        s.mesh          = mesh;
        s.uFun          = LagrangianFunction.create(mesh, mesh.ndim, 'P1');
        s.lambdaFun     = LagrangianFunction.create(bMesh,mesh.ndim, 'P1');
        s.material      = material;
        s.dirichletFun  = createDirichletFunction(bMesh);
        e  = ElasticHarmonicExtension(s);
        [T,lambda,K,Kcoarse] = e.solve();

        if strcmpi(p.Sampling, 'Oversampling')
            [nS,dI]      = defineNumberOfSubdomains(p.Sampling);
            mD           = createMeshDomain(mR,nS);
            [material]     = createMaterial(mD);
            s.mesh           = mR;
            s.material       = material;
            s.domainIndices  = dI;
            s.nSubdomains    = nS;
            m= EIFEMTraining(s);
            data          = m.train();
            [data.material] = createMaterial(mR);

            Kf = 0.5 * (K + K');
            % Normalization
            for i=1:size(T,2)
                phi=T(:,i);
                val = phi' * Kf * phi;
                if val < 1e-12
                    warning(['Modo ', num2str(i), ' casi rígido o mal condicionado']);
                    val = 1e-12;  % evitar división por cero o sqrt negativo
                end
                T(:,i)=phi/sqrt(val);
            end
            T = real(T);
            Kcoarse=T'*K*T;
            Kcoarse=0.5 * (Kcoarse + Kcoarse');
        end

    case 'EIFEM'
        [nS,dI]      = defineNumberOfSubdomains(p.Sampling);
        mD           = createMeshDomain(mR,nS);
        [material]     = createMaterial(mD);
        s.mesh           = mR;
        s.material       = material;
        s.domainIndices  = dI;
        s.nSubdomains    = nS;
        m= EIFEMTraining(s);
        data          = m.train();
        [data.material] = createMaterial(mR);
        data.dirac=false;
        z = OfflineDataProcessor(data);

        EIFEoper = z.computeROMbasis();
        T        = EIFEoper.U;
        mesh     = data.mesh;
        Kcoarse  = EIFEoper.Kcoarse;
end

string = "Airfoil_Isolated2.mat";

% Guarda el .mat per cert radi
FileName=fullfile('AbrilTFGfiles','Data',"Airfoil",string);

switch p.Training
    case 'Multiscale'
        save(FileName,"T","Kcoarse");
    case 'EIFEM'
        save(FileName, "EIFEoper","T","Kcoarse","mesh");
end



%% FUNCTIONS

function mS = createReferenceMesh()

    % file = 'meshAirfoilTetra.m';
    % a.fileName = file;
    % s = FemDataContainer(a);
    % mS = s.mesh;
    filename = 'DEF_Q8_wing_1.mat';
    load(filename);
    s.coord    = EIFEoper.MESH.COOR;
    % s.coord = [s.coord(:,1),s.coord(:,3),s.coord(:,2)];
    s.connec   = EIFEoper.MESH.CN;

    maxC= max(s.coord);
    minC = min(s.coord);

    obj.xmin = minC(1);
    obj.xmax = maxC(1);
    obj.ymin = minC(2);
    obj.ymax = maxC(2);
    obj.zmin = minC(3);
    obj.zmax = maxC(3);

    delta=1E-9;

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)+[0,0,delta];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[0,0,delta];

    s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

    s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

    s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

    s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:)-[0,delta,0];

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:)+[0,delta,0];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:)-[0,delta,0];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:)+[0,delta,0];

    mS = Mesh.create(s);

end

function dF = createDirichletFunction(bMesh)
s.mesh = bMesh;
s.order= 1;
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
            nS = [2 1 1]; %nx ny
            dI = [1 1 1];
    end
end


function [material] = createMaterial(mesh)
    E       = 1;
    nu      = 1/3;
    young   = ConstantFunction.create(E,mesh);
    poisson = ConstantFunction.create(nu,mesh);
    s.type          = 'ISOTROPIC';
    s.ptype         = 'ELASTIC';
    s.ndim          = mesh.ndim;
    s.young         = young;
    s.poisson       = poisson;
    material        = Material.create(s);
end


function meshDom=createMeshDomain(mR,nS)
    if sum(nS > 1)>= 1
        s.nsubdomains   = nS; %nx ny
        s.meshReference = mR;
        s.tolSameNode   = 1e-11;
        m = MeshCreatorFromRVE.create(s);
        [meshDom,~,~,~,~,~,~] = m.create();
    else
        meshDom = mR;
    end
end