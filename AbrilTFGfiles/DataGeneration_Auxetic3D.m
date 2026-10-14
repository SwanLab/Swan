%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS

p.Training   = 'EIFEM';      % 'EIFEM'/'Multiscale', ('EIFisol')
p.Sampling   = 'Isolated';  %'Isolated'/'Oversampling'
p.Inclusion  = 'Hole';        % 'Material'/'Hole'/'MeshRaul'
p.nelem      = 20;
meshName     = p.nelem+"x"+p.nelem;

%% DATA GENERATION

mR           = createReferenceMesh(p);
switch p.Training
    case 'Multiscale'
        p.Sampling      = 'Isolated';
        material  = createMaterial(mR,[1 1]);
        mesh            = mR;
        s.type          = 'discontinuous';
        switch s.type
            case 'continuous'
                bMesh      = mesh.createSingleBoundaryMesh();
                s.lambdaFun     = LagrangianFunction.create(bMesh,mesh.ndim, 'P1');
            case 'discontinuous'
                bMesh        = mesh.createBoundaryMesh();
                for i=1:numel(bMesh)
                    s.lambdaFun{i}=LagrangianFunction.create(bMesh{i}.mesh,mesh.ndim, 'P1');
                end
        end
        s.mesh          = mesh;
        s.bMesh         = bMesh;
        s.material      = material;
        s.uFun          = LagrangianFunction.create(mesh, mesh.ndim, 'P1');
        s.dirichletFun  = createDirichletFunction(bMesh,s.type);
        e  = ElasticHarmonicExtension(s);
        [T,lambda,K,Kcoarse] = e.solve();
        % V  = Integrator.compute(mT.designVariable.fun,mR,2);

    case {'EIFEM','EIFisol'}
        [nS,dI]      = defineNumberOfSubdomains(p.Sampling);
        [material]    = createMaterial(mR,nS);
        s.mesh           = mR;
        s.material       = material;
        s.domainIndices  = dI;
        s.nSubdomains    = nS;
        s.type           = 'continuous';
        s.Coarseorder    = 1;
        m= EIFEMTraining(s);
        data              = m.train();
        data.dirac        = true;
        data.type         = s.type;
        data.Coarseorder  = 1;
        data.material     = createMaterial(mR,[1 1]);
        z = OfflineDataProcessor(data);

        EIFEoper = z.computeROMbasis();
        T        = EIFEoper.U;
        mesh     = data.mesh;
        Kcoarse  = EIFEoper.Kcoarse;
        vT       = mR.computeVolume();
end


string = "Auxetic3DEIFEM";

% Guarda el .mat per cert radi
FileName=fullfile('AbrilTFGfiles','Data',"Auxetic",string);

switch p.Training
    case 'Multiscale'
        save(FileName,"T","Kcoarse","mesh");
    case {'EIFEM','EIFisol'}
        save(FileName, "EIFEoper","T","Kcoarse","mesh");
end
    


%% FUNCTIONS

function mS = createReferenceMesh(p)
    switch p.Inclusion
        case 'Material'
            mS = createStructuredMesh(p);
        case 'Hole'
            mS = createStructuredMesh2D(p);
            lvSet    = createLevelSetFunction(mS);
            uMesh    = computeUnfittedMesh(mS,lvSet);
            mR1       = uMesh.createInnerMesh();
            mR2       = SmoothMesh(mR1);

            [nodes3D,tets]=extrudeTriMeshToTet(mR2.coord,mR2.connec,1,20);
            s.coord=nodes3D;
            s.connec=tets;

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
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[delta,0,0];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)-[delta,0,0];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)+[delta,0,0];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[delta,0,0];

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

        case 'MeshRaul'
            filename = 'DEF_Q4auxL_1.mat';
            load(filename,'EIFEoper');
            s.coord    = EIFEoper.MESH.COOR;
            s.connec   = EIFEoper.MESH.CN;
    
            obj.xmin = min(s.coord(:,1));
            obj.xmax = max(s.coord(:,1));
            obj.ymin = min(s.coord(:,2));
            obj.ymax = max(s.coord(:,2));
    
            delta = 0;
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-delta];
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[delta,-delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[delta,delta];
            mR = Mesh.create(s);

            % mS = ExtrudeMesh(mR);
            [nodes3D,tets]=extrudeTriMeshToTet(mR.coord,mR.connec,1,20);
            s.coord=nodes3D;
            s.connec=tets;
            mS = Mesh.create(s);

    end
end

function mS=createStructuredMesh(p)
    n =p.nelem;
    m = TetraMesh(1.5,1,1,n,n,n);

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

    delta=1E-9;
    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[delta,0,0];

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)-[delta,0,0];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)+[delta,0,0];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[delta,0,0];

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

function mS = createStructuredMesh2D(p)
    n =p.nelem;
    x1      = linspace(-1.5,1.5,n);
    x2      = linspace(-1,1,n);
    [xv,yv] = meshgrid(x1,x2);
    [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
    s.coord  = V(:,1:2);
    s.connec = F;
    % mS=QuadMesh(1.5,1,3/2*n,n);
    mS = Mesh.create(s);
end

function levelSet = createLevelSetFunction(bgMesh)
    gPar.type           = 'Auxetic3D';
    gPar.length         = 2;
    gPar.height         = 2;
    gPar.xCoorCenter    = 0;
    gPar.yCoorCenter    = 0;
    gPar.theta          = 58;
    gPar.beta           = 64;
    gPar.thickness      = 0.3;
    
    g         = GeometricalFunction(gPar);
    phiFun    = g.computeLevelSetFunction(bgMesh);
    levelSet  = phiFun.fValues;
end

function uMesh=computeUnfittedMesh(bgMesh,levelSet)
    sUm.backgroundMesh = bgMesh;
    sUm.boundaryMesh   = bgMesh.createBoundaryMesh();
    uMesh              = UnfittedMesh(sUm);
    uMesh.compute(levelSet);
end

function dF = createDirichletFunction(bMesh,type)
s.mesh = bMesh;
s.order= 1;
s.type = type;
cf = CoarseFunctions(s);
dF = cf.getAnalytical();
end	

function [nS,dI] = defineNumberOfSubdomains(type)
    switch type
        case 'Isolated'
            nS = [1 1 1]; %nx ny
            dI = [1 1 1];
        case 'Oversampling'
            nS = [2 2 1]; %nx ny
            dI = [1 1 1];
    end
end


function [material] = createMaterial(mesh,nSubdomains)
    mD = createMeshDomain(nSubdomains,mesh);
    E       = 1;
    nu      = 1/3;
    young   = ConstantFunction.create(E,mD);
    poisson = ConstantFunction.create(nu,mD);
    s.type          = 'ISOTROPIC';
    s.ptype         = 'ELASTIC';
    s.ndim          = mesh.ndim;
    s.young         = young;
    s.poisson       = poisson;
    material        = Material.create(s);
end

function mD = createMeshDomain(nS,mesh)
    if sum(nS > 1)>= 1
        s.nsubdomains   = nS; %nx ny
        s.meshReference = mesh;
        s.tolSameNode   = 1e-11;
        m = MeshCreatorFromRVE.create(s);
        [mD,~,~,~,~,~,~] = m.create();
    else
        mD = mesh;
    end
end


function mS = SmoothMesh(mR)
    s.coord    = mR.coord;
    s.connec   = mR.connec;

    xmin = min(s.coord(:,1));
    xmax = max(s.coord(:,1));
    ymin = min(s.coord(:,2));
    ymax = max(s.coord(:,2));

    tol=1e-10;

    leftNodes   = find(abs(s.coord(:,1)-xmin) < tol);
    rightNodes  = find(abs(s.coord(:,1)-xmax) < tol);
    bottomNodes = find(abs(s.coord(:,2)-ymin) < tol);
    topNodes    = find(abs(s.coord(:,2)-ymax) < tol);

    [~,iL] = sort(s.coord(leftNodes,2));
    leftNodes = leftNodes(iL);
    [~,iR] = sort(s.coord(rightNodes,2));
    rightNodes = rightNodes(iR);
    [~,iB] = sort(s.coord(bottomNodes,1));
    bottomNodes = bottomNodes(iB);
    [~,iT] = sort(s.coord(topNodes,1));
    topNodes = topNodes(iT);
    y = 0.5*(s.coord(leftNodes,2)+s.coord(rightNodes,2));
    s.coord(leftNodes,2)  = y;
    s.coord(rightNodes,2) = y;

    x = 0.5*(s.coord(bottomNodes,1)+s.coord(topNodes,1));
    s.coord(bottomNodes,1) = x;
    s.coord(topNodes,1)    = x;
    mS = Mesh.create(s);
end



function mS= ExtrudeMesh(mR)
    nelZ = 20; % nº de elementos en Z
    nNodes = size(mR.coord,1);
    nElem  = size(mR.connec,1);
    zLayers = linspace(-1,1,nelZ+1);
    coord3D = zeros(nNodes*(nelZ+1),3);

    for k = 1:(nelZ+1)
        offset = (k-1)*nNodes;
        coord3D(offset+1:offset+nNodes,1:2) = mR.coord;
        coord3D(offset+1:offset+nNodes,3)   = zLayers(k);
    end

    connec3D = zeros(nElem*nelZ,8);
    e = 1;

    for k = 1:nelZ
        lower = (k-1)*nNodes;
        upper = k*nNodes;
        for ie = 1:nElem
            q = mR.connec(ie,:);
            connec3D(e,:) = [ q(1)+lower ...
                              q(2)+lower ...
                              q(3)+lower ...
                              q(4)+lower ...
                              q(1)+upper ...
                              q(2)+upper ...
                              q(3)+upper ...
                              q(4)+upper ];
                              e = e + 1;
        end
    end
    s.coord=coord3D;
    s.connec=connec3D;
    mS = Mesh.create(s);
end

function [nodes3D,tets] = extrudeTriMeshToTet(nodes2D,tri,thickness,nLayers)
    
    N = size(nodes2D,1);
    z = linspace(0,thickness,nLayers+1);
    nodes3D = zeros((nLayers+1)*N,3);
    
    for k=1:nLayers+1
        idx = (k-1)*N + (1:N);
        nodes3D(idx,1:2)=nodes2D;
        nodes3D(idx,3)=z(k);
    end
    
    tets = zeros(3*size(tri,1)*nLayers,4);
    e = 1;
    
    for layer = 1:nLayers
    
        b = (layer-1)*N;
        t = layer*N;
    
        for k=1:size(tri,1)
    
            v = tri(k,:);
            % Bottom
            a=v(1)+b;
            b1=v(2)+b;
            c=v(3)+b;
    
            % Top
            A=v(1)+t;
            B=v(2)+t;
            C=v(3)+t;
    
            % Global ordering
            ids = [a b1 c];
            [~,perm]=sort(ids);
    
            switch perm(1)
                case 1
                    % vertex a is smallest
                    tets(e,:)=[a b1 c C];
                    e=e+1;
    
                    tets(e,:)=[a b1 B C];
                    e=e+1;
    
                    tets(e,:)=[a A B C];
                    e=e+1;
    
                case 2
                    % vertex b is smallest
                    tets(e,:)=[b1 c a A];
                    e=e+1;
    
                    tets(e,:)=[b1 c C A];
                    e=e+1;
    
                    tets(e,:)=[b1 B C A];
                    e=e+1;
    
                case 3
                    % vertex c is smallest
                    tets(e,:)=[c a b1 B];
                    e=e+1;
                    tets(e,:)=[c a A B];
                    e=e+1;
                    tets(e,:)=[c C A B];
                    e=e+1;
            end
        end
    end

    for e = 1:size(tets,1)

        X = nodes3D(tets(e,:),:);

        J = [ (X(2,:)-X(1,:))' ...
              (X(3,:)-X(1,:))' ...
              (X(4,:)-X(1,:))' ];

        if det(J) < 0
            % Intercambiar dos nodos invierte la orientación
            tets(e,[2 3]) = tets(e,[3 2]);
        end
    end

end