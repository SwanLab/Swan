%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS
%r=1e-6:0.05:0.999; 
%r=1e-6:0.1:0.999; 
%r=0:0.05:0.999;
r=0.5;

p.Training   = 'EIFEM';      % 'EIFEM'/'Multiscale'
p.Inclusion  = 'Material';        %'Material'/'Hole'/'HoleRaul'
p.nelem      = 20;
meshName     = p.nelem+"x"+p.nelem;


%% DATA GENERATION

for j = 1:size(r,2)
    radius = r(j);
    mR              = createReferenceMesh(p,radius);
    switch p.Training
        case 'Multiscale'
            material = createMaterialTraining(mR, radius,[1 1],p.Inclusion);
            mesh    = mR;
            bMesh   = mesh.createSingleBoundaryMesh();
            s.mesh=bMesh;
            cf=CoarseFunctions(s);
            f = cf.getAnalytical();
            s.mesh          = mesh;
            s.uFun          = LagrangianFunction.create(mesh, mesh.ndim, 'P1');
            s.lambdaFun     = LagrangianFunction.create(bMesh,mesh.ndim, 'P1');
            s.material      = material;
            s.dirichletFun  = f;
            e  = ElasticHarmonicExtension(s);
            [T,lambda,K,Kcoarse] = e.solve();

        case 'EIFEM'
            samplingType = 'Oversampling'; %'Isolated'/'Oversampling'
            [nS,dI] = defineNumberOfSubdomains(samplingType);
            material = createMaterialTraining(mR, radius,nS,p.Inclusion);
            s.mesh           = mR;
            s.r              = radius;
            s.material       = material;
            s.domainIndices  = dI;
            s.nSubdomains    = nS;            
            m= EIFEMTraining(s);
            data          = m.train();
            data.material = createMaterialTraining(mR, radius,[1 1],p.Inclusion);
            z = OfflineDataProcessor(data);

            EIFEoper = z.computeROMbasis();
            T        = EIFEoper.U;
            mesh     = data.mesh;
            Kcoarse  = EIFEoper.Kcoarse;
            p.Inclusion = fullfile(p.Inclusion,samplingType);
            
    end
    R        = r(j);

    % Initialization for K_all and T_all
    if j==1
        K_all=zeros(8,8,length(r));
        T_all=zeros(mesh.nnodes,19,length(r));
    end

    K_all(:,:,j)=Kcoarse;   
    
    % Reshapes U data and adds coordinates
    t1=reshape(T(:,1).',2,[]).';   % Joins the Tx and Ty coeff at the same line
    t2=reshape(T(:,2).',2,[]).';                                  
    t3=reshape(T(:,3).',2,[]).';                                  
    t4=reshape(T(:,4).',2,[]).';                                  
    t5=reshape(T(:,5).',2,[]).';                                  
    t6=reshape(T(:,6).',2,[]).';   
    t7=reshape(T(:,7).',2,[]).';
    t8=reshape(T(:,8).',2,[]).';

    t_aux=[r(j)*ones(size(mesh.coord,1),1), mesh.coord, t1,t2,t3,t4,t5,t6,t7,t8];  % Adds the radius and coordinates column

    T_all(:,:,j)=t_aux;   % Saves the result for each radius


    %Designa un nom per cada linea corresponent a un radi
    meshName=p.nelem+"x"+p.nelem;
    string = strrep("r"+num2str(r(j), '%.4f'), ".", "_")+"-"+meshName+".mat";

    % Guarda el .mat per cert radi
    FileName=fullfile('AbrilTFGfiles','Data',p.Training,p.Inclusion,meshName,string);


    switch p.Training
        case 'Multiscale'
            save(FileName,"T","Kcoarse","mesh","R"); 
        case 'EIFEM'
            save(FileName, "EIFEoper","T","Kcoarse","mesh","R"); 
    end

    
end



%% Reshapes the T data and saves it in a csv file

% Redimensioning the U_all1
TData=[];
for n=1:size(T_all,3)
    TData=[TData;T_all(:,:,n)];
end

T=array2table(TData,"VariableNames",{'r','x','y','Tx1','Ty1','Tx2','Ty2','Tx3','Ty3','Tx4','Ty4' ...
    'Tx5','Ty5','Tx6','Ty6','Tx7','Ty7','Tx8','Ty8'});

uFileName = fullfile('AbrilTFGfiles','Data',p.Training,p.Inclusion,'DataT.csv');
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
kFileName = fullfile('AbrilTFGfiles','Data',p.Training,p.Inclusion,'dataK.csv');
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
                    mS=mesh_rectangle_via_triangles(r,1,-1,1,-1,15,12,0,0); % 20x20
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