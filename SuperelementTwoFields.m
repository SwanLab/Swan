clc,clear,close all 

%% CONVENTIONAL SOLVER
% Create geometry
fileName = 'CantileverAxialLoad';
s.testName    = [fileName,'.m'];
s.problemCase = 'cantilever';
s.x1          = 1;
s.y1          = 0.5;
s.N           = 49;
s.M           = 100;
s.P           = 1;
s.DoF         = 2;

CantileverAxial = FEMInputWriter(s);
CantileverAxial.createTest();

% Create elasticity problem
a.fileName = fileName;
data = FemDataContainer(a);
fem = FEM.create(data);
fem.solve();
% fem.print(fileName);

% Get main results
uTest = fem.uFun(1,1);
e_test = fem.strainFun(1,1);
sig_test = fem.stressFun(1,1);

%% DOMAIN DECOMPOSITION (N SUBDOMAINS)
subdomains = 6;
mesh = data.mesh;
material = data.material;

% Mesh decomposition
[subMesh, subBoundMesh] = MeshDecomposition(subdomains, mesh);

% Stiffness matrices computation
superMesh = subMesh(1);
[subK] = SubdomainStiffnessMatrix(superMesh,material);
[K,superK,nDofsLeft,nDofsRight,DofsCondensed,DofsBoundary] = SuperelementStiffnessMatrix(subK,subBoundMesh);

% Mass matrices computation
[subC] = SubdomainMassMatrix(subBoundMesh);

% Global LHS matrix assembly
[GlobalK,totalDofsDomain] = AssemblyGlobalStiffness(superK,subdomains);
[GlobalC] = AssemblyGlobalMass(subC);
[GlobalLHS] = AssemblyGlobalLHS(GlobalK,GlobalC);

% Global RHS matrix assembly
[GlobalRHS] = AssemblyGlobalRHS(GlobalLHS,subdomains,1,totalDofsDomain,nDofsLeft,nDofsRight);

% Solver
s.type = "DIRECT";
Sol = Solver.create(s);
GlobalU = Sol.solve(GlobalLHS,GlobalRHS);

% Decompose solution
[u,lambda] = DecomposeSolution(GlobalU,subC,subdomains,nDofsLeft,nDofsRight);

% Reconstruction
[uFull] = FullMeshReconstruction(GlobalU,K,DofsCondensed,DofsBoundary,subdomains);

% Results representation 
%Plots(u,lambda,subC,subMesh,subBoundMesh);
PrintResults(uFull,subMesh);

% Errors
[MaxError] = ConnectionErrorLinf(GlobalU,GlobalLHS,subC,totalDofsDomain)
[Error]    = ConnectionErrorL2(GlobalU,GlobalLHS,subC,totalDofsDomain)
[RelativeError,TipError] = TipErrorL2(u,subBoundMesh,uTest,mesh)



%% FUNCTIONS

function [subMesh, subBoundMesh] = MeshDecomposition(subdomains, mesh)
    subMesh = [];
    subBoundMesh = [];

    xBaricenter = mesh.computeBaricenter();
    xBaricenter = xBaricenter(1,:);

    lengthDomain = max(mesh.coord(:,1));
    lengthSubdomain = lengthDomain/subdomains;

    initial_coord = 0;
    s.coord = mesh.coord;
    figure
    for i=1:subdomains
        isSubdomain = (xBaricenter > initial_coord) & (xBaricenter <= initial_coord + lengthSubdomain);
        subdomainConnec = mesh.connec(isSubdomain,:);
        s.connec = subdomainConnec;
        subdomainMesh = Mesh(s);
        subMesh = cat(1,subMesh,subdomainMesh.computeCanonicalMesh());
        subMesh(i).plot

        boundary = subMesh(i).createBoundaryMesh();
        subBoundMesh = cat(1,subBoundMesh,[boundary{1}, boundary{2}]);
        subBoundMesh(i,1).mesh.plot
        subBoundMesh(i,2).mesh.plot
        initial_coord = initial_coord + lengthSubdomain;
    end
end

function [subK] = SubdomainStiffnessMatrix(superMesh,material)
        s.fValues = zeros(superMesh.nnodes, superMesh.ndim);
        s.mesh    = superMesh;
        p1 = P1Function(s);

        s.type     = 'ElasticStiffnessMatrix';
        s.mesh     = superMesh;
        s.fun      = p1;
        s.material = material;
        lhs = LHSintegrator.create(s);
        subK = lhs.compute();
end

function [K,superK,nDofsLeft,nDofsRight,DofsCondensed,DofsBoundary] = SuperelementStiffnessMatrix(subK,subBoundMesh)
    ndim = subBoundMesh(1,1).mesh.ndim;
    idxNodesRight = [subBoundMesh(1,2).globalConnec(:,1); subBoundMesh(1,2).globalConnec(end,2)];
    idxNodesLeft = [subBoundMesh(1,1).globalConnec(:,1); subBoundMesh(1,1).globalConnec(end,2)];
    DofsRight = zeros(length(idxNodesRight)*ndim,1);
    DofsLeft = zeros(length(idxNodesLeft)*ndim,1);
    for j=1:ndim
        DofsRight(j:2:(end-ndim+j))  = ndim*idxNodesRight+(-ndim+j);
        DofsLeft(j:2:(end-ndim+j)) = ndim*idxNodesLeft+(-ndim+j);
    end
    DofsBoundary = [DofsLeft; DofsRight];
    nDofsLeft = length(DofsLeft);
    nDofsRight = length(DofsRight);
    TotalDofs = 1:1:size(subK,1);
    DofsCondensed = setdiff(TotalDofs,DofsBoundary)';

    K.Kii = subK(DofsCondensed,DofsCondensed);
    K.Kib = subK(DofsCondensed,DofsBoundary);
    K.Kbi = subK(DofsBoundary,DofsCondensed);
    K.Kbb = subK(DofsBoundary,DofsBoundary);

    superK = K.Kbb - K.Kbi*(K.Kii\K.Kib);
end


function [subC] = SubdomainMassMatrix(subBoundMesh)
    subdomains = size(subBoundMesh,1);
    boundaries = size(subBoundMesh,2);
    subC = cell(subdomains,boundaries);
    for i=1:subdomains
        for j=1:boundaries
            % Right now, trial and test functions need to be defined here
            % inside.
            ndim  = subBoundMesh(i,j).mesh.ndim;
            test  = P1Function.create(subBoundMesh(i,j).mesh,ndim); % disp
            trial = P1Function.create(subBoundMesh(i,j).mesh,ndim); % lambda

            s.type    = 'MassMatrix';
            s.mesh    = subBoundMesh(i,j).mesh;
            s.test    = test;
            s.trial   = trial;
            lhs       = LHSintegrator.create(s);
            C         = lhs.compute();
            subC(i,j) = {C};
        end
    end
end

function [GlobalK,totalDofsDomain] = AssemblyGlobalStiffness(superK,subdomains)
   GlobalK = superK;
   for i=2:subdomains
        GlobalK = [                GlobalK                 zeros(size(GlobalK,1),size(superK,2));
                   zeros(size(superK,1),size(GlobalK,2))                superK               ];
   end
   totalDofsDomain = size(GlobalK,1);
end

function [GlobalC] = AssemblyGlobalMass(subC)
    boundaries = size(subC,1);
    % First boundary (Dirichlet)
    GlobalC = subC{1,1};
    
    % Connection boundaries
    for i=2:boundaries
        columnC = [  subC{i-1,2};
                    -subC{i,1} ];
        GlobalC = [               GlobalC                      zeros(size(GlobalC,1),size(columnC,2));
                     zeros(size(columnC,1),size(GlobalC,2))                   columnC                ];
    end

    % Last boundary (Neumann)
    columnC = [zeros(size(subC{end,2},1), size(GlobalC,2))];
    GlobalC = [ GlobalC ;
                columnC ];
end

function [GlobalLHS] = AssemblyGlobalLHS(GlobalK,GlobalC)
    GlobalLHS = [GlobalK          GlobalC       ;
                 GlobalC' zeros(size(GlobalC,2))];
end

function [GlobalRHS] = AssemblyGlobalRHS(GlobalLHS,subdomains,force,totalDofsDomain,nDofsLeft,nDofsRight)
    GlobalRHS = sparse(size(GlobalLHS,1),1);
    dim = 2; % 1 -> Axial force, 2-> Bending force
    dofsNeumann = ((subdomains-1)*(nDofsLeft+nDofsRight)+nDofsLeft+dim):2:totalDofsDomain;

    GlobalRHS(dofsNeumann) = force/length(dofsNeumann);
end

function [U,lambda] = DecomposeSolution(GlobalU,subC,subdomains,nDofsLeft,nDofsRight)
    boundaries = size(subC,2);
    U = cell(subdomains,boundaries);
    lambda = cell(1,subdomains); %Each subdomain has their corresponding Lagrange multipliers on its Left boundary.
    initialDof = 0;
    for i=1:subdomains
        for j=1:boundaries
            if j==1
                nDofsBoundary = nDofsLeft;
            elseif j==2
                nDofsBoundary = nDofsRight;
            end
            DofsBoundary = (1:1:nDofsBoundary) + initialDof;
            subU = GlobalU(DofsBoundary);
            subU = [subU(1:2:end-1), subU(2:2:end)];
            U(i,j) = {subU};
            initialDof = max(DofsBoundary);
        end
    end
    
    for i=1:subdomains
        DofsBoundLagrange = (1:1:size(subC{i,1},2)) + initialDof;
        subLambda = GlobalU(DofsBoundLagrange);
        subLambda = [subLambda(1:2:end-1), subLambda(2:2:end)];
        lambda(i) = {subLambda};
        initialDof = max(DofsBoundLagrange);
    end
end

function [MaxError] = ConnectionErrorLinf(GlobalU,GlobalLHS,subC,totalDofsDomain)
    MaxError = 0;
    totalBoundaries = size(subC,1);
    initialDof = totalDofsDomain;
    for i=1:totalBoundaries
        DofsBoundary = (1:1:+size(subC{i,1},2)) + initialDof;
        Error = full(max(abs(GlobalLHS(DofsBoundary,:)*GlobalU)));
        if Error > MaxError
            MaxError = Error;
        end
        initialDof = max(DofsBoundary);
    end
end


function [MaxError] = ConnectionErrorL2(GlobalU,GlobalLHS,subC,TotalDofsDomain)
    MaxError = 0;
    totalBoundaries = size(subC,1);
    initialDof = TotalDofsDomain;
    for i=1:totalBoundaries
        DofsBoundary = (1:1:+size(subC{i,1},2)) + initialDof;
        Error = full(sqrt(sum((GlobalLHS(DofsBoundary,:)*GlobalU).^2)));
        if Error > MaxError
            MaxError = Error;
        end
        initialDof = max(DofsBoundary);
    end
end

function Plots(U,lambda,subC,subMesh,subBoundMesh)
    % Displacement field
    boundaries = size(subC,2);
    subdomains = size(subMesh,1);
    for i=1:subdomains
        for j=1:boundaries
            s.fValues = full(U{i,j});
            s.mesh    = subBoundMesh(i,1).mesh;
            uField    = P1Function(s);
            uField.plot
        end
    end

    % Lagrage multipliers in P1
    %%%%%% TO BE DONE %%%%%%%

    % Reactions  
    for i=1:subdomains
        VectorLambda = zeros(size(lambda{i},1)*size(lambda{i},2),1);
        VectorLambda(1:2:end-1) = lambda{i}(:,1);
        VectorLambda(2:2:end)   = lambda{i}(:,2);
        React = subC{i,1}*VectorLambda;
        React = [React(1:2:end-1), React(2:2:end)];

        figure
        quiver(subBoundMesh(i).mesh.coord(:,1),subBoundMesh(i).mesh.coord(:,2),...
               React(:,1),zeros(size(React(:,1))))
        title("Reaction X (Boundary"+(i-1)+"-"+(i)+")")
        figure
        quiver(subBoundMesh(i).mesh.coord(:,1),subBoundMesh(i).mesh.coord(:,2),...
               React(:,2), zeros(size(React(:,2))))
        title("Reaction Y (Boundary"+(i-1)+"-"+(i)+")")

        TotalReact = [sum(React(:,1)) sum(React(:,2))]
    end
end

function [RelativeError,TipError] = TipErrorL2(u,subBoundMesh,uTest,mesh)
    IsGloBoundary = mesh.coord(:,1) == max(subBoundMesh(end,2).mesh.coord(:,1));

    ErrorVectorU = full(u{end,2}) - uTest.fValues(IsGloBoundary,:);
    TipError = sqrt(sum(sum(ErrorVectorU.^2)));

    NormUTest = sqrt(sum(sum(uTest.fValues(IsGloBoundary,:).^2)));
    RelativeError = TipError/NormUTest;
end

function PrintResults(u,mesh)
    subdomains = size(mesh,1);
    for i=1:subdomains
        s = [];
        s.fValues = full(u{i});
        s.mesh    = mesh(i);
        p.filename = ['domain',char(string(i))];
        ResultSubDom = P1Function(s);
        ResultSubDom.print(p);
    end
end

function [uFull] = FullMeshReconstruction(GlobalU,K,DofsCondensed,DofsBoundary,subdomains)
    uFull = cell(1,subdomains);
    for i=1:subdomains
        uB = GlobalU((i-1)*length(DofsBoundary)+1:i*length(DofsBoundary));
        uI = K.Kii\(-K.Kib*uB);
    
        u = zeros(length(uB)+length(uI),1);
        u(DofsCondensed) = uI;
        u(DofsBoundary) = uB;

        u = [u(1:2:end-1), u(2:2:end)];
        uFull(i) = {u};
    end
end
