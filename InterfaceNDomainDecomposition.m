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
in.subdomains = 2;
in.mesh       = data.mesh;
in.mat        = data.material;
in.type       = "Three";

% Mesh decomposition
[meshDecomposed] = MeshDecomposition(in);
plotMesh(meshDecomposed);

% Stiffness matrices computation
[subK] = SubdomainStiffnessMatrix(subMesh,material);

% Mass matrices computation
[subC] = SubdomainMassMatrix(subBoundMesh);

% Interface matrices computation
[subD] = SubdomainInterfaceMatrix(subBoundMesh,interMesh);

% Global LHS matrix assembly
[GlobalK] = AssemblyGlobalStiffness(subK);
[GlobalC, TotalDofsDomain] = AssemblyGlobalBoundary(subC,subBoundMesh, subK, GlobalK);
[GlobalD] = AssemblyGlobalInterface(subD);
[GlobalLHS] = AssemblyGlobalLHS(GlobalK,GlobalC,GlobalD);

% Global RHS matrix assembly
[GlobalRHS] = AssemblyGlobalRHS(GlobalLHS,subK,subBoundMesh,1,TotalDofsDomain);

% Solver
s.type = "DIRECT";
Sol = Solver.create(s);
GlobalU = Sol.solve(GlobalLHS,GlobalRHS);

% Decompose solution
[u,lambda] = DecomposeSolution(GlobalU,subK,subC,TotalDofsDomain);
Plots(u,lambda,subC,subMesh,subBoundMesh);
[MaxError] = ConnectionErrorLinf(GlobalU,GlobalLHS,subC,TotalDofsDomain)
[Error]    = ConnectionErrorL2(GlobalU,GlobalLHS,subC,TotalDofsDomain)
[RelativeError,TipError] = TipErrorL2(u,subMesh,uTest,mesh)

PrintResults(u,subMesh);

%% FUNCTIONS

function [subK] = SubdomainStiffnessMatrix(subMesh,material)
    subdomains = length(subMesh);
    subK = cell(subdomains,1);
    for i=1:subdomains
        s.fValues = zeros(subMesh(i).nnodes, subMesh(i).ndim);
        s.mesh    = subMesh(i);
        p1 = P1Function(s);

        s.type     = 'ElasticStiffnessMatrix';
        s.mesh     = subMesh(i);
        s.fun      = p1;
        s.material = material;
        lhs = LHSintegrator.create(s);
        K = lhs.compute();
        subK(i) = {K};
    end
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
            trial = P0Function.create(subBoundMesh(i,j).mesh,ndim); % lambda

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

function [subD] = SubdomainInterfaceMatrix(subBoundMesh,interMesh)
    interfaces = size(interMesh,1);
    interfaceBoundaries = 2; % Is there a better way to get this number?
    subD = cell(interfaces,interfaceBoundaries);
    for i=2:interfaces % There is no interface at Dirichlet boundary.
        for j=1:interfaceBoundaries 
            ndim  = subBoundMesh(i,j).mesh.ndim;
            test  = P0Function.create(subBoundMesh(i+j-2,3-j).mesh,ndim); % lambda (must be the same as the ones for C). Next revision, this sohuld be an input.
            trial = P0Function.create(interMesh(i).mesh,ndim); % interface (ndim should be the one of the Mesh. Correct)
            
            s.type    = 'MassMatrix';
            s.mesh    = interMesh(i).mesh;
            s.test    = test;
            s.trial   = trial;
            lhs       = LHSintegrator.create(s);
            D         = lhs.compute();
            subD(i,j) = {D};
        end
    end
end


function [GlobalK] = AssemblyGlobalStiffness(subK)
   subdomains = size(subK,1);
   GlobalK = subK{1};
   for i=2:subdomains
        GlobalK = [                GlobalK                 zeros(size(GlobalK,1),size(subK{i},2));
                   zeros(size(subK{i},1),size(GlobalK,2))                subK{i}                ];
   end
end

function [GlobalC, TotalDofsDomain] = AssemblyGlobalBoundary(subC,subBoundMesh,subK,GlobalK)
    boundaries = size(subC,1);
    interfaceBoundaries = 2; % Is there a better way to get this number?
    DofsDispGlobal = size(GlobalK,1);
    TotalDofsDomain = size(subK{1},1);
    ndim = subBoundMesh(1).mesh.ndim;

    DofsLagrSub = size(subC{1,1},2);
    columnC = sparse(DofsDispGlobal,DofsLagrSub);

    idxNodes = [subBoundMesh(1,1).globalConnec(:,1); subBoundMesh(1,1).globalConnec(end,2)];
    DofsBoundary = zeros(length(idxNodes)*ndim,1);
    for i=1:ndim
        DofsBoundary(i:2:(end-ndim+i)) = ndim*idxNodes+(-ndim+i);
    end

    columnC(DofsBoundary,:) = subC{1,1};
    GlobalC = columnC;

    for i=2:boundaries
        for j=1:interfaceBoundaries
            DofsLagrSub = size(subC{i,1},2);
            columnC = sparse(DofsDispGlobal,DofsLagrSub);
            idxNodes = [subBoundMesh(i+j-2,3-j).globalConnec(:,1); subBoundMesh(i+j-2,3-j).globalConnec(end,2)];
            DofsBoundary = zeros(length(idxNodes)*ndim,1);
            for j=1:ndim
                DofsBoundary(j:2:(end-ndim+j))  = ndim*idxNodes+(-ndim+j) + TotalDofsDomain - size(subK{i-1},1);
            end
            columnC(DofsBoundary,:) = -(3-2*j)*subC{i+j-2,3-j}; % If j=1 -> subC / j=2 -> -subC
            GlobalC = [GlobalC columnC];

            TotalDofsDomain = TotalDofsDomain + size(subK{i},1);
        end
        TotalDofsDomain = TotalDofsDomain - size(subK{i},1);
    end
end

function [GlobalLHS] = AssemblyGlobalLHS(GlobalK,GlobalC,GlobalD)
    GlobalLHS = [                 GlobalK                        GlobalC            zeros(size(GlobalC,1),size(GlobalD,2));
                                  GlobalC'                zeros(size(GlobalC,2))                   GlobalD                ;
                 zeros(size(GlobalD,2),size(GlobalC,1))          GlobalD'                    zeros(size(GlobalD,2))       ];
end

function [GlobalD] = AssemblyGlobalInterface(subD)
    interfaces = size(subD,1);
    GlobalD = zeros(size(subD{2,2}));
    for i=2:interfaces
        GlobalD = [     GlobalD                              zeros(size(GlobalD,1),size(subD{i,1},2));
                   zeros(size(subD{i,1},1),size(GlobalD,2))                   -subD{i,1}             ;
                   zeros(size(subD{i,2},1),size(GlobalD,2))                    subD{i,2}             ];
    end
end

function [GlobalRHS] = AssemblyGlobalRHS(GlobalLHS,subK,subBoundMesh,Force, TotalDofsK)
    GlobalRHS = sparse(size(GlobalLHS,1),1);
    ndim = subBoundMesh(end).mesh.ndim;
    dim = 2; % 1 -> Axial force, 2-> Bending force

    idxNodesRight = [subBoundMesh(end,2).globalConnec(:,1); subBoundMesh(end,2).globalConnec(end,2)];
    DofsRight = ndim*idxNodesRight+(-ndim+dim) + TotalDofsK - size(subK{end},1);

    GlobalRHS(DofsRight) = Force/length(DofsRight);
end

function [U,lambda] = DecomposeSolution(GlobalU,subK,subC,totalDofsDomain)
    subdomains = size(subK,1);
    totalBoundaries = size(subC,1);
    U = cell(1,subdomains);
    lambda = cell(1,totalBoundaries);
    initialDof = 0;
    for i=1:subdomains
        DofsSubdomain = (1:1:size(subK{i},1)) + initialDof;
        subU = GlobalU(DofsSubdomain);
        subU = [subU(1:2:end-1), subU(2:2:end)];
        U(i) = {subU};
        initialDof = max(DofsSubdomain);
    end

    initialDof = totalDofsDomain;
    for i=1:totalBoundaries
        DofsBoundary = (1:1:+size(subC{i,1},2)) + initialDof;
        subLambda = GlobalU(DofsBoundary);
        subLambda = [subLambda(1:2:end-1), subLambda(2:2:end)];
        lambda(i) = {subLambda};
        initialDof = max(DofsBoundary);
    end
end

function [MaxError] = ConnectionErrorLinf(GlobalU,GlobalLHS,subC,TotalDofsDomain)
    MaxError = 0;
    totalBoundaries = size(subC,1);
    initialDof = TotalDofsDomain;
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
    subdomains = size(subMesh,1);
    for i=1:subdomains
        s.fValues = full(U{i});
        s.mesh    = subMesh(i);
        uField = P1Function(s);
        uField.plot
    end

    % Lagrage multipliers in P1
    %%%%%% TO BE DONE %%%%%%%

    % Reactions  
    totalBoundaries = size(subC,1);
    for i=1:totalBoundaries
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

function [RelativeError,TipError] = TipErrorL2(u,subMesh,uTest,mesh)
    IsSubBoundary = subMesh(end).coord(:,1) == max(subMesh(end).coord(:,1));
    IsGloBoundary = mesh.coord(:,1) == max(subMesh(end).coord(:,1));

    ErrorVectorU = full(u{end}(IsSubBoundary,:)) - uTest.fValues(IsGloBoundary,:);
    TipError = sqrt(sum(sum(ErrorVectorU.^2)));

    NormUTest = sqrt(sum(sum(uTest.fValues(IsGloBoundary,:).^2)));
    RelativeError = TipError/NormUTest;
end

function PrintResults(u,subMesh)
    subdomains = size(subMesh,1);
    for i=1:subdomains
        s = [];
        s.fValues = full(u{i});
        s.mesh    = subMesh(i);
        p.filename = ['domain',char(string(i))];
        ResultSubDom = P1Function(s);
        ResultSubDom.print(p);
    end
end

