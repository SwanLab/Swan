clc,clear,close all 

Meshes = [7 15 31 61 121;
          3 7 15 31  61 ];
RelativeTipError    = zeros(1,1,length(Meshes));
TipError            = zeros(1,1,length(Meshes));
RelativeStressError = zeros(2,2,length(Meshes));
StressError         = zeros(2,2,length(Meshes));
numberElements = [];
for k=1:length(Meshes)
%% CONVENTIONAL SOLVER
% Create geometry

fileName = 'CantileverAxialLoad';
s.testName    = [fileName,'.m'];
s.problemCase = 'cantilever';
s.x1          = 1;
s.y1          = 0.5;
s.N           = Meshes(1,k); %x dimension
s.M           = Meshes(2,k); %y dimension
s.P           = -1;
s.DoF         = 2;

CantileverAxial = FEMInputWriter(s);
CantileverAxial.createTest();

% Create elasticity problem
a.fileName = fileName;
data = FemDataContainer(a);
% fem = FEM.create(data);
% fem.solve();
% fem.print(fileName);

% Get main results
% uTest = fem.uFun(1,1);
% e_test = fem.strainFun(1,1);
% sig_test = fem.stressFun(1,1);

%% Reactions plots fem

% K = fem.LHS;
% u = zeros(length(uTest.fValues)*2,1);
% u(1:2:end-1) = uTest.fValues(:,1);
% u(2:2:end) = uTest.fValues(:,2);
% Reactions = K(fem.boundaryConditions.dirichlet,:)*u;
% Reactions = [Reactions(1:2:end-1) Reactions(2:2:end)];
% idx = data.mesh.coord(:,1) == 0;
% 
% X = data.mesh.coord(idx,1);
% Y = data.mesh.coord(idx,2);

% i=5;
% figure
% plot(Reactions(:,1),Y,'.b')
% hold on
% plot([0 0],[0 0.5],'-')
% quiver(X(5:i:end-5),Y(5:i:end-5),Reactions(5:i:end-5,1),zeros(size(Reactions(5:i:end-5,1))),"off",'b')
% title("Reaction X")
% xlabel('Force magnitude')
% ylabel('Section position')
% 
% figure
% plot(Reactions(:,2),Y,'.b')
% hold on
% plot([0 0],[0 0.5],'-')
% quiver(X(5:i:end-5),Y(5:i:end-5),Reactions(5:i:end-5,2),zeros(size(Reactions(5:i:end-5,2))),"off",'b')
% title("Reaction Y")
% xlabel('Force magnitude')
% ylabel('Section position')

% RootMoment = sum(Reactions(:,1)'*(Y-0.25));
% TestTotalReact = [sum(Reactions(:,1)) sum(Reactions(:,2)) RootMoment]
% 
% 
% idxTip = find(data.mesh.coord(:,1) == 1 & data.mesh.coord(:,2) == 0.25);
% uTipTest = uTest.fValues(idxTip,2);
% 
% Stress = [2*Reactions(1,:)/(s.y1/(s.M-1)); Reactions(2:end-1,:)/(s.y1/(s.M-1)); 2*Reactions(end,:)/(s.y1/(s.M-1))];
% plot(Stress,Y)

refTestValues = load('TestReferenceData.mat');

%% DOMAIN DECOMPOSITION (N SUBDOMAINS)

in.subdomains = 2;
in.mesh       = data.mesh;
in.mat        = data.material;
in.type       = "Two";

% Mesh decomposition
[meshDecomposed] = MeshDecomposer(in);
plotMesh(meshDecomposed);
subMesh = meshDecomposed.subMesh;
subBoundMesh = meshDecomposed.subBoundMesh;

numberElements(end+1)=data.mesh.nelem;

% Stiffness matrices computation
[subK] = SubdomainStiffnessMatrix(subMesh,in.mat);

% Mass matrices computation
[subC] = SubdomainMassMatrix(subBoundMesh);

% Global LHS matrix assembly
[GlobalK] = AssemblyGlobalStiffness(subK);
[GlobalC, TotalDofsDomain] = AssemblyGlobalMass(subC,subBoundMesh, subK, GlobalK);
[GlobalLHS] = AssemblyGlobalLHS(GlobalK,GlobalC);

% Global RHS matrix assembly
[GlobalRHS] = AssemblyGlobalRHS(GlobalLHS,subK,subBoundMesh,s.P,TotalDofsDomain);

% Solver
s.type = "DIRECT";
Sol = Solver.create(s);
GlobalU = Sol.solve(GlobalLHS,GlobalRHS);

% Decompose solution
[u,lambda] = DecomposeSolution(GlobalU,subK,subC,TotalDofsDomain);
[TotalReact,React,Stress] = ComputeReactionsAndStresses(lambda,subC,subBoundMesh);
if k==1
plotStresses(subBoundMesh,Stress,React);
end

[ErrorMax] = ConnectionErrorLinf(GlobalU,GlobalLHS,subC,TotalDofsDomain);
[ErrorL2]    = ConnectionErrorL2(GlobalU,GlobalLHS,subC,TotalDofsDomain);

[RelativeTipError(:,:,k),TipError(:,:,k)] = TipErrorL2(u,subMesh,refTestValues);
[RelativeStressError(:,:,k),StressError(:,:,k)] = StressErrorL2(Stress,subBoundMesh,refTestValues);
if k==1
PrintResults(u,subMesh);
end
end

PlotErrors(RelativeTipError,RelativeStressError,numberElements)

%% FUNCTIONS

function PlotErrors(RelativeTipError,RelativeStressError,numberElements)
    figure
    loglog(numberElements,squeeze(RelativeTipError),'-o')
    grid on
    title('Relative error vs. elements')
    xlabel('Number of elements')
    ylabel('Relative error')
    ytickformat('%.2f')
    legend('Tip displacement (x=1,y=0.25)')

    figure
    loglog(numberElements,squeeze(RelativeStressError(1,1,:)),'-o')
    hold on
    loglog(numberElements,squeeze(RelativeStressError(1,2,:)),'-o')
    loglog(numberElements,squeeze(RelativeStressError(2,1,:)),'-o')
    loglog(numberElements,squeeze(RelativeStressError(2,2,:)),'-o')
    grid on
    title('Relative error vs. elements')
    xlabel('Number of elements')
    ylabel('Relative error')
    ytickformat('%.3f')
    legend('Normal stress (x=0,y=0.05)', ...
            'Normal stress (x=0,y=0.1)', ...
            'Shear stress (x=0,y=0.05)', ...
            'Shear stress (x=0,y=0.1)')
end



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
            trial = P1Function.create(subBoundMesh(i,j).mesh,ndim); % lambda %%%%%%%%%%%%%%%%%%%%%%%%%

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

function [GlobalK] = AssemblyGlobalStiffness(subK)
   subdomains = size(subK,1);
   GlobalK = subK{1};
   for i=2:subdomains
        GlobalK = [                GlobalK                 zeros(size(GlobalK,1),size(subK{i},2));
                   zeros(size(subK{i},1),size(GlobalK,2))                subK{i}                ];
   end
end

function [GlobalC, TotalDofsDomain] = AssemblyGlobalMass(subC,subBoundMesh,subK,GlobalK)
    boundaries = size(subC,1);
    DofsDispGlobal = size(GlobalK,1);
    TotalDofsDomain = size(subK{1},1);
    ndim = subBoundMesh(1).mesh.ndim;

    DofsLagrSub = size(subC{1,1},2);
    columnC = sparse(DofsDispGlobal,DofsLagrSub);

    idxNodesLeft = [subBoundMesh(1,1).globalConnec(:,1); subBoundMesh(1,1).globalConnec(end,2)];
    DofsLeft = zeros(length(idxNodesLeft)*ndim,1);
    for i=1:ndim
        DofsLeft(i:2:(end-ndim+i)) = ndim*idxNodesLeft+(-ndim+i);
    end

    columnC(DofsLeft,:) = subC{1,1};
    GlobalC = columnC;

    for i=2:boundaries
        DofsLagrSub = size(subC{i,1},2);
        columnC = sparse(DofsDispGlobal,DofsLagrSub);

        idxNodesRight = [subBoundMesh(i-1,2).globalConnec(:,1); subBoundMesh(i-1,2).globalConnec(end,2)];
        idxNodesLeft = [subBoundMesh(i,1).globalConnec(:,1); subBoundMesh(i,1).globalConnec(end,2)];
        DofsRight = zeros(length(idxNodesRight)*ndim,1);
        DofsLeft = zeros(length(idxNodesLeft)*ndim,1);
        for j=1:ndim
            DofsRight(j:2:(end-ndim+j))  = ndim*idxNodesRight+(-ndim+j) + TotalDofsDomain - size(subK{i-1},1);
            DofsLeft(j:2:(end-ndim+j)) = ndim*idxNodesLeft+(-ndim+j) + TotalDofsDomain;
        end
        
        columnC(DofsRight,:) = columnC(DofsRight,:) + subC{i-1,2};
        columnC(DofsLeft,:) = columnC(DofsLeft,:) - subC{i,1};
        GlobalC = [GlobalC columnC];

        TotalDofsDomain = TotalDofsDomain + size(subK{i},1);
    end
end

function [GlobalLHS] = AssemblyGlobalLHS(GlobalK,GlobalC)
    GlobalLHS = [GlobalK          GlobalC       ;
                 GlobalC' zeros(size(GlobalC,2))];
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


function  [totalReact,React,Stress] = ComputeReactionsAndStresses(lambda,subC,subBoundMesh)
     boundaries = size(subBoundMesh,1);
     Stress = cell(size(boundaries));
    % Stresses
    for i=1:boundaries
        s.fValues = full(lambda{i});
        s.mesh = subBoundMesh(i,1).mesh;
        if size(subC{1,1},1)==size(subC{1,1},2)
            Stress(i) = {P1Function(s)};
        else
            StressFun = P0Function(s);
            Stress(i) = {StressFun.project('P1')};
        end

    end

    % Reactions  
    totalBoundaries = size(subC,1);
    totalReact = zeros(totalBoundaries,3);
    for i=1:totalBoundaries
        VectorLambda = zeros(size(lambda{i},1)*size(lambda{i},2),1);
        VectorLambda(1:2:end-1) = lambda{i}(:,1);
        VectorLambda(2:2:end)   = lambda{i}(:,2);
        if i==1
            React = -subC{i,1}*VectorLambda;
        else
            React = subC{i,1}*VectorLambda;
        end
        React = [React(1:2:end-1), React(2:2:end)];

        X = subBoundMesh(i).mesh.coord(1,1);
        Y = subBoundMesh(i).mesh.coord(:,2);
        
%         figure
%         plot(React(:,1),Y,'.b')
%         hold on
%         plot([0 0],[0 0.5],'-')
%         quiver(zeros(size(Y(2:2:end-2))),Y(2:2:end-2),React(2:2:end-2,1),zeros(size(React(2:2:end-2,1))),"off",'b')
%         title("Reaction X (Boundary "+(i-1)+"-"+(i)+")"+" [X="+X+"]")
%         xlabel('Force magnitude')
%         ylabel('Section position')
% 
%         figure
%         plot(React(:,2),Y,'.b')
%         hold on
%         plot([0 0],[0 0.5],'-')
%         quiver(zeros(size(Y(2:2:end-2))),Y(2:2:end-2),React(2:2:end-2,2),zeros(size(React(2:2:end-2,2))),"off",'b')
%         title("Reaction Y (Boundary "+(i-1)+"-"+(i)+")"+" [X="+X+"]")
%         xlabel('Force magnitude')
%         ylabel('Section position')
      
        Center = subBoundMesh(i).mesh.coord(end,2)/2;
        RootMoment = sum(React(:,1)'*(subBoundMesh(i).mesh.coord(:,2) - Center));
        totalReact(i,:) = [sum(React(:,1)) sum(React(:,2)) RootMoment];

    end
end

function [RelativeTipError,TipError] = TipErrorL2(u,subMesh,refTestData)
    uTipTest = refTestData.refU;

%     IsSubBoundary = subMesh(end).coord(:,1) == max(subMesh(end).coord(:,1));
%     IsGloBoundary = mesh.coord(:,1) == max(subMesh(end).coord(:,1));
% 
%     ErrorVectorU = full(u{end}(IsSubBoundary,:)) - uTest.fValues(IsGloBoundary,:);
%     FullTipError = sqrt(sum(sum(ErrorVectorU.^2)));
% 
%     NormUTest = sqrt(sum(sum(uTest.fValues(IsGloBoundary,:).^2)));
%     FullTipRelativeError = FullTipError/NormUTest;

    idxTip = find(subMesh(end).coord(:,1) == 1 & subMesh(end).coord(:,2) == 0.25);
    uTip   = full(u{end}(idxTip,2));

    TipError = sqrt((uTipTest-uTip).^2);
    RelativeTipError = TipError/abs(uTipTest);
end


function [RelativeStressError,StressError] = StressErrorL2(Stress,subBoundMesh, refTestData)
    refSigma = refTestData.refSigma;
    refTau = refTestData.refTau;

    idx = subBoundMesh(1,1).mesh.coord(:,1) == 0;
    Y = subBoundMesh(1,1).mesh.coord(idx,2);
    
    Sigma = [interp1(Y,Stress{1}.fValues(:,1),0.05) interp1(Y,Stress{1}.fValues(:,1),0.1)];
    Tau   = [interp1(Y,Stress{1}.fValues(:,2),0.05) interp1(Y,Stress{1}.fValues(:,2),0.1)];
   
    StressError = [sqrt((Sigma-refSigma).^2); sqrt((Tau-refTau).^2)];
    RelativeStressError = [abs(StressError(1,:)./refSigma);  abs(StressError(2,:)./refTau)];
end

function PrintResults(u,subMesh)
    subdomains = size(subMesh,1);
    for i=1:subdomains
        subMesh(i).coord(:,1) = subMesh(i).coord(:,1) + (i-1)*0.001;
        s = [];
        s.fValues = full(u{i});
        s.fValues = [s.fValues zeros(size(s.fValues,1),1)]; % In order to plot deformed (z=0)
        s.mesh    = subMesh(i);
        p.filename = ['domain',char(string(i))];
        ResultSubDom = P1Function(s);
        ResultSubDom.print(p);
    end
end


