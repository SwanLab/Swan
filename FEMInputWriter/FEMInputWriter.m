classdef FEMInputWriter < handle
    
    properties (Access = private)
        fileName
        xmax
        ymax
        Nx
        Ny
        P
        xmesh
        ymesh
        zmesh
        mesh
        nDirichlet
        nNeumann
        DoF
        dispPrescribed
        pointLoads
        edgeNodes
    end
    
    methods (Access = public)
        function obj = FEMInputWriter(cParams)
            obj.init(cParams);
        end
        
        function createTest(obj)
            obj.computeMeshGrid();
            obj.createMesh();
            obj.computeBoundaryConditions();
            obj.computePrescribedDisplacementsMatrix();
            obj.computePrescribedPointLoads();
            obj.computeEdgeNodes();
            obj.writeFile();
        end
    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.fileName = cParams.testName;
            obj.xmax     = cParams.x1;
            obj.ymax     = cParams.y1;
            obj.Nx       = cParams.N;
            obj.Ny       = cParams.M;
            obj.P        = cParams.P;
            obj.DoF      = cParams.DoF;
        end
        
        function computeMeshGrid(obj)
            x1        = linspace(0,obj.xmax,obj.Nx);
            x2        = linspace(0,obj.ymax,obj.Ny);
            [X,Y]     = meshgrid(x1,x2);
            Z         = zeros(size(X));
            obj.xmesh = X;
            obj.ymesh = Y;
            obj.zmesh = Z;
        end
        
        function createMesh(obj)
            [F,V]    = mesh2tri(obj.xmesh,obj.ymesh,obj.zmesh,'f');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
            obj.mesh.plot;
        end
        
        function computeBoundaryConditions(obj)
            t              = 0.3*obj.ymax;
            m              = obj.mesh;
            root           = m.coord(:,1) == 0;
            tipLength      = m.coord(:,1) == obj.xmax;
            tipWidth       = m.coord(:,2) < obj.ymax-t & m.coord(:,2) > t;
            tip            = tipLength & tipWidth;
            obj.nDirichlet = find(root);
            obj.nNeumann   = find(tip);
        end
        
        function computePrescribedDisplacementsMatrix(obj)
            obj.dispPrescribed = obj.computeBoundaryConditionMatrix(obj.DoF,obj.nDirichlet);
        end
        
        function computePrescribedPointLoads(obj)
            Fmat  = obj.computeBoundaryConditionMatrix(obj.DoF,obj.nNeumann);
            nnode = size(Fmat,1)/obj.DoF;
            Pnod  = obj.P/nnode;
            for i = 2:2:size(Fmat,1)
                Fmat(i,3) = Pnod;
            end
            obj.pointLoads = Fmat;
        end
        
        function computeEdgeNodes(obj)
            COOR          = obj.mesh.coord;
            x             = COOR(:,1);
            y             = COOR(:,2);
            edge          = find(x==0 | x==max(x) | y==0 | y==max(y));
            obj.edgeNodes = sort(edge);
        end
        
        function writeFile(obj)
            fileID = fopen([obj.fileName],'w');
            obj.writeProblemData(fileID);
            obj.writeCoordinates(fileID);
            obj.writeConnectivities(fileID);
            obj.writeDirichletData(fileID);
            obj.writeNeumannData(fileID);
            obj.writeNodesSolid(fileID);
            obj.writeExternalBorderNodes(fileID);
            fclose(fileID);
        end
        
        function writeCoordinates(obj,fileID)
            COOR        = zeros(size(obj.mesh.coord,1),4);
            COOR(:,2:3) = obj.mesh.coord;
            COOR(:,1)   = (1:1:size(COOR,1))';
            fprintf(fileID,'%%%% Coordinates\n%% Node\n');
            fprintf(fileID,'gidcoord = [\n');
            for i=1:size(COOR,1)
                fprintf(fileID,'%d %f %f %f;\n',COOR(i,:));
            end
            fprintf(fileID,'];\n');
        end
        
        function writeConnectivities(obj,fileID)
            Tnod        = zeros(size(obj.mesh.connec,1),5);
            Tnod(:,2:4) = obj.mesh.connec;
            Tnod(:,1)   = (1:1:size(Tnod,1))';
            fprintf(fileID,'%%%% Connectivities\n%% Node\n');
            fprintf(fileID,'gidlnods = [\n');
            for i=1:size(Tnod,1)
                fprintf(fileID,'%d %d %d %d %d;\n',Tnod(i,:));
            end
            fprintf(fileID,'];\n');
        end
        
        function writeDirichletData(obj,fileID)
            fprintf(fileID,'%%%% Variable prescribed\n%% Node\n');
            fprintf(fileID,'lnodes = [\n');
            for i=1:size(obj.dispPrescribed,1)
                fprintf(fileID,'%d %d %f;\n',obj.dispPrescribed(i,:));
            end
            fprintf(fileID,'];\n');
        end
        
        function writeNeumannData(obj,fileID)
            fprintf(fileID,'%%%% Point loads\n%% Node\n');
            fprintf(fileID,'pointload_complete = [\n');
            for i=1:size(obj.pointLoads,1)
                fprintf(fileID,'%d %d %f;\n',obj.pointLoads(i,:));
            end
            fprintf(fileID,'];\n');
        end
        
        function writeExternalBorderNodes(obj,fileID)
            fprintf(fileID,'%%%% External Border Nodes\n%% Node\n');
            fprintf(fileID,'External_border_nodes = [\n');
            for i=1:size(obj.edgeNodes,1)
                fprintf(fileID,'%d;\n',obj.edgeNodes(i));
            end
            fprintf(fileID,'];\n');
        end
    end
    
    methods (Access = private, Static)
        function bc = computeBoundaryConditionMatrix(DoF,n)
            bc = zeros(DoF*length(n),3);
            for i = 1:length(n)
                bc(2*i-1,1:2) = [n(i),1];
                bc(2*i,1:2)   = [n(i),2];
            end
        end
        function writeProblemData(fileID)
            fprintf(fileID,'%%%% Data\nData_prb = {\n');
            fprintf(fileID,'''TRIANGLE'';\n''SI'';\n''2D'';\n''Plane_Stress'';\n');
            fprintf(fileID,'''ELASTIC'';\n''MACRO'';\n};\n');
        end
        function writeNodesSolid(fileID)
            fprintf(fileID,'%%%% Nodes solid\n%% Node\n');
            fprintf(fileID,'nodesolid = unique(pointload_complete(:,1));\n');
        end
    end
end
