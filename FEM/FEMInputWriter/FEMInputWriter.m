classdef FEMInputWriter < handle
    
    properties (Access = private)
        fileName
        xmax
        ymax
        zmax
        Nx
        Ny
        Nz
        P
        xmesh
        ymesh
        zmesh
        nDirichlet
        nNeumann
        mesh
        DoF
        dispPrescribed
        pointLoads
        edgeNodes
        pCase
        designVarSettings
    end
    
    methods (Access = public)
        function obj = FEMInputWriter(cParams)
            obj.init(cParams);
        end
        
        function createTest(obj)
            obj.computeMeshGrid();
            obj.createMesh();
            obj.computeCaseBoundaryConditions();
            obj.computePrescribedDisplacementsMatrix();
            obj.computePrescribedPointLoads();
            obj.computeEdgeNodes();
            obj.writeFile();
        end
    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.fileName = cParams.testName;
            obj.pCase    = cParams.problemCase;
            obj.xmax     = cParams.x1;
            obj.ymax     = cParams.y1;
            obj.zmax     = cParams.z1;
            obj.Nx       = cParams.N;
            obj.Ny       = cParams.M;
            obj.Nz       = cParams.O;
            obj.P        = cParams.P;
            obj.DoF      = cParams.DoF;
        end

        function computeMeshGrid(obj)
            x1        = linspace(0,obj.xmax,obj.Nx+1);
            x2        = linspace(0,obj.ymax,obj.Ny+1);
            if obj.DoF == 3
                x3        = linspace(0,obj.zmax,obj.Nz+1);
                [X,Y,Z]     = meshgrid(x1,x2,x3);
            elseif obj.DoF == 2
                [X,Y]     = meshgrid(x1,x2);
                Z         = zeros(size(X));
            end
            obj.xmesh = X;
            obj.ymesh = Y;
            obj.zmesh = Z;
        end
        
        function createMesh(obj)
            switch obj.DoF
                case 2
                    [F,V]    = mesh2tri(obj.xmesh,obj.ymesh,obj.zmesh,'f');
                    s.coord  = V(:,1:2);
                    s.connec = F;
                case 3
                    x = [obj.xmesh];
                    y = [obj.ymesh];
                    z = [obj.zmesh];
                    X = [x(:) y(:) z(:)];
                    Tes = delaunayn(X);
                    s.coord  = X;
                    s.connec = Tes;
            end
            obj.mesh = Mesh(s);
            %             obj.mesh.plot;
        end

        function computeCaseBoundaryConditions(obj)
            switch obj.pCase
                case 'cantilever'
                    obj.computeCantileverBoundaryConditions();
                case 'bridge'
                    obj.computeBridgeBoundaryConditions();
                case 'arch'
                    obj.computeArchBoundaryConditions();
                case 'cantilever3'
                    obj.computeCantilever3DBoundaryConditions();
            end
        end

        function computeCantilever3DBoundaryConditions(obj)
            t = 0.3*obj.ymax;
            m = obj.mesh;
            root = m.coord(:,1) == 0;
            tipLength = m.coord(:,1) == obj.xmax;
            tipLength2 = m.coord(:,2) >= t & m.coord(:,2) <= obj.ymax - t;
            tipLength3 = m.coord(:,3) >= 0 & m.coord(:,3) <= obj.zmax;
            tip            = tipLength & tipLength2 & tipLength3;
            obj.nDirichlet = find(root);
            obj.nNeumann   = find(tip);
        end

        function computeArchBoundaryConditions(obj)
            t              = 0.1*obj.ymax;
            m              = obj.mesh;
            root1          = m.coord(:,2) == 0;
            root2          = m.coord(:,1) <= t | m.coord(:,1) >= obj.xmax-t;
            root           = root1 & root2;
            tipLength      = m.coord(:,2) == 0;
            tipWidth       = m.coord(:,1) >= 0.45*obj.xmax & m.coord(:,1) <= 0.55*obj.xmax;
            tip            = tipLength & tipWidth;
            obj.nDirichlet = find(root);
            obj.nNeumann   = find(tip);
        end
        
        function computeCantileverBoundaryConditions(obj)
            t              = 0.3*obj.ymax;
            m              = obj.mesh;
            root           = m.coord(:,1) == 0;
            tipLength      = m.coord(:,1) == obj.xmax;
            tipWidth       = m.coord(:,2) < obj.ymax-t & m.coord(:,2) > t;
            tip            = tipLength & tipWidth;
            obj.nDirichlet = find(root);
            obj.nNeumann   = find(tip);
        end

        function computeBridgeBoundaryConditions(obj)
            t              = 0.1*obj.xmax;
            iCase          = 4;
            forcePos       = [iCase*10/9,(iCase+1)*10/9];
            m              = obj.mesh;
            sides          = m.coord(:,1) <= t | m.coord(:,1) >= obj.xmax - t;
            bottom         = m.coord(:,2) == 0;
            root           = sides & bottom;
            top            = m.coord(:,2) == obj.ymax;
            pos            = m.coord(:,1) > forcePos(1) & m.coord(:,1) < forcePos(2);
            tip            = top & pos;
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

            for i = 2:obj.DoF:size(Fmat,1)
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
            fileID = fopen(['Input/',obj.fileName],'w');
            obj.writeProblemData(fileID);
            obj.writeCoordinates(fileID);
            obj.writeConnectivities(fileID);
            obj.writeDirichletData(fileID);
            obj.writeNeumannData(fileID);
            obj.writeNodesSolid(fileID);
            obj.writeExternalBorderNodes(fileID);
            obj.writeExternalBorderElements(fileID);
            fclose(fileID);
        end
        
        function writeCoordinates(obj,fileID)
            COOR        = zeros(size(obj.mesh.coord,1),4);
            dims        = size(obj.mesh.coord,2);
            COOR(:,2:dims+1) = obj.mesh.coord;
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
            dims        = size(obj.mesh.connec,2);
            Tnod(:,2:dims+1) = obj.mesh.connec;
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

        function writeProblemData(obj,fileID)
            fprintf(fileID,'%%%% Data\nData_prb = {\n');
            switch obj.DoF
                case 2
                    fprintf(fileID,'''TRIANGLE'';\n''SI'';\n''2D'';\n''Plane_Stress'';\n');
                case 3
                    fprintf(fileID,'''TRIANGLE'';\n''SI'';\n''3D'';\n''Plane_Stress'';\n');
            end
            fprintf(fileID,'''ELASTIC'';\n''MACRO'';\n};\n');
        end
    end

    methods (Access = private, Static)

        function writeExternalBorderElements(fileID)
            fprintf(fileID,'%%%% External Border Elements\n%% Node\n');
            fprintf(fileID,'External_border_elements = [\n');
            fprintf(fileID,'];\n');
        end

        function bc = computeBoundaryConditionMatrix(DoF,n)
            bc = zeros(DoF*length(n),3);
            for i = 1:length(n)
                switch DoF
                    case 2
                        bc(2*i-1,1:2) = [n(i),1];
                        bc(2*i,1:2)   = [n(i),2];
                    case 3
                        bc(3*i-2,1:2) = [n(i),1];
                        bc(3*i-1,1:2) = [n(i),2];
                        bc(3*i,1:2)   = [n(i),3];
                end
            end
        end
        function writeNodesSolid(fileID)
            fprintf(fileID,'%%%% Nodes solid\n%% Node\n');
            fprintf(fileID,'nodesolid = unique(pointload_complete(:,1));\n');
        end
    end
end
