classdef PerimeterVolumeProblem < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        fileName
        inputFile
        topOptSet
        topOptProblem
        mesh
        nx
        ny
        type
    end
    
    methods (Access = public)
        
        function obj = PerimeterVolumeProblem()
            obj.init();
            obj.createMesh();
            obj.createInputFile();
            obj.createSettings();
            obj.solve();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName = 'PerimeterVolumeTest';
            obj.inputFile = 'SquareMacroStrucutred';
            obj.nx = 100;
            obj.ny = 100;
            obj.type = 'QUAD';%'QUAD'%'TRIANGLE';%
        end
        
        function createInputFile(obj)
            s.fileName        = fullfile('Input',obj.inputFile,[obj.inputFile,'.m']);
            s.coord           = obj.mesh.coord;
            s.connec          = obj.mesh.connec;
            s.masterSlave     = [];
            s.isElemInThisSet = [];
            s.corners         = [];
            s.scale           = 'MACRO';
            s.pdim            = '2D';
            s.ptype           = 'ELASTIC';
            s.resultsDir      = fullfile('Input',obj.inputFile);
            s.type            = 'TRIANGLE';
            iPrinter = InputFemFilePrinter(s);
            iPrinter.print();
        end
        
        function createSettings(obj)
           s = SettingsTopOptProblem(obj.fileName);
           s = obj.setInitialAndFinalEpsilon(s);
           obj.topOptSet = s;
        end
        
        function s = setInitialAndFinalEpsilon(obj,s)
           sT = s.incrementalSchemeSettings.targetParamsSettings;
           m  = obj.mesh;
           sT.epsilonPerFinal   = 2*m.computeMeanCellSize();
           sT.epsilonPerInitial = 124*m.computeMeanCellSize();
           %sT.epsilonPerInitial = m.computeCharacteristicLength();
           s.incrementalSchemeSettings.targetParamsSettings = sT;
        end
        
        function createMesh(obj)
            switch obj.type
                case 'TRIANGLE'
                    obj.createTriangularMesh();
                case 'QUAD'
                    obj.createQuadMesh();
            end
        end
        
        function createTriangularMesh(obj)
            xN = linspace(0,1,obj.nx);
            yN = linspace(0,1,obj.ny);
            x = repmat(xN',1,obj.ny);
            y = repmat(yN,obj.nx,1);
            s.coord = [x(:) y(:)];
            s.connec = delaunay(s.coord(:,1),s.coord(:,2));
            obj.mesh = Mesh.create(s);
        end
        
        function createQuadMesh(obj)
            s.xnumintervals = obj.nx;
            s.ynumintervals = obj.ny;
            s.xrange = [0,1];
            s.yrange = [0,1];
            g = rectgrid(s);
            sM.coord(:,1) = g.X;
            sM.coord(:,2) = g.Y;
            sM.connec = g.VI;
            obj.mesh = Mesh.create(sM);
        end

        function solve(obj)
            obj.topOptProblem = TopOpt_Problem(obj.topOptSet);
            obj.topOptProblem.computeVariables();
            obj.topOptProblem.postProcess();
        end
        
    end
    
end