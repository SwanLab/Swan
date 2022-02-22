classdef StressMinimizationWithGivenTopology < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        dataRes
        designVariable
    end
    
    properties (Access = private)
       microCase
       path
       fileResName
       initialIter
       vademecumName
       topOptProblem
    end
    
    methods (Access = public)
        
        function obj = StressMinimizationWithGivenTopology()
            obj.init();
            obj.wrapResAndMesh();
            obj.createTopOpt();
            obj.createDesignVariable();
            obj.solveTopOpt();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.microCase = 'Rectangle';
            obj.vademecumName = 'SuperEllipseQMax';
            meshType  = 'Medium';
            fCase = [obj.microCase,meshType];
            obj.path = ['/media/alex/My Passport/LatticeResults/StressNorm',fCase,'/'];
            obj.fileResName  = 'ExperimentingPlot';
            obj.initialIter = 3;
        end
        
        function createTopOpt(obj)
            s = SettingsTopOptProblem('ExperimentingPlotRectangle');
            s.designVarSettings.creatorSettings.m1 = obj.dataRes.DesignVar1;
            s.designVarSettings.creatorSettings.m2 = obj.dataRes.DesignVar2;
            s.designVarSettings.creatorSettings.alpha0 =  obj.dataRes.AlphaGauss';
            obj.topOptProblem = TopOpt_Problem(s);
        end        
        
        function solveTopOpt(obj)
            obj.topOptProblem.computeVariables();
            obj.topOptProblem.postProcess();
        end
        
        function wrapResAndMesh(obj)
            s.folderPath = obj.path;
            s.fileName   = [obj.fileResName,num2str(obj.initialIter)];
            w = WrapperMshResFiles(s);
            w.compute();
            obj.mesh    = w.mesh;
            obj.dataRes = w.dataRes;
        end
        
        function createDesignVariable(obj)
            s = SettingsDesignVariable();
            s.mesh                  = Mesh_Total(obj.mesh);
            s.type                  = 'MicroParams';
            s.initialCase           = 'full';
            s.scalarProductSettings = obj.createScalarProductParams();
            s.creatorSettings       = obj.createCreatorSettings();
            desVar = DesignVariable.create(s);
            desVar.alpha = obj.dataRes.AlphaGauss';
            obj.designVariable = desVar;
        end
        
        function s = createScalarProductParams(obj)
            s.epsilon = [];
            s.mesh = obj.mesh;
        end
        
        function s = createCreatorSettings(obj)
            s.m1 = obj.dataRes.DesignVar1;
            s.m2 = obj.dataRes.DesignVar2;
            s.homogSettings.fileName = obj.vademecumName;
            s.homogSettings.type     = 'ByVademecum';
        end                
        
    end
    
end