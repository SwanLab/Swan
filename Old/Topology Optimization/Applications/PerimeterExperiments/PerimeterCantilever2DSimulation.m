classdef PerimeterCantilever2DSimulation < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        fileName
        topSet
        epsilonMin
        epsilonMax
        levelSet
        backgroundMesh
        initialVolume
        lambda
        nSteps
        maxIters
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = PerimeterCantilever2DSimulation()
            obj.init();            
            obj.createSettings();
            obj.setEpsilonValue();
            obj.setLevelSetValue();
            obj.solve();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName = 'CantileverTriangleFinePerimeter';
            obj.lambda = 0;
            obj.createBackgroundMesh();
            obj.initLevelSet();
            obj.createEpsilons();
        end
        
        function createBackgroundMesh(obj)
           s = SettingsTopOptProblem(obj.fileName);               
           m = s.designVarSettings.femData.mesh;           
           obj.backgroundMesh = m;
        end
        
        function initLevelSet(obj)           
           d = load('LevelSetTri.mat');    
           obj.levelSet = d.x;
        end
        
        function createEpsilons(obj)
           m = obj.backgroundMesh;            
           obj.epsilonMax = min(32*m.computeMeanCellSize(),m.computeCharacteristicLength());
           obj.epsilonMin = 1*m.computeMeanCellSize();
        end
        
        function createSettings(obj)
           obj.topSet = SettingsTopOptProblem(obj.fileName);     
        end        

        function setEpsilonValue(obj)
             t = obj.topSet.incrementalSchemeSettings.targetParamsSettings;
             t.epsilonPerFinal   = obj.epsilonMin;
             t.epsilonPerInitial = obj.epsilonMax;
             obj.topSet.incrementalSchemeSettings.targetParamsSettings = t;
        end
        
        function setLevelSetValue(obj)
            t = obj.topSet;
            t.designVarSettings.initialCase = 'given';
            t.designVarSettings.creatorSettings.value = obj.levelSet;   
            obj.topSet = t;            
        end
               
        function setNsteps(obj,iter)
           t = obj.topSet;            
           nSt = obj.nSteps(iter);
           t.incrementalSchemeSettings.nSteps = nSt;
           t.incrementalSchemeSettings.targetParamsSettings.nSteps = nSt;           
           obj.topSet = t;              
        end
           
        function solve(obj)
            optProblem = TopOpt_Problem(obj.topSet);
            optProblem.dualVariable.value = obj.lambda;          
            optProblem.computeVariables(); 
            obj.lambda = optProblem.dualVariable.value;
            obj.levelSet = optProblem.designVariable.value;
            optProblem.constraint.shapeFunctions{1}.computeFunction();
            v = optProblem.constraint.shapeFunctions{1}.value;            
            obj.initialVolume = v ;
        end
    end
    
end