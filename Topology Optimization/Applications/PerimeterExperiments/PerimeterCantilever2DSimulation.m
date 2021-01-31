classdef PerimeterCantilever2DSimulation < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        fileName
        topSet
        epsilons
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
            for iter = 1:length(obj.epsilons)
                obj.createSettings();
                obj.setEpsilonValue(iter);
                obj.setLevelSetValue();
                obj.setInitialVolume();
                obj.setNsteps(iter);
                obj.setMaxIter(iter);
                obj.setOptimizer(iter);
                obj.solve();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName = 'CantileverTriangleFinePerimeter';
            obj.lambda = 0;
            obj.createBackgroundMesh();
            obj.initLevelSet();
            obj.createEpsilons();
            obj.createInitialVolume();
            obj.createNsteps();
            obj.createMaxIters();
        end
        
        function createBackgroundMesh(obj)
           s = SettingsTopOptProblem(obj.fileName);               
           m = s.designVarSettings.femData.mesh;           
           obj.backgroundMesh = m;
        end
        
        function initLevelSet(obj)
           m = obj.backgroundMesh;        
         %  obj.levelSet = -ones(size(m.coord,1),1);
         %  d = load('LevelSetQuad.mat');        
            d = load('LevelSetTri.mat');    
           obj.levelSet = d.x;
        end
        
        function createEpsilons(obj)
           m = obj.backgroundMesh;            
          epsMax = 1000*m.computeCharacteristicLength();
          epsMin = 4*m.computeMeanCellSize();
           pow = 6;
           nStep = ceil(log(epsMax/epsMin)/log(pow))+1;
           nStep = 6;
          obj.epsilons = sort(epsMin*pow.^(0:nStep-1),'descend'); 
          obj.epsilons = epsMin*10.^[6:-1:1];
        end
        
        function createSettings(obj)
           obj.topSet = SettingsTopOptProblem(obj.fileName);     
        end
        
        function createInitialVolume(obj)
           s = SettingsTopOptProblem(obj.fileName);               
           vI = s.incrementalSchemeSettings.targetParamsSettings.VfracInitial;
           obj.initialVolume = vI;            
        end
        
        function createNsteps(obj)
            s = SettingsTopOptProblem(obj.fileName);                           
    
            obj.nSteps = ones(length(obj.epsilons),1);
            obj.nSteps(1) = s.incrementalSchemeSettings.nSteps;            
        end
        
        function createMaxIters(obj)
            s = SettingsTopOptProblem(obj.fileName);                           
            obj.maxIters    = 10*ones(length(obj.epsilons),1);
            obj.maxIters(1) = s.optimizerSettings.maxIter;                    
        end
        
        function setEpsilonValue(obj,iter)
%             t = obj.topSet.incrementalSchemeSettings.targetParamsSettings;
%             t.epsilonPerFinal   = obj.epsilons(iter);
%             t.epsilonPerInitial = obj.epsilons(iter);
         %   obj.topSet.incrementalSchemeSettings.targetParamsSettings = t;
        end
        
        function setLevelSetValue(obj)
            t = obj.topSet;
            t.designVarSettings.initialCase = 'given';
            t.designVarSettings.levelSetCreatorSettings.value = obj.levelSet;   
            obj.topSet = t;            
        end
        
        function setInitialVolume(obj)
           t = obj.topSet;            
           vI = obj.initialVolume;
           t.incrementalSchemeSettings.targetParamsSettings.VfracInitial = vI;
           obj.topSet = t;                       
        end
        
        function setNsteps(obj,iter)
           t = obj.topSet;            
           nSt = obj.nSteps(iter);
           t.incrementalSchemeSettings.nSteps = nSt;
           t.incrementalSchemeSettings.targetParamsSettings.nSteps = nSt;           
           obj.topSet = t;              
        end
        
        function setMaxIter(obj,iter)
           t = obj.topSet;            
           mI = obj.maxIters(iter);
           t.optimizerSettings.maxIter = mI;           
           obj.topSet = t;              
        end        
        
        function setOptimizer(obj,iter)
            if iter ~= 1
                obj.topSet.optimizerSettings.type = 'AlternatingPrimalDual';            
            end                
        end
        
        function solve(obj)
            optProblem = TopOpt_Problem(obj.topSet);
            optProblem.dualVariable.value = obj.lambda;          
            optProblem.computeVariables(); 
            obj.lambda = optProblem.dualVariable.value;
            obj.levelSet = optProblem.designVariable.value;
            optProblem.constraint.shapeFunctions{1}.computeCost();
            v = optProblem.constraint.shapeFunctions{1}.value;            
            obj.initialVolume = v ;
        end
    end
    
end