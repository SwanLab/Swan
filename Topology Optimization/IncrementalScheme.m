classdef IncrementalScheme < handle
    
    properties (GetAccess = public, SetAccess = private)
        iStep
        nSteps
        targetParams
    end
    
    properties (Access = private)
        targetParamsManager
        
        epsilonInitial
        epsilonFinal
        epsilonPerInitial
        epsilonPerFinal
        epsilonIsoInitial
        epsilonIsoFinal
        
        shallDisplayStep
    end
    
    methods (Access = public)
        
        function obj = IncrementalScheme(cParams)
            obj.init(cParams);
            obj.createTargetParams(cParams);
        end
        
        function next(obj)
            obj.incrementStep();
            obj.updateTargetParams();
        end
        
        function display(obj)
            disp(['Incremental Scheme - Step: ',int2str(obj.iStep),' of ',int2str(obj.nSteps)]);
        end
        
        function itDoes = hasNext(obj)
            if obj.iStep < obj.nSteps
                itDoes = true;
            else
                itDoes = false;
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.iStep = 0;
            obj.nSteps = cParams.nSteps;
            obj.setupEpsilons(cParams);
            obj.setWhetherShallDisplayStep(cParams.shallPrintIncremental);
        end
        
        function createTargetParams(obj,cParams)
            settingsTargetParams = obj.createTargetSettings(cParams);            
            obj.targetParamsManager = TargetParamsManager(settingsTargetParams);
            obj.targetParams = obj.targetParamsManager.targetParams;
        end
        
        function updateTargetParams(obj)
            obj.targetParamsManager.update(obj.iStep);
        end
        
        function incrementStep(obj)
            obj.iStep = obj.iStep + 1;
            if obj.shallDisplayStep
                obj.display();
            end
        end
        
        function setupEpsilons(obj,cParams)
            L = cParams.mesh.computeCharacteristicLength();
            D = cParams.mesh.computeMeanCellSize();
            obj.assignWithBackup('epsilonInitial',cParams.epsilonInitial,D);
            obj.assignWithBackup('epsilonFinal',cParams.epsilonFinal,obj.epsilonInitial);
            obj.epsilonPerInitial = L;
            obj.epsilonPerFinal = obj.epsilonInitial;
            obj.assignWithBackup('epsilonIsoInitial',cParams.epsilonIsoInitial,nan);
            obj.assignWithBackup('epsilonIsoFinal',cParams.epsilonIsoFinal,nan);
            
        end
        
        function settingsTargetParams = createTargetSettings(obj,cParams)
            settingsTargetParams = SettingsTargetParamsManager();
            
            settingsTargetParams.nSteps = obj.nSteps;
            settingsTargetParams.VfracInitial = cParams.VfracInitial;
            settingsTargetParams.VfracFinal = cParams.VfracFinal;
            settingsTargetParams.constrInitial = cParams.constrInitial;
            settingsTargetParams.constrFinal = cParams.constrFinal;
            settingsTargetParams.optimalityInitial = cParams.optimalityInitial;
            settingsTargetParams.optimalityFinal = cParams.optimalityFinal;
            
            settingsTargetParams.epsilonInitial = obj.epsilonInitial;
            settingsTargetParams.epsilonFinal = obj.epsilonFinal;
            settingsTargetParams.epsilonPerInitial = obj.epsilonPerInitial;
            settingsTargetParams.epsilonPerFinal = obj.epsilonPerFinal;
            settingsTargetParams.epsilonIsotropyInitial = obj.epsilonIsoInitial;
            settingsTargetParams.epsilonIsotropyFinal = obj.epsilonIsoFinal;
        end
        
        function setWhetherShallDisplayStep(obj,flag)
            obj.shallDisplayStep = flag;
            if isempty(flag)
                obj.shallDisplayStep = true;
            end
        end
        
        function assignWithBackup(obj,prop,a,b)
            if ~isempty(a)
                obj.(prop) = a;
            else
                obj.(prop) = b;
            end
        end
        
    end
    
end