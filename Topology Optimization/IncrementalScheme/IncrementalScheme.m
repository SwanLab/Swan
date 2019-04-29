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
            obj.setWhetherShallDisplayStep(cParams.shallPrintIncremental);
        end
        
        function createTargetParams(obj,cParams)
            settingsTargetParams = obj.editTargetParamsSettings(cParams);
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
        
        function setupEpsilons(obj,mesh,settingsTargetParams)
            L = mesh.computeCharacteristicLength();
            D = mesh.computeMeanCellSize();
            obj.assignWithBackup('epsilonInitial',settingsTargetParams.epsilonInitial,D);
            obj.assignWithBackup('epsilonFinal',settingsTargetParams.epsilonFinal,obj.epsilonInitial);
            obj.epsilonPerInitial = L;
            obj.epsilonPerFinal = obj.epsilonInitial;
            obj.assignWithBackup('epsilonIsoInitial',settingsTargetParams.epsilonIsotropyInitial,nan);
            obj.assignWithBackup('epsilonIsoFinal',settingsTargetParams.epsilonIsotropyFinal,nan);
            
        end
        
        function settingsTargetParams = editTargetParamsSettings(obj,cParams)
            obj.setupEpsilons(cParams.mesh,cParams.settingsTargetParams);
            
            settingsTargetParams = cParams.settingsTargetParams;
            settingsTargetParams.nSteps = obj.nSteps;
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