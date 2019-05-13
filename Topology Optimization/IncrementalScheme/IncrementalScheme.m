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
            targetParamsSettings = obj.editTargetParamsSettings(cParams);
            obj.targetParamsManager = TargetParamsManager(targetParamsSettings);
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
        
        function setupEpsilons(obj,mesh,targetParamsSettings)
            L = mesh.computeCharacteristicLength();
            D = mesh.computeMeanCellSize();
            obj.assignWithBackup('epsilonInitial',targetParamsSettings.epsilonInitial,D);
            obj.assignWithBackup('epsilonFinal',targetParamsSettings.epsilonFinal,obj.epsilonInitial);
            obj.epsilonPerInitial = L;
            obj.epsilonPerFinal = obj.epsilonInitial;
            obj.assignWithBackup('epsilonIsoInitial',targetParamsSettings.epsilonIsotropyInitial,nan);
            obj.assignWithBackup('epsilonIsoFinal',targetParamsSettings.epsilonIsotropyFinal,nan);
            
        end
        
        function targetParamsSettings = editTargetParamsSettings(obj,cParams)
            obj.setupEpsilons(cParams.mesh,cParams.targetParamsSettings);
            
            targetParamsSettings = cParams.targetParamsSettings;
            targetParamsSettings.nSteps = obj.nSteps;
            targetParamsSettings.epsilonInitial = obj.epsilonInitial;
            targetParamsSettings.epsilonFinal = obj.epsilonFinal;
            targetParamsSettings.epsilonPerInitial = obj.epsilonPerInitial;
            targetParamsSettings.epsilonPerFinal = obj.epsilonPerFinal;
            targetParamsSettings.epsilonIsotropyInitial = obj.epsilonIsoInitial;
            targetParamsSettings.epsilonIsotropyFinal = obj.epsilonIsoFinal;
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