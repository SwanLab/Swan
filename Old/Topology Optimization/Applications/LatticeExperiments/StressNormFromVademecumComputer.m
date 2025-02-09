classdef StressNormFromVademecumComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        designVariable
        stressNormFunctional
    end
    
    properties (Access = private)
        m1
        m2
        alpha
        mesh
        vademecumName
        fileName
        pNorm
    end
    
    methods (Access = public)
        
        function obj = StressNormFromVademecumComputer(cParams)
            obj.init(cParams)
        end
        
        function v = compute(obj)
           obj.createDesignVariable();
           obj.createStressNormFunctional();
           v = obj.computeStressNorm();
       end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.m1            = cParams.m1;
            obj.m2            = cParams.m2;
            obj.alpha         = cParams.alpha;
            obj.vademecumName = cParams.vademecumName;
            obj.mesh          = cParams.mesh;
            obj.fileName      = cParams.fileName;
            obj.pNorm         = cParams.pNorm;
        end
        
        function valueR = computeStressNorm(obj)
            obj.stressNormFunctional.computeFunction();
            value0 = obj.stressNormFunctional.value0;
            nP = length(obj.pNorm);
            valueR = zeros(nP,1);
            for iStep = 1:nP
                obj.stressNormFunctional.setPnorm(obj.pNorm(iStep))
                obj.stressNormFunctional.computeFunction();
                v = obj.stressNormFunctional.value;
                valueR(iStep) = v*(value0);
            end
        end
        
        function createDesignVariable(obj)
            s = SettingsDesignVariable();
            s.mesh                  = Mesh_Total(obj.mesh);
            s.type                  = 'MicroParams';
            s.initialCase           = 'full';
            s.scalarProductSettings = obj.createScalarProductParams();
            s.creatorSettings       = obj.createCreatorSettings();
            desVar = DesignVariable.create(s);
            desVar.alpha = obj.alpha;
            obj.designVariable = desVar;
        end
        
        function s = createScalarProductParams(obj)
            s.epsilon = [];
            s.mesh = obj.mesh;
        end
        
        function s = createCreatorSettings(obj)
            s.m1 = obj.m1;
            s.m2 = obj.m2;
            s.homogSettings.fileName = obj.vademecumName;
            s.homogSettings.type     = 'ByVademecum';
        end
        
        function createStressNormFunctional(obj)
            s.designVariable = obj.designVariable;
            s.type = 'stressNorm';
            s.targetParameters = [];
            s.filterParams     = obj.createFilterSettings();
            s.femSettings      = obj.createFemSettings();
            s.homogVarComputer = obj.createHomogVarComputer();
            s.targetParameters = obj.createTargetParameters();
            obj.stressNormFunctional = ShapeFunctional.create(s);
        end
        
        function femS = createFemSettings(obj)
            femS.fileName = obj.fileName;
            femS.scale = 'MACRO';
            femS.mesh = obj.mesh;
        end
        
        function s = createFilterSettings(obj)
            s.filterType = 'PDE';
            s.femSettings = obj.createFemSettings();
        end
        
        function targetParams = createTargetParameters(obj)
            s = SettingsTargetParamsManager;
            targetParamsManager = TargetParamsManager(s);
            targetParams = targetParamsManager.targetParams;
        end
        
        function homogVar = createHomogVarComputer(obj)
            s.type = 'ByVademecum';
            s = SettingsHomogenizedVarComputer.create(s);
            s.fileName = obj.vademecumName;
            s.nelem = obj.mesh.nelem;
            homogVar = HomogenizedVarComputer.create(s);
        end
        
    end
    
end