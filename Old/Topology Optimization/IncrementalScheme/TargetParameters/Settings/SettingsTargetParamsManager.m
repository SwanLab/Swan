classdef SettingsTargetParamsManager < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsTargetParamsManager.json'
    end
    
    properties (Access = public)
        nSteps
        constrInitial
        constrFinal
        optimalityInitial
        optimalityFinal
        
        VfracInitial
        VfracFinal
        
        epsilonInitial
        epsilonFinal
        
        epsilonPerInitial
        epsilonPerFinal
        
        epsilonIsoInitial
        epsilonIsoFinal
        
        stressNormExponentInitial
        stressNormExponentFinal
    end
    
    properties (GetAccess = private, SetAccess = public)
        mesh
    end
    
    methods (Access = public)
        
        function obj = SettingsTargetParamsManager(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
            obj.setIfEmpty('stressNormExponentInitial',2);
            obj.setIfEmpty('stressNormExponentFinal',2);                        
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            if ~isempty(obj.mesh)
                obj.initEpsilons();
            end
        end
        
        function initEpsilons(obj)
            L = obj.mesh.computeCharacteristicLength();
            D = obj.mesh.computeMeanCellSize();
            obj.setIfEmpty('epsilonInitial',D);
            obj.setIfEmpty('epsilonFinal',obj.epsilonInitial);
            obj.setIfEmpty('epsilonPerInitial',L);
            obj.setIfEmpty('epsilonPerFinal',obj.epsilonInitial);
            obj.setIfEmpty('epsilonIsoInitial',D);
            obj.setIfEmpty('epsilonIsoFinal',obj.epsilonInitial);
        end
        
        function setIfEmpty(obj,prop,b)
            if isempty(obj.(prop))
                obj.(prop) = b;
            end
        end
        
    end
    
    methods
        
        function set.mesh(obj,m)
            obj.mesh = m;
            obj.init();
        end
        
    end
    
end