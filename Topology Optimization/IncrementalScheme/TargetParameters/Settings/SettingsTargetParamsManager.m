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
    end
    
    properties (GetAccess = private, SetAccess = public)
        mesh
    end
    
    methods (Access = public)
        
        function obj = SettingsTargetParamsManager(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.initEpsilons();
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
    
end