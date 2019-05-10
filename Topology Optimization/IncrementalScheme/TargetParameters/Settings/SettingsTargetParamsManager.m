classdef SettingsTargetParamsManager < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsTargetParamsManager.json'
    end
    
    properties (Access = public)
        nSteps
        VfracInitial
        VfracFinal
        constrInitial
        constrFinal
        optimalityInitial
        optimalityFinal
        
        epsilonInitial
        epsilonFinal
        epsilonPerInitial
        epsilonPerFinal
        epsilonIsotropyInitial
        epsilonIsotropyFinal
    end
    
    methods (Access = public)
        
        function obj = SettingsTargetParamsManager(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
end