classdef SettingsTargetParamsManager < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsTargetParamsManager'
    end
    
    properties (Access = public)
        nSteps
        Vfrac_initial
        Vfrac_final
        constr_initial
        constr_final
        optimality_initial
        optimality_final
        
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
            elseif nargin > 1
                obj.nSteps = varargin{1};
                obj.Vfrac_initial = varargin{2};
                obj.Vfrac_final = varargin{3};
                obj.constr_initial = varargin{4};
                obj.constr_final = varargin{5};
                obj.optimality_initial = varargin{6};
                obj.optimality_final = varargin{7};
                
                obj.epsilonInitial = varargin{8};
                obj.epsilonFinal = varargin{9};
                obj.epsilonPerInitial = varargin{10};
                obj.epsilonPerFinal = varargin{11};
                obj.epsilonIsotropyInitial = varargin{12};
                obj.epsilonIsotropyFinal = varargin{13};
            end
        end
        
    end
    
end