classdef SettingsFilter < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsFilter.json'
    end
    
    properties (Access = public)
        filterType
        domainType
        designVariable
        quadratureOrder
        femSettings
        designVarType
        mesh
    end
    
    methods (Access = public)
        
        function obj = SettingsFilter(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
%     methods (Access = private)
%         
%         function setFemSettingsMesh(obj)
%             if ~isempty(obj.designVar)
%                 obj.femSettings.mesh = obj.designVar.mesh;
%             end
%         end
%         
%     end
%     
%     methods
%         
% %         function set.designVar(obj,dV)
% %             obj.designVar = dV;
% %             obj.setFemSettingsMesh();
% %         end
%         
%     end
    
end