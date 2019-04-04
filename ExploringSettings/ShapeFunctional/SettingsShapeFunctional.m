classdef SettingsShapeFunctional < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsShapeFunctional'
    end
    
    properties
        filterParams 
        filename 
        domainType 
        materialInteporlationParams
    end
    
     methods (Access = public)
        
        function obj = SettingsShapeFunctional(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
%                 case 2
%                     obj.unfittedType = varargin{1};
%                     obj.meshBackground = varargin{2};
%                     obj.interpolationBackground = Interpolation.create(obj.meshBackground,'LINEAR');
%                 case 3
%                     obj.unfittedType = varargin{1};
%                     obj.meshBackground = varargin{2};
%                     obj.interpolationBackground = varargin{3};
%                 case 4
%                     obj.unfittedType = varargin{1};
%                     obj.meshBackground = varargin{2};
%                     obj.interpolationBackground = varargin{3};
%                     obj.includeBoxContour = varargin{4};
            end
        end        
     end
end