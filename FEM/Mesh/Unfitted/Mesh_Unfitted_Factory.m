classdef Mesh_Unfitted_Factory < handle
    
    properties (Access = private)
        includeBoxContour
        nargin
    end
    
    methods (Access = public, Static)
        
        function meshUnfitted = create(cParams,PropertyName,PropertyValue)
            obj = Mesh_Unfitted_Factory();
            obj.nargin = nargin;
            
            obj.determineFlagState(PropertyName,PropertyValue);
            
            if obj.shallBeComposite(cParams)
                meshUnfitted = Mesh_Unfitted_Composite(cParams);
            else
                meshUnfitted = Mesh_Unfitted(cParams);
            end
        end
        
    end
    
    methods (Access = private)
        
        function determineFlagState(obj,PropertyName,PropertyValue)
            switch obj.nargin
                case 1
                    obj.includeBoxContour = false;
                case 3
                    obj.assignPropertyValue(PropertyName,PropertyValue);
            end
        end
        
        function assignPropertyValue(obj,PropertyName,PropertyValue)
            obj.checkProperty(PropertyName,PropertyValue);
            obj.includeBoxContour = PropertyValue;
        end
        
        function shall = shallBeComposite(obj,cParams)
            type = cParams.unfittedType;
            if obj.includeBoxContour && strcmp(type,"INTERIOR")
                warning('Contours are always included for INTERIOR mesh type.')
                shall = false;
            else
                shall = obj.includeBoxContour;
            end
        end
        
    end
    
    methods (Access = private, Static)
        
        function checkProperty(PropertyName,PropertyValue)
            if ~strcmp(PropertyName,'includeBoxContour')
                error('Invalid property. Valid property names: includeBoxContour')
            else
                if ~islogical(PropertyValue)
                    error('includeBoxContour must be a logical value');
                end
            end
        end
        
    end
end