classdef Mesh_Unfitted_Factory < handle
    
    properties (Access = private)
        includeBoxContour
    end
    
    methods (Access = public, Static)
        
        function meshUnfitted = create(cParams)
            obj = Mesh_Unfitted_Factory();
            obj.init(cParams);
            
            if obj.shallBeComposite(cParams)
                meshUnfitted = Mesh_Unfitted_Composite(cParams);
            else
                meshUnfitted = Mesh_Unfitted(cParams);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.includeBoxContour = cParams.includeBoxContour;
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