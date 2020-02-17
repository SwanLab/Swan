classdef Mesh_Unfitted_Factory < handle
        
    methods (Access = public)
        
        function mesh = create(obj,cParams)            
            if obj.shallBeComposite(cParams)
                mesh = Mesh_Unfitted_Composite(cParams);
            else
                mesh = Mesh_Unfitted_Single(cParams);
            end
        end
        
    end
    
    methods (Access = private)
                
        function shall = shallBeComposite(obj,cParams)
            type = cParams.unfittedType;
            if cParams.includeBoxContour && strcmp(type,"INTERIOR")
               % warning('Contours are always included for INTERIOR mesh type.')
                shall = false;
            else
                shall = cParams.includeBoxContour;
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