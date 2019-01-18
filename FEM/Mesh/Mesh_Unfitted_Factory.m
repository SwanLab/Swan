classdef Mesh_Unfitted_Factory < handle
    properties (Access = private)
        includeBoxContour
        shallBeComposite
        nargin
    end
    
    methods (Access = public, Static)
        function mesh_unfitted = create(meshType,meshBackground,interpolation_background,PropertyName,PropertyValue)
            obj = Mesh_Unfitted_Factory;
            obj.nargin = nargin;
            
            obj.determineFlagState(PropertyName,PropertyValue);
            obj.checkFlagStateConsistency(meshType);
            
            if obj.shallBeComposite
                mesh_unfitted = Mesh_Unfitted_Composite(meshType,meshBackground,interpolation_background);
            else
                mesh_unfitted = Mesh_Unfitted(meshType,meshBackground,interpolation_background);
            end
        end
    end
    
    methods (Access = private)
        function determineFlagState(obj,PropertyName,PropertyValue)
            if obj.nargin == 3
                obj.includeBoxContour = false;
            elseif obj.nargin == 5
                obj.assignPropertyValue(PropertyName,PropertyValue);
            end
        end
        
        function assignPropertyValue(obj,PropertyName,PropertyValue)
            obj.checkProperty(PropertyName,PropertyValue);
            obj.includeBoxContour = PropertyValue;
        end
        
        function checkFlagStateConsistency(obj,meshType)
           if obj.includeBoxContour && strcmp(meshType,"INTERIOR")
               warning('Contours are always included for INTERIOR mesh type.')
               obj.shallBeComposite = false;
           else
               obj.shallBeComposite = obj.includeBoxContour;
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