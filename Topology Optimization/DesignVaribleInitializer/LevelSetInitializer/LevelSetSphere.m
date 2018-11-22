classdef LevelSetSphere < LevelSetSphereNdim
    

    
    methods (Access = public)
        
        function obj = LevelSetSphere(input)
            obj.fracRadius = 1-1e-6;
            obj.compute(input);
        end
    end    
    
    methods (Access = protected)

        
%         function computeDesignVariable(obj)
%             phi = obj.levelSet;
%              switch obj.optimizerName
%                 case {'SLERP','HAMILTON-JACOBI'}
%                     obj.x = phi;
%                 otherwise
%                     initial_holes = ceil(max(phi,0))>0;
%                     obj.x = obj.ini_design_value*ones(obj.lsSize);
%                     obj.x(initial_holes) = obj.hole_value;
%              end                        
%         end
        
    end
    
end


