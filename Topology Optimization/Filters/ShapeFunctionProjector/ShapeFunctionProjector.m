classdef ShapeFunctionProjector < handle
    
    properties (Access = protected)
        mesh
     end

    methods (Access = public, Static)
       
        function obj = create(cParams)
            f = ShapeFunctionProjectorFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end
        
    end
    
  methods (Access = public, Abstract)
      project(obj)
  end
    
end