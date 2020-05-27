classdef MeshPlotter < handle
    
    methods (Access = public, Static, Abstract)
       plot(mesh,ax)
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams) 
            f = MeshPlotterFactory();
            obj = f.create(cParams);
        end
    end
    
end

