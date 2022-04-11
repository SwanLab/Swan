classdef FemPrinter < handle
    
    properties (Access = public)
        mesh
        quad
        iter
        variables
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = FemPrinter(cParams)
            obj.mesh = cParams.mesh;
            obj.quad = cParams.quad;
            obj.iter = cParams.iter;
            obj.variables = cParams.variables;
        end
        
        function print(obj,fileName)
            dI = obj.createPostProcessDataBase(fileName);
            postprocess = Postprocess('Elasticity',dI);
            q = obj.getQuadrature();
            d.fields = obj.variables;
            d.quad = q;
            postprocess.print(obj.iter,d);
        end
        
    end
    
    methods (Access = private)
        
        function d = createPostProcessDataBase(obj,fileName)
            dI.mesh    = obj.mesh;
            dI.outName = fileName;
            ps = PostProcessDataBaseCreator(dI);
            d = ps.getValue();
        end
        
    end
    
end