classdef FemPrinter < handle
    
    properties (Access = private)
        mesh
        quad
        iter
        variables
        ndim
        pdim
        ptype
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = FemPrinter(cParams)
            obj.mesh = cParams.mesh;
            obj.quad = cParams.quad;
            obj.iter = cParams.iter;
            obj.variables = cParams.variables;
            obj.ndim  = cParams.ndim;
            obj.ptype = cParams.ptype;
            obj.pdim  = cParams.pdim;
        end
        
        function print(obj,fileName)
            dI = obj.createPostProcessDataBase(fileName);
            dI.ndim   = obj.ndim;
            dI.pdim   = obj.pdim;
            dI.ptype  = obj.ptype;            
            p = Postprocess('ElasticityMicro',dI);
            q        = obj.quad;
            d.fields = obj.variables;
            d.quad   = q;
            p.print(obj.iter,d);
        end
        
    end
    
    methods (Access = private)
        
        function d = createPostProcessDataBase(obj,fileName)
            dI.mesh        = obj.mesh;
            dI.outFileName = fileName;
            ps = PostProcessDataBaseCreator(dI);
            d = ps.create();
        end
        
    end
    
end