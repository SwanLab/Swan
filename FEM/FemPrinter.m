classdef FemPrinter < handle
    
    properties (Access = private)
        mesh
        quad
        iter
        fields
        ndim
        pdim
        ptype
        type
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = FemPrinter(cParams)
            obj.mesh   = cParams.mesh;
            obj.quad   = cParams.quad;
            obj.iter   = cParams.iter;
            obj.fields = cParams.fields;
            obj.ndim   = cParams.ndim;
            obj.ptype  = cParams.ptype;
            obj.pdim   = cParams.pdim;
            obj.type   = cParams.type;
        end
        
        function print(obj,fileName)
            dI = obj.createPostProcessDataBase(fileName);
            dI.ndim   = obj.ndim;
            dI.pdim   = obj.pdim;
            dI.ptype  = obj.ptype;
            dI.name = 'field';
            p = Postprocess(obj.type,dI);
            q        = obj.quad;
            d.fields = obj.fields;
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