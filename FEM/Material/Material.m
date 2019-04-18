classdef Material < handle
    
    properties (Access = public)
        nelem
    end
    
    
    methods (Access = public, Static)
        function material = create(geometry,mesh)
            cParams.ptype = mesh.ptype;
            cParams.pdim  = mesh.pdim;
            cParams.nelem = geometry(1).interpolation.nelem;
            cParams.geometry = geometry;
            cParams.mesh  = mesh;
            f = MaterialFactory();
            material = f.create(cParams);
        end
    end
    
end