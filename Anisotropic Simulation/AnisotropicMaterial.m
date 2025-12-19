classdef AnisotropicMaterial < handle

    properties (Access = private)
        mesh
    end

    methods (Access = public)

        function obj = AnisotropicMaterial(cParams)
            obj.mesh = cParams.mesh;
        end

        
        function C = evaluate(obj,xV,C_voigt)
            nGauss = size(xV,2);
            nElem  = obj.mesh.nelem;
            C_tensor = convert2Tensor(C_voigt,'Constitutive');
            C = repmat(C_tensor,[1 1 1 1 nGauss nElem]);
        end

        function plot(obj,mesh)
            s.mesh = mesh;
            s.projectorType = 'P1D';
            proj = Projector.create(s);
            p1fun = proj.project(obj);
            p1fun.plot();
        end

    end


end
