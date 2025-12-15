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

            dim=2;
            if dim == 2
                pairs = [1 1; 2 2; 1 2];
            elseif dim == 3
                pairs = [1 1; 2 2; 3 3; 2 3; 1 3; 1 2];
            else
                error('Tensor must be 2D or 3D.');
            end

            dimVoigt = size(pairs,1);
            C_tensor = zeros(dim, dim, dim, dim);
        
            for m = 1:dimVoigt
                for n = 1:dimVoigt
                    i = pairs(m,1); j = pairs(m,2);
                    k = pairs(n,1); l = pairs(n,2);
        
                    C_tensor(i,j,k,l) = C_voigt(m,n);
                    C_tensor(j,i,k,l) = C_voigt(m,n);
                    C_tensor(i,j,l,k) = C_voigt(m,n);
                    C_tensor(j,i,l,k) = C_voigt(m,n);
                end
            end
        
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
