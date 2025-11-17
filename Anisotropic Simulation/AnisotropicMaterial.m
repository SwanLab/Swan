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

            % Matriz de rigidez en notación Voigt (2D, 3x3)
            % (ejemplo: valores arbitrarios, cambia por tu matriz homogeneizada)

            % Mapeo Voigt 2D
            % 1->11, 2->22, 3->12
            map = [1 1; 2 2; 1 2];

            % Inicializar tensor de 4º orden (2D)
            C_tensor = zeros(2,2,2,2);

            for I = 1:3
                for J = 1:3
                    i = map(I,1); j = map(I,2);
                    k = map(J,1); l = map(J,2);

                    % Factores de corte (Voigt usa ingenieril γ12, con factor 2)
                    facI = 1; if I==3, facI = sqrt(2); end
                    facJ = 1; if J==3, facJ = sqrt(2); end

                    val = C_voigt(I,J) / (facI*facJ);

                    % Asignar respetando simetrías
                    C_tensor(i,j,k,l) = val;
                    C_tensor(j,i,k,l) = val;
                    C_tensor(i,j,l,k) = val;
                    C_tensor(j,i,l,k) = val;

                    C_tensor(k,l,i,j) = val;
                    C_tensor(l,k,i,j) = val;
                    C_tensor(k,l,j,i) = val;
                    C_tensor(l,k,j,i) = val;
                end
            end
            C_tensor(1,2,2,1) = 0;
            C_tensor(2,1,1,2) = 0;
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
