classdef RHSintegrator_Stokes < handle

    properties (Access = private)
        mesh
        velocityField
        pressureField
        forcesFormula
    end

    methods (Access = public)

        function obj = RHSintegrator_Stokes(cParams)
            obj.init(cParams);
        end

        function rhs = integrate(obj)
            rhs = obj.computeRHS();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.velocityField = cParams.velocityField;
            obj.pressureField = cParams.pressureField;
            obj.forcesFormula = cParams.forcesFormula;
        end

        function RHS_elem = computeRHS(obj)
            Fext = obj.computeVolumetricFext();
            g = obj.computeVelocityDivergence();
            RHS_elem{1,1} = Fext;
            RHS_elem{2,1} = g;
%             RHS = AssembleVector(obj,RHS_elem);
        end

        function Fext = computeVolumetricFext(obj)
            geometry = obj.velocityField.geometry;
            shapesV  = obj.velocityField.interpolation.shape;
            dvol = geometry.dvolu;
            ngaus = size(dvol,2);
            nnode = obj.velocityField.dim.nnodeElem;
            nunkn = obj.velocityField.dim.ndimf;
            Fext = zeros(nnode*nunkn,1,obj.mesh.nelem);

            f = obj.calculateForcesFromExpression();

            for igaus=1:ngaus
                for inode=1:nnode
                    for iunkn=1:nunkn
                        elemental_dof = inode*nunkn-nunkn+iunkn; %% dof per guardar el valor de la integral
                        shape = shapesV(inode,igaus);
                        fvalue = f(iunkn,igaus,:);
                        v= squeeze(shape.*fvalue);
                        Fext(elemental_dof,1,:)= squeeze(Fext(elemental_dof,1,:)) + v(:).*dvol(:,igaus);

                    end
                end
            end
        end

        function f = calculateForcesFromExpression(obj)
            ngaus  = size(obj.velocityField.interpolation.shape,2);
            xGauss = obj.velocityField.xGauss;
            nelem = obj.mesh.nelem;
            for ielem = 1:nelem
                ind=1;
                for igaus = 1:ngaus
                    xG = xGauss(:,igaus,ielem);
                    pos_node = num2cell(xG);
                    fCell = obj.forcesFormula(pos_node{:});
                    fMat = cell2mat(fCell);
                    F(:,igaus,ielem) = fMat;
                    ind=ind+length(fMat);
                end
            end
            f = F;
        end

        function g = computeVelocityDivergence(obj)
            nunkn = obj.velocityField.dim.ndimf;
            nnode = obj.pressureField.dim.nnodeElem;
            g = zeros(nnode*nunkn,1,obj.mesh.nelem);
        end

    end

end