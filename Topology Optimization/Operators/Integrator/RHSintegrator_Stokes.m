classdef RHSintegrator_Stokes < handle

    properties (Access = private)
        mesh
        velocityFun
        pressureFun
        forcesFormula
        quadrature
    end

    methods (Access = public)

        function obj = RHSintegrator_Stokes(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function rhs = integrate(obj)
            rhs = obj.computeRHS();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.velocityFun   = cParams.velocityFun;
            obj.pressureFun   = cParams.pressureFun;
            obj.forcesFormula = cParams.forcesFormula;
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC'); % ehhh
            obj.quadrature = q;
        end

        function RHS = computeRHS(obj)
            Fext = obj.computeVolumetricFext();
            g = obj.computeVelocityDivergence();
            rhs{1,1} = Fext;
            rhs{2,1} = g;
            RHS = obj.assemble(rhs);
        end

        function Fext = computeVolumetricFext(obj)
            shapesV = obj.velocityFun.computeShapeFunctions(obj.quadrature);
            dvol = obj.mesh.computeDvolume(obj.quadrature)';
            ngaus = size(dvol,2);
            nNode = size(shapesV, 1);
            nDimV = obj.velocityFun.ndimf;
            Fext = zeros(nNode*nDimV,1,obj.mesh.nelem);

            f = obj.calculateForcesFromExpression();

            for igaus=1:ngaus
                for inode=1:nNode
                    for iunkn=1:nDimV
                        elemental_dof = inode*nDimV-nDimV+iunkn; %% dof per guardar el valor de la integral
                        shape = shapesV(inode,igaus);
                        fvalue = f(iunkn,igaus,:);
                        v= squeeze(shape.*fvalue);
                        Fext(elemental_dof,1,:)= squeeze(Fext(elemental_dof,1,:)) + v(:).*dvol(:,igaus);

                    end
                end
            end
        end

        function f = calculateForcesFromExpression(obj)
            nGaus  = obj.quadrature.ngaus;
            xV = obj.quadrature.posgp;
            xGauss = obj.mesh.computeXgauss(xV);
            nElem = obj.mesh.nelem;
            nDimf = obj.velocityFun.ndimf;
            F = zeros(nDimf,nGaus,nElem);
            for iElem = 1:nElem
                ind=1;
                for iGaus = 1:nGaus
                    xG = xGauss(:,iGaus,iElem);
                    pos_node = num2cell(xG);
                    fCell = obj.forcesFormula(pos_node{:});
                    fMat = cell2mat(fCell);
                    F(:,iGaus,iElem) = fMat;
                    ind=ind+length(fMat);
                end
            end
            f = F;
        end

        function g = computeVelocityDivergence(obj)
            shp = obj.pressureFun.computeShapeFunctions(obj.quadrature);
            nDofE = size(shp,1);
            g = zeros(nDofE,1,obj.mesh.nelem);
        end

        function RHS = assemble(obj, rhs)
            s.fun = [];
            assembler = AssemblerFun(s);
            RHS = assembler.assembleVectorStokes(rhs, obj.velocityFun, obj.pressureFun);
        end

    end

end