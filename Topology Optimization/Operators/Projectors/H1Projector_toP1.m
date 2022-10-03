classdef H1Projector_toP1 < Projector

    properties (Access = private)
        fieldMass
        fieldStiffness
    end
    
    methods (Access = public)

        function obj = H1Projector_toP1(cParams)
            obj.init(cParams);
            obj.createFieldMass();
            obj.createFieldStiffness();
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.type    = obj.mesh.type;
            s.connec  = obj.mesh.connec;
            s.fValues = xProj;
            xFun = P1Function(s);
        end

    end

    methods (Access = private)

        function createFieldMass(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1; 
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            obj.fieldMass = Field(s);
        end

        function createFieldStiffness(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1; 
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'CONSTANT';
            obj.fieldStiffness = Field(s);
        end
        
        function LHS = computeLHS(obj)
            LHSM    = obj.computeMassMatrix();
            LHSK    = obj.computeStiffnessMatrix();
            epsilon = obj.mesh.computeMeanCellSize();
            LHS     = LHSM + epsilon^2*LHSK;
        end

        function LHSM = computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.fieldMass;
            lhs = LHSintegrator.create(s);
            LHSM = lhs.compute();
        end

        function LHSK = computeStiffnessMatrix(obj)
            s.type  = 'StiffnessMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.fieldStiffness;
            lhs = LHSintegrator.create(s);
            LHSK = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);
            obj.mesh.interpolation.computeShapeDeriv(xV);
            shapes = permute(obj.mesh.interpolation.shape,[1 3 2]);
            conne = obj.mesh.connec;

            nGaus = quad.ngaus;
            nFlds = fun.ndimf;
            nElem = obj.mesh.nelem;
            nNode = size(conne,2);
            nDofs = obj.mesh.nnodes;

            fLoc  = zeros(nNode,nElem,nFlds);
            fGaus = fun.evaluate(xV);
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(iField,igaus,:));
                    for inode = 1:nNode
                        Ni = shapes(inode,igaus);
                        fL(:,1) = squeeze(fLoc(inode,:,iField));
                        int = Ni*fG.*dVg;
                        fLoc(inode,:,iField) = fL + int;
                    end
                end
                for iElem = 1:nElem
                    dofs = conne(iElem,:);
                    f(dofs,iField) = f(dofs,iField) + fLoc(:,iElem,iField);
                end
            end
            RHS = f;
        end

        function q = createRHSQuadrature(obj, fun)
            ord = obj.determineQuadratureOrder(fun);
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end
        
    end

end