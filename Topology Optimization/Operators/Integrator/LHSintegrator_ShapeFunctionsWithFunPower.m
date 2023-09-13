classdef LHSintegrator_ShapeFunctionsWithFunPower < handle

    properties (Access = private)
        mesh
        fun
        exponent
        designVarType
    end

    properties (Access = private)
        interpolationLHS
        quadrature
    end

    methods (Access = public)

        function obj = LHSintegrator_ShapeFunctionsWithFunPower(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            switch obj.designVarType
                case 'Density'
                    lhs = obj.computeElementalLHS();
                    LHS = obj.assembleIntegrand(lhs);
                case 'LevelSet'

            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh             = cParams.mesh;
            obj.fun              = cParams.fun;
            obj.exponent         = cParams.exponent;
            obj.designVarType    = cParams.designVarType;
            obj.interpolationLHS = Interpolation.create(obj.mesh,'LINEAR');
        end
        
        function createQuadrature(obj)
            quadratureOrder = 'CUBIC';
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadratureOrder);
            obj.quadrature = quad;
        end

        function lhs = computeElementalLHS(obj)
            p      = obj.exponent;
            quad   = obj.quadrature;
            dVolu  = obj.mesh.computeDvolume(quad);
            nGaus  = obj.quadrature.ngaus;
            nElem  = size(dVolu,2);
            nNodes = obj.mesh.nnodeElem;
            nDim   = nNodes*obj.fun.ndimf;
            xg     = quad.posgp;
            obj.interpolationLHS.computeShapeDeriv(xg);
            shapes = obj.interpolationLHS.shape;
            funVals = obj.fun.evaluate(xg);
            lhs = zeros(nDim, nDim, nElem);

            for igauss = 1 :nGaus
                for inode= 1:nNodes
                    Ni = shapes(inode,igauss,:);
                    for jnode= 1:nNodes
                        Nj = shapes(jnode,igauss,:);
                        for iunkn= 1:obj.fun.ndimf
                            idof = obj.fun.ndimf*(inode-1)+iunkn;
                            jdof = obj.fun.ndimf*(jnode-1)+iunkn;
                            dvol = dVolu(igauss,:);
                            f  = squeeze(funVals(iunkn,igauss,:));
                            lhs(idof, jdof, :)= squeeze(lhs(idof,jdof,:)) ...
                                + (Ni^p*Nj*f.^(p-1)).*dvol';
                        end
                    end
                end
            end

        end

        function LHS = assembleIntegrand(obj,lhsElem)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assembleFunctions(lhsElem, obj.fun, obj.fun);
        end

    end

end