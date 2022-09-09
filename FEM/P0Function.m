classdef P0Function < FeFunction
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        quadrature
        dim
    end
    
    properties (Access = private)
       mesh 
       fValues
       quadOrder
    end
    
    methods (Access = public)
        
        function obj = P0Function(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createFvaluesByElem();
        end
        
        function fNodal = projectToLinearNodalFunction(obj)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS();
            fNodal = (LHS\RHS);
        end

        function plot(obj, m, f)
%             disMesh = m.createDiscontinousMesh();
%             p0funct = obj.preparePlottingFunction(disMesh, f);
%             coor = disMesh.coord;
%             conn = disMesh.connecM;
%             trisurf(conn, coor(:,1), coor(:,2), p0funct)
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.fValues   = cParams.fValues;
            obj.quadOrder = 'LINEAR';
%             obj.connec    = cParams.connec;
        end
        
        function q = createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadOrder);
            obj.quadrature = q;
        end
        
        function LHS = computeLHS(obj)
            sf.mesh  = obj.mesh;
            sf.ndimf = 1; 
            sf.interpolationOrder = 'LINEAR';
            sf.quadratureOrder = 'QUADRATIC';
            field = Field(sf);

            s.mesh         = obj.mesh;
            s.type         = 'MassMatrix';
            s.field        = field;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end
        
        function RHS = computeRHS(obj)
            s.type      = 'ShapeFunction';
            s.mesh      = obj.mesh;
            s.meshType  = obj.mesh.type;
            s.fType     = 'Gauss';
            s.fGauss    = obj.computeFgauss();
            s.xGauss    = obj.computeXgauss();
            s.quadOrder = obj.quadOrder;
            s.npnod     = obj.mesh.nnodes;
            s.globalConnec = obj.mesh.connec;
            rhs = RHSintegrator.create(s);

            fG = obj.computeFgauss();
            xG = obj.computeXgauss();
            RHS = rhs.compute(fG,xG);
        end
        
        function x = computeXgauss(obj)
            xG = obj.quadrature.posgp;
            x = repmat(xG,[1,1,obj.mesh.nelem]);
        end
        
        function f = computeFgauss(obj)
            ngaus = obj.quadrature.ngaus;
            fV(1,:) = obj.fValues;
            f = repmat(fV,[ngaus,1]);
        end

        function p0F = preparePlottingFunction(obj, dM, f)
            nnodeElem = dM.nnodeElem;
            fRepeated = zeros(size(f,1), nnodeElem);
            for iNode = 1:nnodeElem
                fRepeated(:,iNode) = f;
            end
            fRepeated = transpose(fRepeated);
            p0F = fRepeated(:);
        end

        function createFvaluesByElem(obj)
        end

    end
    
end