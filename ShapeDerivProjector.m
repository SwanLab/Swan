classdef ShapeDerivProjector < handle

    properties (Access = private)
        mesh
    end
    
    methods (Access = public)

        function obj = ShapeDerivProjector()
            obj.init();
        end

        function xFun = project(obj,m,f,df)
            obj.mesh = m;
            LHS     = obj.computeLHS();
            for i = 1:length(df)
                RHSdfun(:,i) = obj.computeRHSdFun(df{i});
            end
            % 
            RHSfun = obj.computeRHSFun(f);
            RHS     = RHSfun + RHSdfun;
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = xProj;
            xFun = P1Function(s);
        end

    end

    methods (Access = private)

        function init(obj)
            
        end
        
        function LHS = computeLHS(obj)
            LHSM    = obj.computeMassMatrix();
            LHSK    = obj.computeStiffnessMatrix();
            epsilon = 1.5*obj.mesh.computeMeanCellSize();
            LHS     = LHSM + epsilon^2*LHSK;
        end

        function LHSM = computeMassMatrix(obj)
            s.test  = P1Function.create(obj.mesh,1);
            s.trial = P1Function.create(obj.mesh,1);
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrix';
            s.quadratureOrder = 'QUADRATIC';
            lhs = LHSintegrator.create(s);
            LHSM = lhs.compute();
        end

        function LHSK = computeStiffnessMatrix(obj)
            s.type = 'StiffnessMatrix';
            s.mesh = obj.mesh;
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.quadratureOrder = 'CONSTANT';
            lhs = LHSintegrator.create(s);
            LHSK = lhs.compute();
        end

        function RHS = computeRHSFun(obj,fun)
            s.type     = 'ShapeDivergence';
            s.mesh     = obj.mesh;
            s.quadratureOrder = 'QUADRATIC';            
            rhs = RHSintegrator.create(s);
            rhs = rhs.compute(fun);
            RHS = rhs.fValues;
        end

        function RHS = computeRHSdFun(obj,dfun)
            s.type     = 'ShapeFunction';
            s.mesh     = obj.mesh;
            s.test    = P1Function.create(obj.mesh,1); 
            s.quadratureOrder = 'QUADRATIC';            
            rhs = RHSintegrator.create(s);
            RHS = rhs.compute(dfun);
        end  

    end

end