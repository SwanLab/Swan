classdef PDEShapeDerivProjector < handle

    properties (Access = private)
        mesh
    end
    
    methods (Access = public)

        function obj = PDEShapeDerivProjector()
            obj.init();
        end

        function xFun = project(obj,m,f,df,u,p)
            obj.mesh = m;
            LHS     = obj.computeLHS();
            RHS     = obj.computeRHSterms(f,df,u,p);
            % 
            RHS     = RHSfun + RHSdfun;
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = xProj;
            xFun = P1Function(s);
        end

        function RHS = computeRHSterms(obj,f,df,u,p)
           RHS1 = obj.computeFirstRHS(f);
           RHS2 = obj.computeSecondRHS(u,p);
           RHS3 = obj.computeThirdRHS();
           RHS4 = obj.computeFourthRHS(p);
           RHS = RHS1 + RHS2 + RHS3 + RHS4;
        end

        function rhs = computeFirstRHS(obj,f)
            rhs = obj.computeRHSFunTimesDivergence(f);
        end        

        function rhs = computeSecondRHS(obj,u,p)
            s.f1 = u;
            s.f2 = p;
            s.operation = @(f1,f2) f1.*f2;
            s.ndimf = 1 ;
            s.quad  = obj.createQuadrature();
            f = ComposedGradFunction(s);            
            rhs = obj.computeRHSFunTimesDivergence(f);    
        end

        function rhs = computeThirdRHS(obj,u,p)
             
        end

        function rhs = computeFourthRHS(obj,p)
          
        end        


        function q = createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('CUBIC');
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

        function RHS = computeRHSFunTimesDivergence(obj,fun)
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