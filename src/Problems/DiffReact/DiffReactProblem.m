classdef DiffReactProblem < handle
    
    properties (GetAccess = public, SetAccess = protected)
        variables
        x
    end
    
    properties (Access = private)
        mesh
        solver
        epsilon
        LHStype
        problemData
    end

    methods (Access = public)
        
        function obj = DiffReactProblem(cParams)
            obj.init(cParams);
            obj.createSolver();
        end

        function computeVariables(obj,rhs)
            RHS = rhs;
            LHS = obj.computeLHS(obj.epsilon);
            obj.variables.x = obj.solver.solve(LHS,RHS);

            a.mesh = obj.mesh;
            a.fValues = obj.variables.x;
            a.order = 'P1';
            obj.x = LagrangianFunction(a);
        end
        
        function LHS = computeLHS(obj, epsilon)
            obj.epsilon = epsilon;
            vF = LagrangianFunction.create(obj.mesh,1,'P1');
            uF = LagrangianFunction.create(obj.mesh,1,'P1');
            ndof  = uF.nDofs;
            Mr    = sparse(ndof,ndof);
            if strcmp(obj.LHStype, "StiffnessMassBoundaryMass")
                Mr = IntegrateLHS(@(u,v) DP(v,u),vF,uF,obj.mesh,'Boundary',3);
            end
            K = IntegrateLHS(@(u,v) DP(Grad(v),Grad(u)),vF,uF,obj.mesh,'Domain');            
            M = IntegrateLHS(@(u,v) DP(v,u),vF,uF,obj.mesh,'Domain',3);
            LHS = (obj.epsilon^2).*K + M + obj.epsilon*Mr;
        end
       
        function print(obj,filename)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            s.quad = quad;
            s.mesh = obj.mesh;
            s.iter = 0;
            s.fields    = obj.variables.x;
            s.ptype     = 'DIFF-REACT';
            s.ndim      = 3;
            s.pdim      = obj.problemData.pdim;
            s.type      = 'ScalarNodal';
            fPrinter = FemPrinter(s);
            fPrinter.print(filename);
        end

    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh    = cParams.mesh;
            obj.LHStype = cParams.LHStype;
        end
        
        function createSolver(obj)
            s.type = 'DIRECT';
            obj.solver = Solver.create(s);
        end
        
    end

end