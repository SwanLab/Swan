classdef MappingComputer < handle


    properties (Access = private)
        mesh
        interpolator
        dilatedOrientation
        testFunction
    end

    methods (Access = public)

        function obj = MappingComputer(cParams)
            obj.init(cParams);
        end

        function uF = compute(obj)
            LHS = obj.computeStiffnessMatrix();
            In  = obj.interpolator;
            for iDim = 1:obj.mesh.ndim
                RHS = obj.computeRHS(iDim);
                uC  = obj.solveSystem(LHS,RHS);
                uC  = In*uC;                          
                uV(iDim,:,:) = reshape(uC,obj.mesh.nnodeElem,[]);
            end
            s.mesh    = obj.mesh;
            s.fValues = reshape(uV,obj.mesh.ndim,[])';
            s.order   ='P1D';                                  
        
            uF = LagrangianFunction(s);     

            uF = project(abs(uF),'P1D');            

            uFf = uF.getVectorFields();
            uFf{1} = uFf{1} - Mean(uFf{1},2);
            uFf{2} = uFf{2} - Mean(uFf{2},2);
            s.mesh    = obj.mesh;
            s.fValues(:,1) = uFf{1}.fValues;
            s.fValues(:,2) = uFf{2}.fValues;
            s.order   ='P1D';                                  
            uF = LagrangianFunction(s);             

   
        end
        %% 

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.dilatedOrientation = cParams.dilatedOrientation;
            obj.interpolator       = cParams.interpolator;
            obj.testFunction       = LagrangianFunction.create(obj.mesh,1,'P1D');
        end

        function K = computeStiffnessMatrix(obj)
            s.mesh  = obj.mesh;
            s.type  = 'StiffnessMatrix';
            s.test  = obj.testFunction;
            s.trial = obj.testFunction;
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
        end

        function RHS = computeRHS(obj,iDim)
            aI = obj.dilatedOrientation{iDim};
            s.mesh            = obj.mesh;
            s.quadratureOrder = 3;
            s.type            = 'ShapeDerivative';
            test = obj.testFunction;
            rhs  = RHSintegrator.create(s);
            rhsV = rhs.compute(aI,test);
            In   = obj.interpolator;
            RHS = In'*rhsV;          
        end

        function u = solveSystem(obj,LHS,RHS)
            In = obj.interpolator;            
            LHS = In'*LHS*In;

            I = ones(size(LHS,1),1);
            LHS = [LHS,I;I',0];  
            RHS = [RHS;0];


            a.type = 'DIRECT';
            s = Solver.create(a);
            u = s.solve(LHS,RHS);  
           
            u = u(1:end-1);
        end

    end

end
