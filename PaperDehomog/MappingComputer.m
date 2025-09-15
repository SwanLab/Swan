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
                uC  = obj.solveSaddleSystem(LHS,RHS);
                uD  = In*uC;                          
                uV(iDim,:,:) = reshape(uD,obj.mesh.nnodeElem,[]);
                uCF(:,iDim) = uC;
            end

            % 
           %  s.mesh    = obj.mesh;
           %  s.fValues = uCF;
           %  s.order   ='P1';                                          
           %  uFC = LagrangianFunction(s);    
           % 
           %   s.mesh    = obj.mesh;
           %   s.fValues = (uFC.fValues);
           %   s.order   ='P1';                                  
           %   uFC = LagrangianFunction(s);   
           % 
           % %  uFCD = project(uFC,'P1D');

          
            s.mesh    = obj.mesh;
            s.fValues = reshape(uV,obj.mesh.ndim,[])';
            s.order   ='P1D';                                          
            uF = LagrangianFunction(s);                

             
            s.mesh    = obj.mesh;
            s.fValues = (uF.fValues);
            s.order   ='P1D';                                  
            uF = LagrangianFunction(s);  


            % uFf = uF.getVectorFields();
            % uFf{1} = uFf{1} - Mean(uFf{1},2);
            % uFf{2} = uFf{2} - Mean(uFf{2},2);
            %  s.mesh    = obj.mesh;
            %  s.fValues(:,1) = uFf{1}.fValues;
            %  s.fValues(:,2) = uFf{2}.fValues;
            %  s.order   ='P1D';
            %  uF = LagrangianFunction(s);
            % 
            % 
        %    dif = uFC-uF;
        %    I   = ConstantFunction.create(1,obj.mesh);
       %     a = L2norm(dif)/L2norm(I);
            
             % 
             % 
            uF = uF;
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
            s.quadratureOrder = 4;
            lhs = LHSIntegrator.create(s);
            K = lhs.compute();
        end

        function RHS = computeRHS(obj,iDim)
            aI = obj.dilatedOrientation{iDim};
            s.mesh            = obj.mesh;
            s.quadratureOrder = 4;
            s.type            = 'ShapeDerivative';
            s.test = obj.testFunction;
            rhs  = RHSIntegrator.create(s);
            rhsV = rhs.compute(aI);
            In   = obj.interpolator;
            RHS = In'*rhsV;          
        end

        function u = solveSaddleSystem(obj,LHS,RHS)
            M = obj.createMassMatrix();
            eta = 1e-2;%1e-2;%1e-15;%1e-2;%0.00000000000001;
            In = obj.interpolator;            
            LHS = In'*LHS*In+eta*In'*M*In;
            %LHS = In'*LHS*In+eta*M;
            %I = ones(size(LHS,1),1);
            %LHS = [LHS,I;I',0];  
           %RHS = [RHS;0];
            a.type = 'DIRECT';
            s = Solver.create(a);
            u = s.solve(LHS,RHS);             
            %u = u(1:end-1);
        end

        function M = createMassMatrix(obj)
            s.type            = 'MassMatrix';
            s.mesh            = obj.mesh;
            s.test            = obj.testFunction;            
            s.trial           = obj.testFunction;%LagrangianFunction.create(obj.mesh,1,'P1');
           % s.test            = LagrangianFunction.create(obj.mesh,1,'P1');
           % s.trial           = LagrangianFunction.create(obj.mesh,1,'P1');
            s.quadratureOrder = 4;
            LHS               = LHSIntegrator.create(s);
            M = LHS.compute();
        end

    end

end
