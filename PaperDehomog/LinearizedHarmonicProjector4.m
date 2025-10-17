classdef LinearizedHarmonicProjector4 < handle

    properties (Access = private)
        eta
        internalDOFs
        massMatrixSS
        massMatrixGG
        massMatrixBB
        stiffnessMatrixBB
        fB
        fBV
        fS
        perimeter
    end

    properties (Access = private)
        mesh
        boundaryNodes
        density
    end

    methods (Access = public)

        function obj = LinearizedHarmonicProjector4(cParams)
            obj.init(cParams);
            obj.fB = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.fBV = LagrangianFunction.create(obj.mesh, 2, 'P1');
            obj.fS = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.createInternalDOFs();                        
            obj.eta = (10*obj.mesh.computeMeanCellSize)^2;  
            obj.perimeter = obj.density.*(1-obj.density);%ConstantFunction.create(1,obj.mesh);%
            obj.massMatrixBB      = IntegrateLHS(@(u,v) DP(v,obj.perimeter.*u),obj.fB,obj.fB,obj.mesh,'Domain',4);
            obj.massMatrixSS      = IntegrateLHS(@(u,v) DP(v,u),obj.fS,obj.fS,obj.mesh,'Domain',4);
            obj.stiffnessMatrixBB = IntegrateLHS(@(u,v) DP(Grad(v),(1-obj.perimeter).*Grad(u)),obj.fB,obj.fB,obj.mesh,'Domain');
        end

        function b = solveProblem(obj,bBar,b)
            b    = project(b,obj.fB.order);
            %bBar = project(bBar,obj.fB.order);
            RHS = obj.computeRHS(bBar);
            LHS = obj.computeLHS(b);
            nInt = size(obj.internalDOFs,2);
            x = [b.fValues(:);zeros(nInt,1)];
            res = norm(LHS*x - RHS)/norm(x);
            [resL,resH,resB,resG] = obj.evaluateResidualNorms(bBar,b);
            i = 1;
            theta = 0.5;
            while res(i) > 1e-12
                xNew   = LHS\RHS;
                x = theta*xNew + (1-theta)*x;
                b   = obj.createVectorFromSolution(x);
                b   = obj.projectInUnitBall(b);
                LHS = obj.computeLHS(b);
                i   = i+1;
                res(i) = norm(LHS*x - RHS)/norm(x);
                [resL(i),resH(i),resB(i),resG(i)] = obj.evaluateResidualNorms(bBar,b);
                disp(['iter ',num2str(i),' residual ',num2str(res(i))])
                close all
                plotVector(b)
                %fig = figure; set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
                drawnow                
            end
            figure()
            plot(1:i,log10([res; resL; resH; resB; resG]))
            legend('LHS*x-RHS','resDistance','resHarmonic','resUnitBall','resGradient')

            % Mbb  = obj.massMatrixBB;
            % Kbb  = obj.stiffnessMatrixBB;            
            % A    = Mbb + obj.eta*Kbb;
            % rhsB = IntegrateRHS(@(v) DP(v,obj.perimeter.*bBar),bBar,obj.mesh,'Domain',4);   
            % Z  = sparse(obj.fB.nDofs,obj.fB.nDofs);            
            % xV    = [A Z; Z' A]\rhsB;
            % s.fValues = reshape(full(xV),[],2);
            % s.mesh    = obj.mesh;
            % s.order   = obj.fB.order;
            % b = LagrangianFunction(s);
            b = project(b,'P1');
        end

        function [resLnorm,resHnorm,resBnorm,resGnorm] = evaluateResidualNorms(obj,bBar,b)
            [resL,resH,resB,resG] = obj.evaluateAllResiduals(bBar,b);
            resLnorm = Norm(resL,'L2');
            resHnorm = Norm(resH,'L2');
            resBnorm = Norm(resB,'L2');
            resGnorm = Norm(resG,'L2');
        end

        function b = createVectorFromSolution(obj,x)
            nB = obj.fB.nDofs;
            bV = x(1:2*nB);
            s.fValues = reshape(bV,[],2);
            s.mesh    = obj.mesh;
            s.order   = obj.fB.order;
            b = LagrangianFunction(s);
        end

        function [resL,resH,resB,resG] = evaluateAllResiduals(obj,bBar,b)
            resL = project(DP(b-bBar,obj.perimeter.*(b-bBar)),'P1');
            resH = obj.evaluateHarmonicResidual(b);                
            resB = project(norm(b)- 1,'P1');
            resG = project(DDP(Grad(b),(1-obj.perimeter).*Grad(b)),'P1');
        end

        function resH = evaluateHarmonicResidual(obj,b)
            bs = b.getVectorFields();            
            f = (-Grad(bs{1}).*bs{2}+Grad(bs{2}).*bs{1});
            rhsV = IntegrateRHS(@(v) DP(Grad(v),f),obj.fS,obj.mesh,'Domain',4);
            rhsV(obj.boundaryNodes) = 0;
            Mss = obj.massMatrixSS;
            Kss = IntegrateLHS(@(u,v) DP(Grad(v),Grad(u)),obj.fS,obj.fS,obj.mesh,'Domain');
            eps = (1*obj.mesh.computeMeanCellSize)^2;
            LHS = Mss+eps*Kss;
            hf = LHS\rhsV;
            resH = obj.createFunction(full(hf),obj.fS.order);
        end


    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh             = cParams.mesh;
            obj.boundaryNodes    = cParams.boundaryMesh;
            obj.density          = cParams.density;            
        end
  
        function createInternalDOFs(obj)
            bNodes = obj.boundaryNodes;
            iDOFs  = setdiff(1:obj.fS.nDofs,bNodes);
            obj.internalDOFs = iDOFs;
        end

        function f = createFunction(obj,fV,order)
            s.fValues = fV;
            s.mesh    = obj.mesh;
            s.order   = order;
            f = LagrangianFunction(s);
        end

        function LHS = computeLHS(obj,b)
            Mbb  = obj.massMatrixBB;
            Kbb  = obj.stiffnessMatrixBB;
            nInt = size(obj.internalDOFs,2);
            bs = b.getVectorFields();
            normB  = (norm(b) - 1);
            MnormB = IntegrateLHS(@(u,v) DP(v,2*normB.*u),obj.fB,obj.fB,obj.mesh,'Domain',4);
            Kb1 = IntegrateLHS(@(u,v) DP(Grad(v),bs{1}.*Grad(u)),obj.fB,obj.fS,obj.mesh,'Domain',4);
            Kb2 = IntegrateLHS(@(u,v) DP(Grad(v),bs{2}.*Grad(u)),obj.fB,obj.fS,obj.mesh,'Domain',4);
            
            %Nb1 = computeAdvection(bs{1},obj.fB,obj.fB,obj.mesh,4);
            %Nb2 = computeAdvection(bs{2},obj.fB,obj.fB,obj.mesh,4);

            %v = LagrangianFunction.create(obj.mesh, 2, 'P1');

            %Nb1t = IntegrateLHS(@(u,v) u.*DP(Grad(v),Grad(bs{1})),obj.fB,obj.fS,obj.mesh,'Domain',4)';            
            %Nb2t = IntegrateLHS(@(u,v) u.*DP(Grad(v),Grad(bs{2})),obj.fB,obj.fS,obj.mesh,'Domain',4)'; 

            Nb1 = IntegrateLHS(@(u,v) DP(Grad(u),Grad(bs{1})).*v,obj.fB,obj.fS,obj.mesh,'Domain',4);
            Nb2 = IntegrateLHS(@(u,v) DP(Grad(u),Grad(bs{2})).*v,obj.fB,obj.fS,obj.mesh,'Domain',4);

            iDOFs = obj.internalDOFs;
            Kb1 = Kb1(:,iDOFs);
            Kb2 = Kb2(:,iDOFs);
            Nb1 = Nb1(:,iDOFs);
            Nb2 = Nb2(:,iDOFs);
            Z  = sparse(obj.fB.nDofs,obj.fB.nDofs);
            Zh = sparse(nInt,nInt);
            eps = 0;
            A  = Mbb + obj.eta*Kbb + eps*MnormB;
            LHS = [A          ,          Z, (-Kb2+Nb2);...
                   Z'         ,          A,  (Kb1-Nb1);...
                   (-Kb2+Nb2)', (Kb1-Nb1)',         Zh];
        end

        function RHS = computeRHS(obj,bBar)
            rhsB = IntegrateRHS(@(v) DP(v,obj.perimeter.*bBar),obj.fBV,obj.mesh,'Domain',3);
            rhsB = reshape(rhsB,2,[])';
            rhsB = rhsB(:);
            rhsH = zeros(size(obj.internalDOFs,2),1);
            RHS  = [rhsB;rhsH];
        end

    end

    methods (Access = private, Static)

        function vNF = projectInUnitBall(vF)
            v    = vF.fValues;
            norm = vecnorm(v,2,2);
            v    = v./norm;
            vNF = LagrangianFunction.create(vF.mesh,vF.ndimf,vF.order);
            vNF.setFValues(v);
        end

    end    

end