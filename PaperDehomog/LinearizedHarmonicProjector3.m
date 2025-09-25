classdef LinearizedHarmonicProjector3 < handle

    properties (Access = private)
        eta
        internalDOFs
        massMatrixSS
        massMatrixGG
        massMatrixBB
        stiffnessMatrixBB
        fB
        fS
        fG
        perimeter
    end

    properties (Access = private)
        mesh
        boundaryNodes
        density
    end

    methods (Access = public)

        function obj = LinearizedHarmonicProjector3(cParams)
            obj.init(cParams);
            obj.fB = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.fS = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.fG = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.createInternalDOFs();                        
            obj.eta = (20*obj.mesh.computeMeanCellSize)^2;  
            obj.perimeter = ConstantFunction.create(1,obj.mesh);%obj.density.*(1-obj.density);
            obj.massMatrixBB      = IntegrateLHS(@(u,v) DP(v,obj.perimeter.*u),obj.fB,obj.fB,obj.mesh,'Domain');
            obj.massMatrixGG      = IntegrateLHS(@(u,v) DP(v,u),obj.fG,obj.fG,obj.mesh,'Domain');
            obj.massMatrixSS      = IntegrateLHS(@(u,v) DP(v,u),obj.fS,obj.fS,obj.mesh,'Domain');
            obj.stiffnessMatrixBB = IntegrateLHS(@(u,v) DP(Grad(v),(obj.perimeter).*Grad(u)),obj.fB,obj.fB,obj.mesh,'Domain');
        end

        function b = solveProblem(obj,bBar,b)
            RHS = obj.computeRHS(bBar);
            LHS = obj.computeLHS(b);
            nInt = size(obj.internalDOFs,2);
            x = [b.fValues(:);zeros(nInt,1);zeros(obj.fG.nDofs,1)];
            res = norm(LHS*x - RHS)/norm(x);
            [resL,resH,resB,resG] = obj.evaluateResidualNorms(bBar,b);
            i = 1;
            theta = 0.5;
            while res(i) > 1e-6
                xNew   = LHS\RHS;
                x = theta*xNew + (1-theta)*x;
                b   = obj.createVectorFromSolution(x);
                LHS = obj.computeLHS(b);
                i   = i+1;
                res(i) = norm(LHS*x - RHS)/norm(x);
                [resL(i),resH(i),resB(i),resG(i)] = obj.evaluateResidualNorms(bBar,b);
                disp(['iter ',num2str(i),' residual ',num2str(res(i))])
            end
            figure()
            plot(1:i,log([res; resL; resH; resB; resG]))
            legend('LHS*x-RHS','resDistance','resHarmonic','resUnitBall','resGradient')
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
            resG = project(DDP(Grad(b),Grad(b)),'P1');
        end

        function resH = evaluateHarmonicResidual(obj,b)
            bs = b.getVectorFields();            
            f = (-Grad(bs{1}).*bs{2}+Grad(bs{2}).*bs{1});
            rhsV = IntegrateRHS(@(v) DP(Grad(v),f),obj.fS,obj.mesh,4);
            rhsV(obj.boundaryNodes) = 0;
            Mss = obj.massMatrixSS;
            hf = Mss\rhsV;
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
            Mb1 = IntegrateLHS(@(u,v) DP(v,bs{1}.*u),obj.fB,obj.fG,obj.mesh,'Domain');
            Mb2 = IntegrateLHS(@(u,v) DP(v,bs{2}.*u),obj.fB,obj.fG,obj.mesh,'Domain');
            Kb1 = IntegrateLHS(@(u,v) DP(Grad(v),bs{1}.*Grad(u)),obj.fB,obj.fS,obj.mesh,'Domain');
            Kb2 = IntegrateLHS(@(u,v) DP(Grad(v),bs{2}.*Grad(u)),obj.fB,obj.fS,obj.mesh,'Domain');
            Nb1 = IntegrateLHS(@(u,v) DP(Grad(bs{1}),v.*Grad(u)),obj.fB,obj.fS,obj.mesh,'Domain');
            Nb2 = IntegrateLHS(@(u,v) DP(Grad(bs{2}),v.*Grad(u)),obj.fB,obj.fS,obj.mesh,'Domain'); 
            iDOFs = obj.internalDOFs;
            Kb1 = Kb1(:,iDOFs);
            Kb2 = Kb2(:,iDOFs);
            Nb1 = Nb1(:,iDOFs);
            Nb2 = Nb2(:,iDOFs);
            Z  = sparse(obj.fB.nDofs,obj.fB.nDofs);
            Zh = sparse(nInt,nInt);
            Zsg = sparse(nInt,obj.fG.nDofs);
            Zgg = sparse(obj.fG.nDofs,obj.fG.nDofs);
            A  = Mbb + obj.eta*Kbb;
            LHS = [A          ,          Z, (-Kb2+Nb2),Mb1;...
                   Z          ,          A,  (Kb1-Nb1),Mb2;...
                   (-Kb2+Nb2)', (Kb1-Nb1)',         Zh,Zsg;...
                   Mb1'       ,       Mb2',       Zsg',Zgg];
        end

        function RHS = computeRHS(obj,bBar)
            rhsB = IntegrateRHS(@(v) DP(v,obj.perimeter.*bBar),bBar,obj.mesh,3);
            rhsH = zeros(size(obj.internalDOFs,2),1);
            rhsU = IntegrateRHS(@(v) v,obj.fG,obj.mesh,2);
            RHS  = [rhsB;rhsH;rhsU];
        end

    end

end