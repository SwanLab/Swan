classdef LinearizedHarmonicProjector < handle
    
    properties (Access = private)
        dim
        massMatrix
        advectionMatrixX
        advectionMatrixY
        reducedAdvectionMatrixX
        reducedAdvectionMatrixY
        LHS
        solver
        dofsInElem
    end    

    properties (Access = private)
       mesh   
       boundaryMesh
    end
    
    methods (Access = public)
        
        function obj = LinearizedHarmonicProjector(cParams)
            obj.init(cParams);            
            obj.createDimension();
            obj.computeDofConnectivity();            
            obj.computeMassMatrix();
            obj.createSolver()
        end
        
        function lambda0 = computeInitalLambda(obj)
            b    = obj.boundaryMesh;
            nInt = setdiff(1:obj.dim.npnod,b);
            lambda0 = zeros(length(nInt),1);
        end

        function [v,lambda] = solveProblem(obj,v0,vH)
            
            vH = obj.projectUnitBall(vH);
            obj.plotOrientation(vH,1);
            i = 1;
            isErrorLarge = true;
            while isErrorLarge
                obj.plotOrientation(vH,2);
                [Ax,Ay] = obj.computeAdvectionMatrix(vH);
                for idim = 1:2
                    switch idim
                        case 1
                            A = Ax;
                        case 2
                            A = Ay;
                    end
                    Ared = obj.computeReducedAdvectionMatrix(A);
                    obj.computeLHS(Ared);
                    lhs    = obj.LHS;
                    rhs    = obj.computeRHS(v0(:,idim));
                    sol    = obj.solver.solve(lhs,rhs);
                    v(:,idim) = sol(1:obj.dim.npnod,1);
                    lambda(:,idim) = sol(obj.dim.npnod+1:end,1);
                end  
                err(i) = norm(vH(:)-v(:))/norm(v(:));
                isErrorLarge = err(i) > 1e-2;
                i = i + 1;                
                figure(99)
                plot(log10(err))                               
                vH = v;                
            end



        end

        function tp = projectUnitBall(obj,t)
            uP = UnitBallProjector([]);
            tp = uP.project(t);
        end

        function plotOrientation(obj,t,iFigure)
            figure(100)
            subplot(2,2,iFigure)            
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = t(:,1);
            ty = t(:,2);
            quiver(x,y,tx,ty);
            tp = obj.projectUnitBall(t);
            tx = tp(:,1);
            ty = tp(:,2);
            subplot(2,2,2*(2-1)+iFigure)            
            quiver(x,y,tx,ty);            
        end

        function optPrim = computePrimalOptimaility(obj,lambda,v,vH)
            M        = obj.massMatrix;  
            Ared     = obj.reducedAdvectionMatrix;  
            gradPrim = M*(v-vH) + (Ared'*lambda);
            optPrim  = norm(gradPrim);
        end

        function vH = projectByDual(obj,v)
            Kred = obj.reducedAdvectionMatrix;            
            M    = obj.massMatrix;            
            lhs = Kred*(M\Kred');
            rhs = Kred*v;            
            lambda  = obj.solver.solve(lhs,rhs);
            vH = v - M\(Kred'*lambda);
        end    

        function cost = computeCost(obj,v,vH)
            M    = obj.massMatrix;     
            cost = 0;
            for i = 1:size(v,2)
                vI   = v(:,i);
                incV = vI- vH(:,i);
                cost = cost + incV'*M*incV/(vI'*M*vI);            
            end
        end

        function optDual = computeDualOptimality(obj,vH)
            Ared = obj.reducedAdvectionMatrix;            
            optDual = norm(Ared*vH);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.mesh             = cParams.mesh;
           obj.boundaryMesh     = cParams.boundaryMesh;
        end
        
        function createDimension(obj)
            q = Quadrature();
            q = q.set(obj.mesh.type);
            s.mesh = obj.mesh;
            s.pdim = '1D';
            s.ngaus = q.ngaus;
            d = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end        
        
        function computeMassMatrix(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod        = obj.mesh.npnod;
            s.type         = 'MassMatrix';
            s.dim          = obj.dim;
            s.quadType     = 'QUADRATIC';
            s.dofsInElem   = obj.dofsInElem;    
            lhs = LHSintegrator.create(s);
            M = lhs.compute();
            %M = diag(sum(M));
            %M = eye(size(M,1));
            obj.massMatrix = M;
        end

        function [LHSX,LHSY] = computeAdvectionMatrix(obj,vH)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod        = obj.mesh.npnod;
            s.type         = 'AdvectionMatrix';
            s.dim          = obj.dim;
            s.dofsInElem   = obj.dofsInElem; 
            s.b            = vH;
            lhs = LHSintegrator.create(s);
            [LHSX,LHSY] =  lhs.compute();
       end        

       function Ared = computeReducedAdvectionMatrix(obj,A)
           b    = obj.boundaryMesh;
           nInt = setdiff(1:obj.dim.npnod,b);
           Ared = A(nInt,:);
       end

       function createSolver(obj)
            s = Solver.create();
            obj.solver = s;
       end

       function computeDofConnectivity(obj)
           connec = obj.mesh.connec;
           ndimf  = obj.dim.ndimField;
           nnode  = obj.dim.nnode;
           dofsElem  = zeros(nnode*ndimf,size(connec,1));
           for inode = 1:nnode
               for iunkn = 1:ndimf
                   idofElem   = obj.nod2dof(inode,iunkn);
                   globalNode = connec(:,inode);
                   idofGlobal = obj.nod2dof(globalNode,iunkn);
                   dofsElem(idofElem,:) = idofGlobal;
               end
           end
           obj.dofsInElem = dofsElem;
       end

       function  computeLHS(obj,A)
           M    = obj.massMatrix;
           Z    = obj.computeZeroFunction();
           lhs  = [M,A';A,Z];
           obj.LHS = lhs;
       end

       function Z = computeZeroFunction(obj)
           b    = obj.boundaryMesh;
           nInt = setdiff(1:obj.dim.npnod,b);           
           Z    = zeros(length(nInt),length(nInt));           
       end

        function rhs = computeRHS(obj,v)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            s.mesh  = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod = obj.dim.npnod;
            s.dim   = obj.dim;
            s.type = 'SIMPLE';
            int = Integrator.create(s);
            rhs = int.integrateFnodal(v,q.order);
            b = obj.boundaryMesh;
            nInt = setdiff(1:obj.dim.npnod,b);
            Z   = zeros(length(nInt),1);
            rhs = [rhs;Z];
        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = obj.dim.ndimField;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end        

    end
    
end