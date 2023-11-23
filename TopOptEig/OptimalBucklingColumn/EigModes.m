classdef EigModes < handle
       
    properties (Access = private)
        mesh
        V
        D
        v1t
        v2t
        eigModesPlotter
        lambda
        mode1Disp
        mode2Disp
        freeNodes        
    end

    properties (Access = private)
        dim
        sectionVariables
        stiffnessMatComputer
        bendingMatComputer
        inertiaMoment
        youngModulus
    end

    methods (Access = public)
        
        function obj = EigModes(cParams)
            obj.init(cParams)
            obj.createDimensions();
            obj.createBoundaryConditions();
            obj.createStiffnessMatrix();
            obj.createBendingMatrix();
        end             

        function l = provideEigenValue(obj)
            obj.computeEigenModesAndValues();            
            obj.computeLambda();   
            l = obj.lambda;
        end

       function grad = provideDerivative(obj,eigNum)
            obj.reorderModes(obj.lambda,obj.V,obj.D);
            Belem =  obj.bendingMatComputer.elementalBendingMatrix;
            eigV1 = obj.D(1,1)
            eigV2 = obj.D(2,2)
            difEigs = abs(eigV2-eigV1);
            if difEigs > 1 
                dfdx = obj.computeSimpleEig(Belem);
            else 
                dfdx = obj.computeDoubleEig(Belem);
            end
            grad = dfdx(eigNum,:);
        end

    end

    methods (Access = private)
        
        function dfdx = computeSimpleEig(obj,Belem)
            [v1,v2] = obj.computeEigenVectorByDof();         
            vBv(:,1) = obj.productMatrix(v1,v1,Belem);
            vBv(:,2) = obj.productMatrix(v2,v2,Belem);            
            nVar     = obj.sectionVariables.nDesVarElem;            
            vBv      = repmat(vBv,nVar,1);            
            dI       = obj.sectionVariables.computeInertiaDerivative();            

            dfdx = -dI.*vBv;
            dfdx(end+1,:) = 1;
            dfdx = dfdx';
        end

        function dfdx = computeDoubleEig(obj,Belem)
            [v1,v2] = obj.computeEigenVectorByDof();            
            vBv11 = obj.productMatrix(v1,v1,Belem);
            vBv12 = obj.productMatrix(v1,v2,Belem);
            vBv22 = obj.productMatrix(v2,v2,Belem);
            vBv   = obj.getEigenValues(vBv11,vBv22,vBv12,vBv12); 
            nVar  = obj.sectionVariables.nDesVarElem;            
            vBv    = repmat(vBv,nVar,1);

            dI   = obj.sectionVariables.computeInertiaDerivative();
            dfdx = -vBv.*dI;
            dfdx(end+1,:) = 1;
            dfdx = dfdx';
        end

        function [v1t,v2t] = computeEigenVectorByDof(obj)
            free  = obj.freeNodes;
            ndofe = obj.dim.ndofsElem;
            nnodeElem = obj.dim.nnodeElem;
            ndofn = obj.dim.ndofsElem/nnodeElem;
            nElem = obj.mesh.nelem;
            v1F    = zeros(obj.dim.ndofs,1);
            v2F    = zeros(obj.dim.ndofs,1);
            v1F(free,1) = obj.v1t;
            v2F(free,1) = obj.v2t;
            gElemt = ndofn*((1:nElem)-1); 
            v1t = zeros(nElem,ndofe);
            v2t = zeros(nElem,ndofe);
            for kDof = 1:ndofe
                indexK = gElemt(:)+kDof;
                v1t(:,kDof) = v1F(indexK);  
                v2t(:,kDof) = v2F(indexK);                          
            end
        end        

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh                 = cParams.mesh;
            %obj.designVariable       = cParams.designVariable;
            obj.sectionVariables    = cParams.sectionVariables;
            obj.inertiaMoment        = cParams.inertiaMoment;
            obj.youngModulus         = cParams.youngModulus;
        end
        
        function createDimensions(obj)
            s.mesh = obj.mesh;
            s.pdim = '2D';
            s.ngaus = 2;
            s.type = 'Vector';
            s.name = 'x';
            s.ndimf = 2;
            s.fieldName = 'u';
            d = DimensionVariables.create(s);
            d.compute();
            obj.dim = d;
        end
        
        function createBoundaryConditions(obj)
            d = obj.dim;
            fixnodes = union([1,2], [d.ndofs-1,d.ndofs]);
            %fixnodes = [1,d.ndofs-1];
            nodes = 1:d.ndofs;
            free  = setdiff(nodes,fixnodes);
            obj.freeNodes = free;
        end

        function createStiffnessMatrix(obj)
            s.type = 'StiffnessMatrixColumn';
            s.dim = obj.dim;
            s.mesh = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.freeNodes      = obj.freeNodes;
            K = LHSintegrator.create(s);
            obj.stiffnessMatComputer = K;
            obj.stiffnessMatComputer.compute();            
        end

        function createBendingMatrix(obj)
            s.type         = 'BendingMatrix';
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.inertiaMoment  = obj.inertiaMoment;
            s.youngModulus   = obj.youngModulus;
            %s.designVariable = obj.designVariable;
            s.sectionVariables = obj.sectionVariables;
            s.freeNodes      = obj.freeNodes;
            B = LHSintegrator.create(s);
            obj.bendingMatComputer = B;
        end

        function computeEigenModesAndValues(obj) 
            obj.bendingMatComputer.compute();            
            [Kfree,free]  = obj.stiffnessMatComputer.provideFreeStiffnessMatrix();
            obj.freeNodes = free;
            Bfree  = obj.bendingMatComputer.provideFreeBendingMatrix();
            obj.computeEigenFunctionAndValues(Bfree,Kfree);         
        end

        function computeLambda(obj)
            l = sort(diag(obj.D));
            obj.lambda = l;
        end

        function computeEigenFunctionAndValues(obj,B,K)
            [v,d] = eigs(B,K,2,'SM');
            obj.V  = v;
            obj.D  = d; 
        end
        
        function S = getEigenValues(obj,A1,A2,A12,A21)
            a = 1;
            b = -(A1 + A2);
            c = (A1.*A2)-(A12.*A21);
            lambd1 = (-b+sqrt(b.^2-4*a.*c))/(2*a);
            lambd2 = (-b-sqrt(b.^2-4*a.*c))/(2*a);
            S = sort([lambd1, lambd2],2);
        end

        function Wab = productMatrix(obj,Wa,Wb,Belem)
            d     = obj.dim;
            ndofe = d.ndofsElem;
            Wab   = zeros(obj.mesh.nelem,1);
            for iDof = 1:ndofe
                for jDof = 1:ndofe
                    Bij(:,1)   = squeeze(Belem(iDof,jDof,:));
                    w          = Wa(:,iDof).*Bij.*Wb(:,jDof);
                    Wab        = Wab + w;
                end
            end
        end

        function reorderModes(obj,lambda,V,D)
            if lambda(1)==D(1,1)
                V1=V(:,1);
                V2=V(:,2);
            else
                V1=V(:,2);
                V2=V(:,1);
            end
            obj.v1t = V1;
            obj.v2t = V2;             
        end

        function [m1disp,m2disp] = computeBucklingModes(obj,v1,v2)
            N = obj.mesh.nelem;
            Mode1=zeros(2*(N+1),1);
            Mode2=zeros(2*(N+1),1);
            for i=3:2*N
                Mode1(i)=v1(i-2);
                Mode2(i)=v2(i-2);
            end
            m1 = Mode1;
            m2 = Mode2;
            m1disp = m1(1:2:end);
            m2disp = m2(1:2:end);
        end             

    end
    
end