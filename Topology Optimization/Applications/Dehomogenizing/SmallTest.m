classdef SmallTest < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        mesh
        alpha
        orientation
    end
    
    methods (Access = public)
        
        function obj = SmallTest()
            obj.init();
            obj.createMesh();
            obj.createAngle();
            obj.createOrientation();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end

        function createMesh(obj)
            s.coord = [0 0; 0.5 0.5; 0 1; -0.5 0.5];
            s.connec = [1 2 3; 1 3 4];
            obj.mesh = Mesh.create(s); 
        end

        function createAngle(obj)
            s.fHandle = @(x) obj.createAlphaValues(x);
            s.mesh    = obj.mesh;
            s.ndimf   = 1;
            a         = AnalyticalFunction(s);
            obj.alpha = a;%project(a,'P1');%a             
        end

        function f = createAlphaValues(obj,x)
            x1 = x(1,:,:);
            x2 = x(2,:,:);
            x10 = (max(x1(:))+min(x1(:)))/2;
            x20 = 0;            
            f = atan2(x2-x20 +0.0*(max(x2(:))),x1-x10);
          %  isLeft = x1 < (min(x1(:))+ max(x1(:)))/2;
          %  f(isLeft) = f(isLeft) + pi;
        end

        function createOrientation(obj)
            nDim = obj.mesh.ndim;
            for iDim = 1:nDim

                a1 = project(cos(project(obj.alpha,'P1D')),'P1D'); 
                a2 = project(sin(project(obj.alpha,'P1D')),'P1D'); 
               % aV = reshape(',1,[])';
                s.fValues = [a1.fValues, a2.fValues];
                s.mesh    = obj.mesh;
                s.order   = 'P1D';
                s.ndimf   = 2;
                aF = LagrangianFunction(s); 

                coord= obj.mesh.computeBaricenter();
                x1 = coord(:,1);
                x2 = coord(:,2);                
                isLeft = x1 < (min(x1(:))+ max(x1(:)))/2;                

                aFelem = aF.getFvaluesByElem;
                tElem = aFelem;
                tElem(2,:,isLeft) = aFelem(1,:,isLeft);
                tElem(1,:,isLeft) = -aFelem(2,:,isLeft);

                
                s.operation = @(xV) obj.createOrientationFunction(iDim,xV);
                s.ndimf     = 2;
                s.mesh      = obj.mesh;
                aF = DomainFunction(s);     

                %better change the sign than adding 180degrees
                obj.orientation{iDim} = aF;
            end
        end        

        function or = createOrientationFunction(obj,iDim,xV)
            alphaV = obj.alpha.evaluate(xV);
            if iDim == 1                
                or(1,:,:) = cos(alphaV);
                or(2,:,:) = sin(alphaV);

                x = obj.mesh.computeXgauss(xV);
                x1 = x(1,:,:);
                x2 = x(2,:,:);                
                isLeft = x1 < (min(x1(:))+ max(x1(:)))/2;
            %  f(isLeft) = f(isLeft) + pi;               
              %  or(1,isLeft)  = -sin(alphaV(isLeft));
              %  or(2,isLeft) = cos(alphaV(isLeft));
            else
                alphaV = obj.alpha.evaluate(xV);
                or(1,:,:) = -sin(alphaV);
                or(2,:,:) = cos(alphaV);
            end
        end        
        
    end
    
end