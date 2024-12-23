classdef SmallTest < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        meshSize
        xmin
        xmax
        ymin
        ymax
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
            obj.createOrientations();
            obj.createOrientedMapping();
        end
        
    end
    
    methods (Access = private)
        
         function init(obj)
            obj.meshSize = 0.032;
            obj.xmin = 0;
            obj.xmax = 2;
            obj.ymin = 0;
            obj.ymax = 1;
        end

        function createMesh(obj)
            h = obj.meshSize;
            xv = obj.xmin:h:obj.xmax;
            yv = obj.ymin:h:obj.ymax;
            [X,Y] = meshgrid(xv,yv);
            s.coord(:,1) = X(:);
            s.coord(:,2) = Y(:);
            [F,V] = mesh2tri(X,Y,zeros(size(X)),'f');
            s.coord  = V(:,1:2);
            s.connec = F;

            m = Mesh.create(s);
            obj.mesh = m;

            % s.coord = [0 0; 0.5 0.5; 0 1; -0.5 0.5];
            % s.connec = [1 2 3; 1 3 4];
            % obj.mesh = Mesh.create(s); 
        end
     
        function createOrientations(obj)
            coord = obj.mesh.coord;
            x1 = coord(:,1);
            x2 = coord(:,2);
            x10 = (max(x1(:))+min(x1(:)))/2;
            x20 = -0.5*max(x2(:));                                    
            r = sqrt((x1-x10).^2+(x2-x20).^2);
            fR = obj.normalize([(x1-x10)./r,(x2-x20)./r]);            
            fT = obj.normalize([-(x2-x20)./r,(x1-x10)./r]);
            obj.orientation{1} = obj.createOrientationField(fR);
            obj.orientation{2} = obj.createOrientationField(fT);
        end 

        function fN = normalize(obj,f)
            fN = f./vecnorm(f,2,2);
        end

        function f = createOrientationField(obj,fV)
            fD = obj.createP1DiscontinousOrientation(fV);
            f  = obj.computeAngleOrientation(fD);            
        end

        function fS = computeAngleOrientation(obj,f)
            aV = f.fValues;
            aV1 = aV(:,1);
            aV2 = aV(:,2);
            fN(:,1) = aV1.^2-aV2.^2;
            fN(:,2) = 2*aV1.*aV2;
            s.fValues = fN;
            s.mesh    = obj.mesh;
            s.order   = 'P1D';
            s.ndimf   = 2;
            fS = LagrangianFunction(s);            
        end

        function aF = createP1DiscontinousOrientation(obj,fV)
            s.fValues = fV;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            s.ndimf   = 2;
            aF = LagrangianFunction(s);
            aF = project(aF,'P1D');
        end

       
        function createOrientedMapping(obj)
            s.mesh = obj.mesh;
            s.orientation = obj.orientation;
            oM = OrientedMappingComputer(s);
            dCoord = oM.computeDeformedCoordinates()
        end
        
    end
    
end