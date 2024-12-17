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
            obj.createOrientations();
            obj.createOrientedMapping();
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
     
        function createOrientations(obj)
            coord = obj.mesh.coord;
            x1 = coord(:,1);
            x2 = coord(:,2);
            x10 = (max(x1(:))+min(x1(:)))/2;
            x20 = -0.5*max(x2(:));                                    
            r = sqrt((x1-x10).^2+(x2-x20).^2);
            fR = [x1-x10./r,x2-x20./r];
            fT = [-(x2-x20)./r,x1-x10./r];
            obj.orientation{1} = obj.createOrientationField(fR);
            obj.orientation{2} = obj.createOrientationField(fT);
        end 

        function f = createOrientationField(obj,fV)
            fD = obj.createP1DiscontinousOrientation(fV);
            f = obj.computeDoubleAngleOrientation(fD);            
            f  = obj.computeOppositeSignInLeftPart(f);
        end

        function fS = computeDoubleAngleOrientation(obj,f)
            aV = f.fValues;
            aV1 = aV(:,1);
            aV2 = aV(:,2);
            fN(:,1) = 2*aV1.^2-1;
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

        function fS = computeOppositeSignInLeftPart(obj,fS)
            fElem = fS.getFvaluesByElem;
            fSElem = fElem;
            isLeft = obj.isLeft;
            fSElem(:,:,isLeft) = -fElem(:,:,isLeft);           
            connec = fS.getDofConnec;
            for iNode = 1:3
                iDof = (iNode-1)*fS.ndimf+1;
                node = (connec(:,iDof)-1)/fS.ndimf+1;
                fV(node,:) = squeeze(fSElem(:,iNode,:))';
            end
            s.fValues = fV;
            s.mesh    = obj.mesh;
            s.order   = 'P1D';
            s.ndimf   = 2;
            fS = LagrangianFunction(s);
        end

        function itIs = isLeft(obj)
            coord= obj.mesh.computeBaricenter();
            x1 = coord(:,1);
            x2 = coord(:,2);
            itIs = x1 < (min(x1(:))+ max(x1(:)))/2;
        end
        
        function createOrientedMapping(obj)
            s.mesh = obj.mesh;
            s.orientation = obj.orientation;
            oM = OrientedMappingComputer(s);
            dCoord = oM.computeDeformedCoordinates()
        end
        
    end
    
end