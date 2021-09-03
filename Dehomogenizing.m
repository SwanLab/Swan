classdef Dehomogenizing < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        backgroundMesh
        boundaryMesh
        levelSet       
        uMesh
    end
    
    properties (Access = private)
        nx1
        nx2
        coord
        epsilon
        Ncell
    end
    
    methods (Access = public)
        
        function obj = Dehomogenizing()
            obj.init();
            obj.createCoord();            
            obj.createEpsilon();
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();
            obj.createLevelSet();            
            obj.createUnfittedMesh();
            figure()
            obj.uMesh.plotStructure();
            figure()
            obj.uMesh.plot();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nx1   = 252;
            obj.nx2   = 252;
            obj.Ncell = 5;
        end
        
        function createCoord(obj)
            x1 = linspace(0,1,obj.nx1);
            x2 = linspace(0,1,obj.nx2);
            x1T = repmat(x1,obj.nx2,1);
            x2T = repmat(x2',1,obj.nx1);
            obj.coord = [x1T(:),x2T(:)];
        end
        
        function createEpsilon(obj)
            x1max = max(obj.coord(:,1));
            x1min = min(obj.coord(:,1));
            x2max = max(obj.coord(:,2));
            x2min = min(obj.coord(:,2));            
            xMax = max(x1max,x2max);
            xMin = min(x1min,x2min);
            obj.epsilon = (xMax-xMin)/obj.Ncell;            
        end
        
        function createBackgroundMesh(obj)
            x1 = obj.coord(:,1);
            x2 = obj.coord(:,2);
            %[x1,x2] = obj.rotateCoordinates(x1,x2);                        
            connec   = delaunay(x1,x2);
            s.coord  = obj.coord;
            s.connec = connec;
            obj.backgroundMesh = Mesh(s);  
            x = linspace(min(x1),max(x1),obj.nx1);
            y = obj.computeMicroCoordinate(x);
            y = obj.periodicFunction(y);
           % figure(1)
           % hold on
           % plot(x,y)
        end
        
        function createBoundaryMesh(obj)
            sB.backgroundMesh = obj.backgroundMesh;
            sB.dimension = 1:3;
            sB.type = 'FromReactangularBox';
            bMc = BoundaryMeshCreator.create(sB);
            obj.boundaryMesh  = bMc.create();            
        end
        
        function createLevelSet(obj)
            [fx1,fx2] = obj.createMacroFunctions();
            [fy1,fy2] = obj.createMicroFunctions();     
            q = obj.createMacroExponent();            
            obj.levelSet = obj.createLevelSetValue(fx1,fx2,fy1,fy2,q);
        end
               
        function [yp1,yp2] = createMicroFunctions(obj)
            x1 = obj.coord(:,1);
            x2 = obj.coord(:,2);              
            [y1,y2] = obj.createFastCoordinates(x1,x2);    
            yp1 = obj.periodicFunction(y1);
            yp2 = obj.periodicFunction(y2);                        
            %[fy1,fy2] = obj.rotateCoordinates(fy1,fy2);            
        end
        
        function [y1,y2] = rotateCoordinates(obj,x1,x2)
            theta = obj.computeTheta(x1,x2);
            R(1,1,:) = cos(theta);
            R(1,2,:) = -sin(theta);
            R(2,1,:) = sin(theta);
            R(2,2,:) = cos(theta);
            y1 = (squeeze(R(1,1,:)).*x1 + squeeze(R(1,2,:)).*x2);
            y2 = (squeeze(R(2,1,:)).*x1 + squeeze(R(2,2,:)).*x2);           
        end
        
        function [y1,y2] = createFastCoordinates(obj,x1,x2)
            [x1,x2] = obj.rotateCoordinates(x1,x2);            
            y1 = obj.computeMicroCoordinate(x1);
            y2 = obj.computeMicroCoordinate(x2);
        end        
        
        function q = createMacroExponent(obj)
            x1 = obj.coord(:,1);
            x2 = obj.coord(:,2);                        
            q = obj.macroExponent(x1,x2);
        end
        
        function [fx1,fx2] = createMacroFunctions(obj)
            x1 = obj.coord(:,1);
            x2 = obj.coord(:,2);                        
            fx1 = obj.macroFunction(x1);
            fx2 = obj.macroFunction(x2);            
        end
        
        function createUnfittedMesh(obj)
            s.boundaryMesh   = obj.boundaryMesh;
            s.backgroundMesh = obj.backgroundMesh;
            obj.uMesh = UnfittedMesh(s);
            obj.uMesh.compute(obj.levelSet)
        end        
        
        function y = computeMicroCoordinate(obj,x)
            eps = obj.epsilon;        
            y = (x-min(x))/eps;
        end        
                   
    end
    
    methods (Access = private, Static)
                
        function phi = createLevelSetValue(fx1,fx2,fy1,fy2,q)
           %phi = (fy1.*fy2)./(fx1.*fx2)-1;            
           phi = 1-(abs((fy1-1/2)./fx1).^q + abs((fy2-1/2)./fx2).^q);            
        end
        
        function f = periodicFunction(y)
            %f = abs(cos(pi/2*y)).^2;
            f = y - floor(y);
        end        
        
        function f = computeTheta(x1,x2)
            x = x1;%max(x1,x2);            
            thetaMin = -1*pi/32;
            thetaMax = 0;
            xmin = min(x);
            xmax = max(x);
            theta = thetaMax +(thetaMin-thetaMax)*(x-xmin)/(xmax-xmin);
            f = theta;         
        end
        
        function f = macroExponent(x1,x2)
            x = max(x1,x2);            
            qMin = 36;
            qMax = 36;
            xmin = min(x);
            xmax = max(x);
            q = qMax +(qMin-qMax)*(x-xmin)/(xmax-xmin);
            f = q;
        end
       
        function f = macroFunction(x)
            %f = (1-0.5*x);
            mxmin = 0.905;
            mxmax = 0.205;
            xmin = min(x);
            xmax = max(x);
            m = mxmax +(mxmin-mxmax)*(x-xmin)/(xmax-xmin);
            f = (m)/2;
        end
        
    end
    
end