classdef StartingWithDehomogb < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        mesh
        bField
    end

    methods (Access = public)

        function obj = StartingWithDehomogb(bDesignVariable)
            obj.init()
            obj.mesh   = bDesignVariable.fun.mesh;
            obj.bField = bDesignVariable.fun; 
            obj.createLevelSet(0.2)
        end

    end

    methods (Access = private)

        function init(obj)

        end
        

        function ls = createLevelSet(obj,eps)
            %       s.operation  = @(xV) obj.geometricalFunction(xV,eps);
            %        s.ndimf      = 1;
            %        s.mesh       = obj.fineMesh;
            %        f  = DomainFunction(s);
            %        ls = project(f,'P1');
            close all

            x = obj.mesh.coordElem;

            fValues   = obj.geometricalFunction(x,eps);
            s.mesh    = obj.mesh;
            s.fValues = fValues(:);
            s.order   = 'P1D';
            ls = LagrangianFunction(s);
            ls = project(ls,'P1');

            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls.fValues);
            uMesh.plot()
        end


        function fH = geometricalFunction(obj, xV, eps)
            s.fHandle = @(x) obj.coordFun(x);
            s.ndimf = 2;
            s.mesh = obj.mesh;
            xCoord = AnalyticalFunction(s);
            s.operation = @(xV) obj.evaluateCellCoord(xV, eps, xCoord);
            s.mesh = obj.mesh;
            s.ndimf = xCoord.ndimf;
            txi = DomainFunction(s);
            xiV = txi.evaluate(xV);

            
            b = obj.bField.evaluate(xV);
            
            

            
            a = exp(b.^2);
            d = (1 + b.^2) ./ a;          
            
            xi1 = xiV(1,:,:);
            xi2 = xiV(2,:,:);

            
            phi1 = a .* xi1 + b .* xi2;
            phi2 = b .* xi1 + d .* xi2;
            phiV = [phi1; phi2];

            


            s_type.type        = 'SquareDeformedLS';
            s_type.m_x1        = 0.5 * ones(size(b));   
            s_type.m_x2        = 0.5 * ones(size(b));   
            s_type.xCoorCenter = 0.5;
            s_type.yCoorCenter = 0.5;
            g  = GeometricalFunction(s_type);
            f  = g.getHandle;
            fH = f(phiV);
        end

        function f = coordFun(obj,x)
            x1 = x(1,:,:);
            x2 = x(2,:,:);
            f = [x1;x2];
        end

        function txi = evaluateCellCoord(obj,xV,eps,xCoord)
            x = xCoord.evaluate(xV);

            %  x = obj.deformedCoord.evaluate(xV);
            y = obj.computeMicroCoordinate(x,eps);
            txi = obj.periodicFunction(y);
        end


        function y = computeMicroCoordinate(obj,x,eps)
            % nDim = size(x,1);
            % y = zeros(size(x,1),size(x,2),obj.mesh.nelem);
            % for iDim = 1:nDim
            %     xI    = x(iDim,:,:);
            %     %xImin = min(xI(:))
            %    % xImin = eps(iDim)
            %    % xImin = 0;
            %     y(iDim,:,:) = (xI)/(eps(iDim));
            % end
            y = x/eps(1);
            %y = (x-min(x(:))-eps)/eps;
        end  
    end

    methods (Access = private, Static)

        function f = periodicFunction(y)
            %f = ((cos(1*pi*(y)))).^2;           
        %    f = (cos(2*pi*y));   

            %f = abs(y-floor(y)-0.5);
            f = abs(y-floor(y)-0.5);
        end        

    end

end