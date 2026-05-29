classdef StartingWithDehomogSquare < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        mesh
    end

    methods (Access = public)

        function obj = StartingWithDehomogSquare()
            obj.init()
            obj.mesh = UnitTriangleMesh(100,100);            
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

         function fH = geometricalFunction(obj,xV,eps)
            
            s.fHandle = @(x) obj.coordFun(x);
            s.ndimf   = 2;
            s.mesh    = obj.mesh;
            xCoord     = AnalyticalFunction(s);   


            s.operation = @(xV) obj.evaluateCellCoord(xV,eps,xCoord);
            s.mesh      = obj.mesh;
            s.ndimf     = xCoord.ndimf;
            txi = DomainFunction(s);

            txiV = txi.evaluate(xV);
            
            x1_macro = xV(1,:,:);
            x2_macro = xV(2,:,:);

            % m(x) constante por agora — depois virá da optimização
            m_x1 = 0.5 * ones(size(x1_macro));
            m_x2 = 0.5 * ones(size(x2_macro));

            s.m_x1        = m_x1;
            s.m_x2        = m_x2;
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.5;
            s.type        = 'SquareVariableSize';



            g = GeometricalFunction(s);
            f = g.getHandle;
            fH = f(txiV);

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