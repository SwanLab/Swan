classdef StartingWithDehomogTwoVaria < handle

    properties (Access = private)
        mesh
        bField
        rhoField
    end

    methods (Access = public)

        function obj = StartingWithDehomogTwoVaria(bDesignVariable, rhoDesignVariable)
            obj.init()
            obj.mesh     = bDesignVariable.mesh;
            obj.bField   = bDesignVariable;
            obj.rhoField = rhoDesignVariable;
            obj.createLevelSet(0.035)
        end

    end

    methods (Access = private)

        function init(obj)
        end

        function ls = createLevelSet(obj, eps)
            close all

            x = obj.mesh.coordElem;
            fValues = obj.geometricalFunction(x, eps);
            
            s.mesh    = obj.mesh;
            s.fValues = fValues(:);
            s.order   = 'P1D';
            ls = LagrangianFunction(s);

            sF.trial = LagrangianFunction.create(obj.mesh, 1, 'P1');
            sF.mesh  = obj.mesh;
            filter   = FilterLump(sF);
            ls       = filter.compute(ls, 3);

            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(-ls.fValues);
            uMesh.plot()
            
            
            
           
        end

        % function fH = geometricalFunction(obj, xV, eps)
        %     s.fHandle = @(x) obj.coordFun(x);
        %     s.ndimf   = 2;
        %     s.mesh    = obj.mesh;
        %     xCoord    = AnalyticalFunction(s);
        % 
        %     s.operation = @(xV) obj.evaluateCellCoord(xV, eps, xCoord);
        %     s.mesh      = obj.mesh;
        %     s.ndimf     = xCoord.ndimf;
        %     txi         = DomainFunction(s);
        %     xiV         = txi.evaluate(xV);            
        % 
        % 
        %     b   = obj.bField.evaluate(xV);
        %     rho = obj.rhoField.evaluate(xV);
        % 
        % 
        %     a = exp(b.^2);
        %     d = (1 + b.^2) ./ a;
        % 
        % 
        %     xi1  = xiV(1,:,:);
        %     xi2  = xiV(2,:,:);
        %     xic1 = xi1 - 0.5;
        %     xic2 = xi2 - 0.5;
        % 
        % 
        %     phi1 = a.*xic1 + b.*xic2;
        %     phi2 = b.*xic1 + d.*xic2;
        % 
        %     % rhoMin = 0.1;
        %     % rhoMax = 0.7;
        %     % rhoEff = min(max(rho,rhoMin),rhoMax);
        %
        %     fH = max(abs(phi1)./(rho/2),abs(phi2)./(rho/2)) - 1;
        % end
        function fH = geometricalFunction(obj, xV, eps)
            s.fHandle   = @(x) obj.coordFun(x);
            s.ndimf     = 2;
            s.mesh      = obj.mesh;
            xCoord      = AnalyticalFunction(s);

            s.operation = @(xV) obj.evaluateCellCoord(xV, eps, xCoord);
            s.mesh      = obj.mesh;
            s.ndimf     = xCoord.ndimf;
            txi         = DomainFunction(s);
            xiV         = txi.evaluate(xV);

            b   = obj.bField.evaluate(xV);
            rho = obj.rhoField.evaluate(xV);

            % Projecta rho para 0/1 com threshold 0.5
            % Suaviza a transicao com sigmoid para evitar descontinuidades
            beta     = 50;                          % controla a nitidez
            rhoProj  = 1 ./ (1 + exp(-beta*(rho - 0.5)));
            rhoProj  = min(max(real(rhoProj), 0), 1 - 1e-6);

            a = exp(b.^2);
            d = (1 + b.^2) ./ a;

            xic1 = xiV(1,:,:) ;
            xic2 = xiV(2,:,:) ;

            e1 =  d.*xic1 - b.*xic2;
            e2 = -b.*xic1 + a.*xic2;

            holeSize = sqrt(1 - rhoProj) / 2;
            holeSize = max(holeSize, 1e-6);

            fH = max(abs(e1)./holeSize, abs(e2)./holeSize) - 1;
        end

        function f = coordFun(obj, x)
            x1 = x(1,:,:);
            x2 = x(2,:,:);
            f  = [x1; x2];
        end

        function txi = evaluateCellCoord(obj, xV, eps, xCoord)
            x   = xCoord.evaluate(xV);
            y   = obj.computeMicroCoordinate(x, eps);
            txi = obj.periodicFunction(y);
        end

        function y = computeMicroCoordinate(obj, x, eps)
            y = x / eps(1);
        end

    end

    methods (Access = private, Static)

        function f = periodicFunction(y)
            f = y - floor(y);
            f = f - 0.5;
        end

    end

end