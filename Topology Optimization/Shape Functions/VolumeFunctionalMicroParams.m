classdef VolumeFunctionalMicroParams < handle

    properties (Access = private)
        quadrature
        totalVolume
    end

    properties (Access = private)
        fileName
        microParams
    end

    methods (Access = public)
        function obj = VolumeFunctionalMicroParams(cParams)
            obj.init(cParams);
            obj.loadVademecum();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            J  = obj.computeFunction(xD);
            dJ = obj.computeGradient(xD);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
           obj.fileName    = cParams.fileName;
           obj.microParams = cParams.microParams;
        end

         function [mxV,myV,rho] = loadVademecum(obj)
            fName = [obj.fileName,'WithAmplificators'];
            matFile   = [fName,'.mat'];
            file2load = fullfile('Vademecums',matFile);
            v = load(file2load);
            var = v.d;
            mxV = var.domVariables.mxV;
            myV = var.domVariables.myV;
             for imx = 1:length(mxV)
                 for imy = 1:length(myV)
                     rho(imx,imy) = var.variables{imx,imy}.('Density');
                 end
             end
         end 

        function [mL,cells] = obtainLocalCoord(obj,xV)
            mx = obj.microParams{1};
            my = obj.microParams{2};
            mxG = mx.evaluate(xV);
            myG = my.evaluate(xV);
            mG(:,1) = mxG(:);
            mG(:,2) = myG(:);
            [mL,cells] = obj.sMesh.obtainLocalFromGlobalCoord(mG);
        end         

        function createTotalVolume(obj)
            dV = obj.mesh.computeDvolume(obj.quadrature);
            obj.totalVolume = sum(dV(:));
        end

        function J = computeFunction(obj,x)
            volume = Integrator.compute(x,obj.mesh,obj.quadrature.order);
            J      = volume/obj.totalVolume;
        end

        function dJ = computeGradient(obj,x)
            test    = obj.gradientTest;
            fValues = ones(test.nDofs,1)/obj.totalVolume;
            dJ      = FeFunction.create(test.order,fValues,obj.mesh);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume';
        end
    end
end

