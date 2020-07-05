classdef OneOptimalExponentComputerAndFunctionVariation < handle
    
    properties (Access = public)
        qOptIter
        fOptIter
        fValues
        qValues
        qMin
        qMax        
    end
    
    properties (Access = private)
        stressProblem        
        print
        hasToCaptureImage        
    end
    
    properties (Access = private)
        rho
        txi
        phi
        pNorm
        fileName
        hMesh
    end
    
    methods (Access = public)
        
        function obj = OneOptimalExponentComputerAndFunctionVariation(cParams)
            obj.init(cParams);
            obj.createStressProblem();
            obj.computeQbounds();               
        end
                
        function compute(obj)         
            obj.computeOptimalExponent();
            %obj.printOptimalMicroStructure();
            nValues = 50;
            obj.computeStressNormRelationWithQ(nValues);
        end
                
        function printOptimalMicroStructure(obj)
            obj.print = true;
            obj.hasToCaptureImage = true;
            obj.createStressProblem();
            qOpt = obj.qOptIter(end);
            obj.stressProblem.objective(qOpt);
        end
        
        function f = computeValue(obj,q)
             f = obj.stressProblem.objective(q);
        end
        
        function computeStressNormRelationWithQ(obj,nValues)
            obj.print = true;
            obj.createStressProblem();
            obj.computeQbounds();                 
            
            qVal = linspace(obj.qMin,obj.qMax,nValues);
            fVal = zeros(nValues,1);
            for ivalues = 1:nValues
                obj.stressProblem.iter = ivalues;
                qV = qVal(ivalues);
                fVal(ivalues) = obj.stressProblem.objective(qV);
                disp([num2str(ivalues/nValues*100),'%'])
            end
            obj.qValues = qVal;
            obj.fValues = fVal;
        end        
        
        function computeOptimalExponent(obj)
            f = @(q) obj.stressProblem.objective(q);
            xIter = [];
            fIter = [];
            options = optimset('Display','iter','OutputFcn', @myoutput);
            options.MaxIter = 200;
            options.TolX = 1e-5;
            x = fminbnd(f,obj.qMin,obj.qMax,options);
            function stop = myoutput(x,optimvalues,state)
                stop = false;
                if isequal(state,'iter')
                    xIter = [xIter; x];
                    fIter = [fIter;optimvalues.fval];
                end
            end
            obj.qOptIter = xIter;
                         obj.fOptIter = fIter;
        end
        
        function cVariables = obtainCellVariables(obj,q)
            v = obj.stressProblem.computeCellVariables(q);
            cVariables = v;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
            obj.rho   = cParams.rho;
            obj.txi   = cParams.txi;
            obj.phi   = cParams.phi;
            obj.pNorm = cParams.pNorm;
            obj.hMesh = cParams.hMesh;
            obj.print = false;
            obj.hasToCaptureImage = false;
        end
        
        function computeQbounds(obj)
            obj.qMin = obj.stressProblem.qMin;
            obj.qMax = obj.stressProblem.qMax;
        end
        
        function createStressProblem(obj)
            s.rho                  = obj.rho;
            s.xi                   = obj.txi;
            s.fileName             = obj.fileName;
            s.phi                  = obj.phi;
            s.hMesh                = obj.hMesh;
            s.pNorm                = obj.pNorm;
            s.print                = obj.print;
            s.hasToCaptureImage = obj.hasToCaptureImage;
            sProblem   = StressNormVsQproblemCreator(s);
            sProblem.compute();
            obj.stressProblem = sProblem;
        end           
        
    end
    
end


