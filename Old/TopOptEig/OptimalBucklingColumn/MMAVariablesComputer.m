classdef MMAVariablesComputer < handle

    properties (Access = private)
        designVariable
        xMin
        xMax
        xOld1
        xOld2
        lOW
        uPP
        a0Val
        aMMA
        dVal
        cVal
        x0
    end
    
    properties (Access = private) 
        mesh
        nConstraints
        nValues
        minThick
        maxThick
    end
    
    methods (Access = public)
        
        function obj = MMAVariablesComputer(cParams)
            obj.init(cParams)
            obj.computeInitialVariablesMMA();
        end

        function xmma = compute(obj,nValues,nConstraints,iter,xval,f0val,df0dx,df0dx2,fval,dfdx,dfdx2)
            n_val = nValues;
            m = nConstraints;
            xmin = obj.xMin;
            xmax = obj.xMax;
            xold1 = obj.xOld1;
            xold2 = obj.xOld2;
            low = obj.lOW;
            upp = obj.uPP;
            a0 = obj.a0Val;
            a_mma = obj.aMMA;
            d = obj.dVal;
            c = obj.cVal;            
            [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
                mmasub(m,n_val,iter,xval,xmin,xmax,xold1,xold2, ...
                f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a_mma,c,d);
            obj.lOW = low;
            obj.uPP = upp;
            obj.xOld1 = xold1;

            obj.xOld2 = obj.xOld1;
            obj.xOld1 = xval;            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nConstraints = cParams.nConstraints;           
            obj.minThick      = cParams.minThick;
            obj.maxThick      = cParams.maxThick;
            obj.nValues      =  cParams.nValues;
            obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
        end

        function obj = computeInitialVariablesMMA(obj)
            x = obj.designVariable.value;
            N = obj.mesh.nelem;
            n_val = obj.nValues;
            m = obj.nConstraints;
            alpha = obj.minThick;
            beta = obj.maxThick;
            obj.xMin=alpha*ones(N,1);
            obj.xMin=[obj.xMin; 0];
            obj.xMax=beta*ones(N,1);
            obj.xMax=[obj.xMax; 1000];
            obj.xOld1=x;
            obj.xOld2=x;
            obj.lOW = zeros(n_val,1);
            obj.uPP = ones(n_val,1);
            obj.a0Val = 1;
            obj.aMMA = zeros(m,1);
            obj.dVal = zeros(m,1);
            obj.cVal = 1000*ones(m,1);        
        end        

    end
    
end