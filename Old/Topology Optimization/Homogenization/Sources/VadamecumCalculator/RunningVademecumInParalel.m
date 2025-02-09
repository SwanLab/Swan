classdef RunningVademecumInParalel < handle

    properties (Access = private)
        samplePoints
        cellVar
        fileName
    end

    properties (Access = private)
        nCase
        nMx
        nMy
        nPhi
        phiMin
        phiMax
    end
    
    methods (Access = public)
        
        function obj = RunningVademecumInParalel(cParams)
            obj.init(cParams);    
            obj.compute();
        end
        
    end
    

     methods (Access = private)

        function init(obj,nC)
            obj.nCase  = nC;
            obj.nMx    = 20;
            obj.nMy    = 20;
            obj.nPhi   = 10;
        end         
        
        function compute(obj)
            obj.computeSameStressSign();
            obj.computeDifferentStressSign();                        
        end
        
        function computeSameStressSign(obj)
            obj.phiMin = 0;
            obj.phiMax = pi/4;      
            obj.fileName = 'OptimalSuperEllipseSameStressSign';
            obj.computeVademecums();
        end
        
        function computeDifferentStressSign(obj)
            obj.phiMin = pi/2;
            obj.phiMax = 3*pi/4;      
            obj.fileName = 'OptimalSuperEllipseDifferentStressSign';
            obj.computeVademecums();                     
        end
        
        function computeVademecums(obj)
            obj.createSamplePoints();                 
            obj.cellVar = cell(obj.nMx,obj.nMy,obj.nPhi);
            imx = obj.nCase;
            for imy = 1:obj.nMy
                obj.displayPercentatge(imx,imy);
                for iphi = 1:obj.nPhi
                    c = obj.computeOneVademecum(imx,imy,iphi);
                    obj.cellVar{imx,imy,iphi} = c;
                end
            end
            c = obj.cellVar;
            save([obj.fileName,num2str(imx)],'c');            
        end
        
        
        function cellV = computeOneVademecum(obj,imx,imy,iphi)
            iter = obj.computeGlobalIteration(imx,imy);
            s.rho = obj.samplePoints.rhoV(iter);
            s.xi = obj.samplePoints.txiV(iter);
            s.phi = obj.samplePoints.phiV(iphi);   
            s.nCase = [num2str(imx),num2str(imy),num2str(iphi)];
            v = VademecumCalculator(s);
            v.compute();            
            %v.cellVariables = s.nCase;
            cellV = v.cellVariables;
        end
         
       function createSamplePoints(obj)
            s.type = 'FromMxMy';
            s.nMx  = obj.nMx;
            s.nMy  = obj.nMy;
            s.nPhi = obj.nPhi;
            s.phiMin = obj.phiMin;
            s.phiMax = obj.phiMax;
            sample = SamplePointsCreatorForOptimalExponentComputer.create(s);
            sample.compute();
            obj.samplePoints = sample;
        end        
        
        function iter = computeGlobalIteration(obj,imx,imy)
            iter = (imy + obj.nMx*(imx -1));
        end

        function displayPercentatge(obj,imx,imy)
            perc = imy/(obj.nMy)*100;
            %disp([num2str(perc),'% done']);
            fid = fopen('Output.txt','a+');
            fprintf(fid,'%6.2f %%\n',perc);
            fclose(fid);
        end       
        
     end
     
end