classdef testFourthOrderAmplificatorTensor < testShowingError
    
    properties (Access = protected)
        tol = 5e-2;
        testName = 'testShapeStressWithAmplificator';
    end
    
    properties (Access = private)
        fileName   
        printingDir
        gmsFile   
        homogenizer
        mx
        my
        mxV
        myV
        gPtensor
        ncomp
        alpha
        PmaxNorm
        Pl2Norm        
    end
    
    
    methods (Access = public)
        
        function obj = testFourthOrderAmplificatorTensor()
            obj.init();
            obj.computeVademecumAmplificatorTensor();
            obj.computeAmplificatorTensorComponentsNorm();
            obj.plotNorms();
        end
        
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            obj.error = 0;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName    = 'FourthOrderAmplificator';
            obj.printingDir = fullfile(pwd,'Output',obj.fileName);  
            obj.gmsFile     = [fullfile(obj.printingDir,obj.fileName),'.msh'];  
            obj.ncomp = 21;          
            nMx = 20;
            nMy = 20;
            obj.mxV = linspace(0.01,0.99,nMx);
            obj.myV = linspace(0.01,0.99,nMy);     
            obj.gPtensor = zeros(nMx,nMy,obj.ncomp);
            obj.computeMultiIndex();
        end
        
        function computeMultiIndex(obj)
            p = 2;
            n = 6;
            obj.alpha = multinomial_expand(p,n);              
        end
        
        function computeVademecumAmplificatorTensor(obj)           
            for imx = 1:length(obj.mxV)
                for imy = 1:length(obj.myV)
                    obj.mx = obj.mxV(imx);
                    obj.my = obj.myV(imy);   
                    obj.generateMeshFile();
                    obj.createFemMatOoInputData();
                    obj.createNumericalHomogenizer();
                    obj.obtainGeneralizedPtensor(imx,imy);
                end                
            end           
        end
        
        function generateMeshFile(obj)
            d.mxV         = obj.mx;
            d.myV         = obj.my;
            d.fileName    = obj.fileName;
            d.printingDir = obj.printingDir;
            d.freeFemFileName = 'SmoothRectangle';
            fG = FreeFemMeshGenerator(d);
            fG.generate();
        end       
        
        function createFemMatOoInputData(obj)
            oD = fullfile('Output',obj.fileName);
            fullOutFile = fullfile(oD,[obj.fileName,'.m']);
            c = GmsFile2FemMatOoFileConverter(obj.gmsFile,oD,fullOutFile);
            c.convert();
        end           
        
        function createNumericalHomogenizer(obj)
            defaultDB = NumericalHomogenizerDataBase([obj.fileName,'.m']);
            dB = defaultDB.dataBase;
            dB.outFileName                   = obj.fileName;
            dB.print                         = true;
            dB.levelSetDataBase.levelSetType = 'full';
            obj.homogenizer = NumericalHomogenizer(dB);                                            
        end
        
        function obtainGeneralizedPtensor(obj,imx,imy)
            obj.gPtensor(imx,imy,:) = obj.homogenizer.generalizedPtensor();
        end
        
        function computeAmplificatorTensorComponentsNorm(obj)
            for icomp = 1:obj.ncomp
                pcomp = obj.gPtensor(:,:,icomp);
                obj.PmaxNorm(icomp) = max(abs(pcomp(:)));
                obj.Pl2Norm(icomp) = obj.computeL2norm(pcomp);
            end            
        end
        
        function l2norm = computeL2norm(obj,fval)
            l2norm = sqrt(trapz(obj.myV,trapz(obj.mxV,fval.*fval,2)));                
        end
        
        function plotNorms(obj)
            obj.plotNorm(obj.PmaxNorm,'PtensorCompLogInfNorm');
            obj.plotNorm(obj.Pl2Norm,'PtensorCompLogL2Norm');
        end
        
        function plotNorm(obj,norm,plotName)
            fig = figure();
            h = bar(log(norm));
            alphaLeg = obj.computeAlphaLegend();
            set(gca, 'XTickLabel',alphaLeg, 'XTick',1:numel(alphaLeg));
            set(gca,'XTickLabelRotation',45);
            p = barPrinter(fig,h);
            path = '/home/alex/Dropbox/Amplificators/Images/FourthOrderAmplificator/';
            p.print([path,plotName]);
        end
        
        function a = computeAlphaLegend(obj)
            for ia = 1:size(obj.alpha,1)
               a{ia} = mat2str(obj.alpha(ia,:)); 
            end
        end
        
    end
    
    
end