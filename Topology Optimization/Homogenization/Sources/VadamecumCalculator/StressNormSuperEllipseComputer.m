classdef StressNormSuperEllipseComputer < handle
    
    properties (Access = private)
        microProblem
        fName
        outputFolder
    end
    
    properties (Access = private)
        mx
        my
        q
        phi
        pNorm
        print
        fileName
        hasToCaptureImage
        testName
        iter
    end
    
    methods (Access = public)
        
        function obj = StressNormSuperEllipseComputer(cParams)
            obj.init(cParams);
        end
        
        function sPnorm = compute(obj)
            obj.createMicroProblem()
            sPnorm = obj.computePstressNorm();
        end
        
        function var = computeCellVariables(obj)
            obj.compute();
            mProblem = obj.microProblem;
            var.Ctensor = mProblem.variables.varFromCh.Chomog;
            var.tstress = mProblem.variables.varFromCh.tstress;
            var.tstrain = mProblem.variables.varFromCh.tstrain;
            var.displ   = mProblem.variables.varFromCh.tdisp;    
            var.volume  = mProblem.mesh.computeVolume();
            var.mesh    = mProblem.mesh;
            
            quad = Quadrature.set(var.mesh.geometryType);
            quad.computeQuadrature('CONSTANT');
            volume = var.mesh.computeDvolume(quad);  
            var.dV    = volume;            
            var.nstre = size(var.tstress,1);
            var.ngaus = size(var.tstress,2);
        end
        
        function printImage(obj)
            microP = obj.microProblem;
            outputName = [obj.fileName,'Print'];
            if obj.print
                outName = outputName;
                dI.mesh            = microP.mesh;
                dI.outName         = outName;
                dI.pdim            = '2D';
                dI.ptype           = 'MICRO';
                ps = PostProcessDataBaseCreator(dI);
                dB = ps.getValue();
                dB.printers = 'Density';
                postCase = 'Density';
                postProcess = Postprocess(postCase,dB);
                d.x = ones(size(microP.mesh.coord(:,1),1),1);
                it = 0;
                postProcess.print(it,d);
                if obj.hasToCaptureImage
                    f = outputName;
                    imagePath = '/home/alex/git-repos/MicroStructurePaper/';
                    outPutNameImage = fullfile(imagePath,[f,num2str(it)]);
                    inputFileName = fullfile('Output',f,[f,num2str(it),'.flavia.res']);
                    sI.fileName = f;
                    sI.outPutImageName = outPutNameImage;
                    sI.inputFileName = inputFileName;
                    imageCapturer = GiDImageCapturer(sI);
                    imageCapturer.capture();
                end
            end
        end
        
        function printStress(obj)
            if obj.print
                microP = obj.microProblem;
                outputName = [obj.fileName,'Print'];
                dI.mesh    =  microP.mesh;
                dI.outName = outputName;
                dI.pdim    = '2D';
                dI.ptype   = 'MICRO';
                ps = PostProcessDataBaseCreator(dI);
                dB = ps.getValue();
                postCase = 'ElasticityMicro';
                postProcess = Postprocess(postCase,dB);
                d.fields = microP.variables;
                d.quad   = microP.element.quadrature;
                postProcess.print(obj.iter,d);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mx       = cParams.mx;
            obj.my       = cParams.my;
            obj.q        = cParams.q;
            obj.phi      = cParams.phi;
            obj.pNorm    = cParams.pNorm;
            obj.print    = cParams.print;
            obj.fileName = cParams.fileName;
            obj.iter     = cParams.iter;
            obj.hasToCaptureImage = cParams.hasToCaptureImage;
            obj.computeFileNameAndOutputFolder();
        end
        
        function computeFileNameAndOutputFolder(obj)
            obj.fName = [obj.fileName];
            obj.outputFolder = fullfile(pwd,'Output',obj.fName);
        end
        
        function sPnorm2 = computePstressNorm(obj)
            stress = [sin(obj.phi) cos(obj.phi) 0];
            p = obj.pNorm;
            sPnorm2 = obj.microProblem.computeStressPnorm(stress,p);
        end
        
        function [coord, connec] = readCoordConnec(obj)
            obj.testName = 'RVE_Square_Triangle_FineFine';
            %obj.testName = 'RVE_Square_Triangle_Fine';
            run(obj.testName)
        end
        
        function createMicroProblem(obj)
            [coord, connec] = obj.readCoordConnec();
            s.coord = coord(:,2:end-1);
            s.connec = connec(:,2:end);
            mesh = Mesh_Total(s);
            
            
            s.widthH = obj.mx;
            s.widthV = obj.my;
            s.pnorm = obj.q;
            s.type = 'smoothRectangle';
            
            s.levelSetCreatorSettings = s;
            s.type = 'LevelSet';
            s.mesh = mesh;
            s.scalarProductSettings.epsilon = 1;
            
            ls = LevelSet(s);
            
            uMesh = ls.getUnfittedMesh;
            
            
            cInner = uMesh.meshBackground.connec(uMesh.backgroundFullCells,:);
            connec = [cInner;uMesh.innerCutMesh.connec];
            coord = uMesh.innerCutMesh.coord;
            
            [sM.coord,sM.connec] = obj.computeUniqueCoordConnec(coord,connec);
                        
            bMesh = uMesh.meshBackground;
            % bMesh.plot()
            
            
            coord = bMesh.coord;
            coord(:,3) = 0;
            mshG = msh(coord,bMesh.connec);
            meshMmg  = mmg(mshG,1e-3);
           
            hausd(meshMmg,0.01);            
            hmin(meshMmg,0.001);
            hmax(meshMmg,0.01);
            
            %hausd(meshMmg,0.005);            
            %hmin(meshMmg,0.0005);
            %hmax(meshMmg,0.005);            
            
            map(meshMmg,ls.value);
           
            verbose(meshMmg,-1);
            meshMmg.oldFileName = [obj.fileName,'Out'];
            meshMmg.newFileName = [obj.fileName,'In'];
            
            mesh1 = runLs(meshMmg);
            
            it = mesh1.col == 3;
            connec = mesh1.elt(it,:);
            coord  = mesh1.vtx(:,1:2);
            [s.coord,s.connec] = obj.computeUniqueCoordConnec(coord,connec);
            
            
            mF = Mesh().create(s);
              
           figure(10)
           clf()
           mF.plot();
           drawnow
            
            femSolver = Elastic_Problem_Micro.create(obj.testName);
            femSolver.setMesh(mF);
            
            props.kappa = .75;
            props.mu    = .375;
            femSolver.setMatProps(props);
            obj.microProblem = femSolver;
        end
        
        function [newCoord,newConnec] = computeUniqueCoordConnec(obj,coord,connec)
            allNodes = connec(:);
            [uNodes,ind,ind2] = unique(allNodes,'rows','stable');
            
            allCoords    = coord;
            uniqueCoords = allCoords(uNodes,:);
            newCoord    = uniqueCoords;
            
            nnode = size(connec,2);
            nCell = size(connec,1);
            newConnec = reshape(ind2,nCell,nnode);
        end
        
        
    end
    
    
end