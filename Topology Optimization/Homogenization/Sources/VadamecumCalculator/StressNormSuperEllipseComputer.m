classdef StressNormSuperEllipseComputer < handle
    
    properties (Access = private)
        microProblem
        microProblem2
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
        hMesh
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
            %   obj.createMesh();
            %   obj.createMicroProblem();
            obj.createMicroProblem2()
            sPnorm = obj.computePstressNorm();
        end
        
        function var = computeCellVariables(obj)
            obj.compute();
            mProblem = obj.microProblem2;
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
            microP = obj.microProblem2;
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
                microP = obj.microProblem2;
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
            obj.hMesh    = cParams.hMesh;
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
            %sPnorm2 = obj.microProblem.computeStressPnorm(stress,p);
            sPnorm2 = obj.microProblem2.computeStressPnorm(stress,p);
        end
        
        function createMicroProblem(obj)
            d.gmsFile = [fullfile(obj.outputFolder,obj.fName),'.msh'];
            d.outFile = obj.fName;
            d.print   = obj.print;
            d.iter = obj.iter;
            d.hasToCaptureImage = obj.hasToCaptureImage;
            nH = NumericalHomogenizerCreatorFromGmsFile(d);
            homog = nH.getHomogenizer();
            mProblem = homog.getMicroProblem();
            
            mN = mProblem.mesh;
            coord = mN.coord;
            coord(:,3) = 0;
            mshG = msh(coord,mN.connec);
            refine  = mmg(mshG,1e-3);
            hsiz(refine,0.007);
            nosurf(refine);
            [mesh1] = run(refine);
            
            s.coord = mesh1.vtx(:,1:2);
            s.connec = mesh1.elt;
            
            mN1 = Mesh().create(s);
            
            % mN1.plot()
            drawnow
            
            mProblem.setMesh(mN1);
            props.kappa = .75;
            props.mu    = .375;
            mProblem.setMatProps(props);
            
            obj.microProblem = mProblem;
            
        end
        
        function createMesh(obj)
            d = SettingsFreeFemMeshGenerator();
            d.freeFemFileName = 'SmoothRectangle';
            d.hMax  = obj.hMesh;%0.002;%0.0025;
            d.mxV             = obj.mx;
            d.myV             = obj.my;
            d.fileName        = obj.fileName;
            d.printingDir     = obj.outputFolder;
            d.qNorm           = obj.q;
            fG = FreeFemMeshGenerator(d);
            fG.generate();
        end
        
        function [coord, connec] = readCoordConnec(obj)
            obj.testName = 'RVE_Square_Triangle_FineFine';
            %obj.testName = 'RVE_Square_Triangle_Fine';
            run(obj.testName)
        end
        
        function createMicroProblem2(obj)
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
            hminBmesh = bMesh.computeMinCellSize();
            hmeanBmesh = bMesh.computeMeanCellSize();
            %hmin(meshMmg,hminBmesh);
            % hmax(meshMmg,10*hminBmesh);
            hausd(meshMmg,hminBmesh)
            %nosurf(meshMmg);
            
            map(meshMmg,ls.value);
            %hsiz(refine,0.007);
            %nosurf(meshMmg);
            %hgrad(refine,10)
            
            verbose(meshMmg,-1);
            [mesh1] = runLs(meshMmg);
            
            %             plot(mesh1)
            %
            
            
            it = mesh1.col == 3;
            connec = mesh1.elt(it,:);
            coord  = mesh1.vtx(:,1:2);
            [s.coord,s.connec] = obj.computeUniqueCoordConnec(coord,connec);
            
            
            mN2 = Mesh().create(s);
            %  mN2.plot()
            
            
            %
            %
            %     mN = Mesh().create(sM);
            %             mN2.computeMinCellSize
            %             mN2.computeMeanCellSize
            %  mN.plot();
            
            
            %
            %             coord = mN.coord;
            %             coord(:,3) = 0;
            %             mshG = msh(coord,mN.connec);
            %             meshMmg  = mmg(mshG,1e-3);
            %             hsiz(meshMmg,0.007);
            %             nosurf(meshMmg);
            %             verbose(meshMmg,-2);
            %             [mesh1] = run(meshMmg);
            %
            %             coord = mesh1.vtx(:,1:2);
            %             connec = mesh1.elt;
            %
            %             [s.coord,s.connec] = obj.computeUniqueCoordConnec(coord,connec);
            %
            %             mN1 = Mesh().create(s);
            %             mN1.plot()
            
            %mF = mN1;
            
            
            mF = mN2;
            
           figure(10)
           clf()
           mF.plot();
           drawnow
            
            femSolver = Elastic_Problem_Micro.create(obj.testName);
            femSolver.setMesh(mF);
            
            props.kappa = .75;
            props.mu    = .375;
            femSolver.setMatProps(props);
            obj.microProblem2 = femSolver;
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