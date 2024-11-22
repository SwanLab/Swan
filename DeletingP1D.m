classdef DeletingP1D < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
       mesh
       fun
       funP1
    end
    
    methods (Access = public)
        
        function obj = DeletingP1D()
            obj.init()
            obj.createMesh();
            obj.createAnalyticalFunciton();
            obj.createP1ContinousFun();
            obj.createP1DiscontinousFun();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end
        
        function createMesh(obj)
            m = TriangleMesh(1,1,20,20);
            obj.mesh = m;
        end
        
        function createAnalyticalFunciton(obj)
            s.fHandle = @(x) cos(100*x(1,:,:)+x(2,:,:));
            s.ndimf  = 1;
            s.mesh   = obj.mesh;
            obj.fun = AnalyticalFunction(s);
        end

        function createP1ContinousFun(obj)
            obj.funP1 = obj.fun.project('P1');
        end

        function createP1DiscontinousFun(obj)
            funP1D = obj.funP1.project('P1D');

            s.ndimf = 1;
            s.order = 'P1';
            s.mesh      = obj.mesh;
            s.fValues   = funP1D.fValues;
            s.dofConnec = funP1D.dofConnec;
            s.dofCoord  = funP1D.dofCoord;
            s.dofs.getNumberDofs = size(funP1D.dofCoord,1);
            funP1DC = LagrangianFunction(s);


            obj.funP1.computeL2norm()
            funP1DC.computeL2norm()

            obj.funP1.plot()
            funP1DC.plot()

            gradP1DC = Grad(funP1DC).project('P1',obj.mesh);
            gradP1DC.plot()

            gradP1 = Grad(obj.funP1).project('P1',obj.mesh);
            gradP1.plot()

            dif = gradP1DC - gradP1;
            dif.computeL2norm()/gradP1.computeL2norm()
        end        
        
        
        
        
        
    end
    
end