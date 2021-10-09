classdef MatlabtestPerimeterCircleTriangle < testUnfitted
    
    properties (Access = protected)
        testName = 'test_circle_triangle';
        analyticalValue = 2*pi;
        meshType = 'BOUNDARY';
    end

    %% Heredat de testUnfittedPerimeterIntegration
    properties (Access = protected)
        meshIncludeBoxContour = false
    end
    
    methods (Access = public)

        
        function obj = MatlabtestPerimeterCircleTriangle()
            obj.createTopOpt();
            obj.integrateSurface();
        end
    end
    %% Heredat de testUnfittedIntegration_ExternalIntegrator
    methods (Access = protected)
        
        function totalIntegral = computeGeometricalVariable(obj)
            switch obj.meshType  
                case 'INTERIOR'
                    totalIntegral = obj.unfittedMesh.computeMass();
                case 'BOUNDARY'
                    totalIntegral = obj.unfittedMesh.computePerimeter();
            end

        end
        
    end
    
    %% Heredat de testShowingError
   
    properties (Access = protected)
       error
    end
        
    methods (Access = protected)
        function printTestPassed(obj)
           cprintf('green',obj.FileName);                                    
           cprintf('green',' PASSED.');
           cprintf('black',['Error: ',num2str(obj.error),'\n']);
        end
        
        function printTestNotPassed(obj)
            cprintf('red',obj.FileName);                        
            cprintf('red',' FAILED.');
            cprintf('red',['Error: ',num2str(obj.error),'\n']);
        end
        
        function hasPassed = hasPassed(obj)
            obj.computeError()
            hasPassed = obj.error < obj.tol();
        end
    end

    %% Heredat de testUnfittedIntegration
    properties (Access = private)
        varAdim
    end
    properties (Access = protected)
        tol = 6e-2;
    end
    methods (Access = protected)
        
        function integrateSurface(obj)
            obj.createMesh();
            
            geomVar = obj.computeGeometricalVariable();
            obj.varAdim = geomVar/obj.analyticalValue;
        end
        
        function computeError(obj)
            obj.error = abs(obj.varAdim - 1);
        end
        
    end
end

