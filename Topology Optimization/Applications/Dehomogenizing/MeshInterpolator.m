classdef MeshInterpolator < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
       fineMesh 
       coarseMesh
       coarseMeshDisc
       funToRemesh
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = MeshInterpolator()
            close all
            obj.createCoarseMesh();
            obj.createCoarseDiscontinousMesh();            
            obj.createP1ContinousFunction();
            obj.remeshP1ContinousFunction();
            obj.createP1DiscontinousFunction();
            obj.remeshP1DiscontinousFunction();
        end
        
    end
    
    methods (Access = private)
               




        






        function remeshP1ContinousFunction(obj)
            m  = obj.coarseMesh;
            mF = m.remesh();        
            fFineC = fC.refine(mCoarse,mFine);   
            fFineD.plot(mFine)            
        end
        
    end
    
end