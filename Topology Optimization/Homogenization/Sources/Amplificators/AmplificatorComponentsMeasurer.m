classdef AmplificatorComponentsMeasurer < handle
    
    properties (Access = public)
       PmaxNorm
       Pl2Norm
       monomials
    end
    
    properties (Access = private)
        microCase
        pNorm
        vademecumData        
        Ptensor
        mxV
        myV
        
        chi
        mxT
        myT
        volume
        feasibleIndex
        
    end
    
    methods (Access = public)
        
        function obj = AmplificatorComponentsMeasurer(d)
            obj.init(d);
        end
        
        function compute(obj)
            obj.computeVademecumAmplificatorData();
            obj.obtainMonomials();
            obj.obtainDomainVariables();
            obj.obtainCellVariables();
            obj.computeAmplificatorTensorComponentsNorm(); 
            obj.computeComponentTxiRhoNorm();
            obj.plotComponent();
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj,d)
            obj.microCase = d.microCase;
            obj.pNorm     = d.pNorm;
        end
        
         function computeVademecumAmplificatorData(obj)
            d.fileName = obj.microCase;
            d.pNorm    = obj.pNorm;
            va = VademecumPtensorComputer(d);            
            va.compute();
            obj.vademecumData = va.vademecumData;
         end
        
         function obtainMonomials(obj)
             obj.monomials = obj.vademecumData.monomials;
         end
         
         function obtainDomainVariables(obj)
              obj.mxV = obj.vademecumData.domVariables.mxV;            
              obj.myV = obj.vademecumData.domVariables.myV;                         
         end
         
        function obtainCellVariables(obj)          
           for imx = 1:length(obj.mxV)
               for imy = 1:length(obj.myV)
                    P = obj.vademecumData.variables{imx,imy}.Ptensor;
                    v = obj.vademecumData.variables{imx,imy}.volume;
                    obj.Ptensor(imx,imy,:) = P;
                    obj.volume(imx,imy) = v;
               end
           end            
        end        
        
        function computeAmplificatorTensorComponentsNorm(obj)
            ncomp = size(obj.Ptensor,3);
            for icomp = 1:ncomp
                pcomp = obj.Ptensor(:,:,icomp);
                obj.PmaxNorm(icomp,1) = max(abs(pcomp(:)));
                obj.Pl2Norm(icomp,1) = obj.computeL2norm(pcomp);
            end
        end
        
        function l2norm = computeL2norm(obj,fval)
            l2norm = sqrt(trapz(obj.myV,trapz(obj.mxV,fval.*fval,2)));
        end        
        
        function computeComponentTxiRhoNorm(obj)
            obj.computeTxiMxMyVariables();
            obj.computeFeasibleIndex();                        
        end
        
        function computeTxiMxMyVariables(obj)
            for i = 1:length(obj.mxV)
                for j = 1:length(obj.myV)
                    mx = obj.mxV(i);
                    my = obj.myV(j);
                    obj.chi(i,j) = mx/my;
                    obj.mxT(i,j) = mx;
                    obj.myT(i,j) = my;
                end
            end
        end
                  
        function computeFeasibleIndex(obj)
            d.mx = obj.mxV;
            d.my = obj.myV;
            d.chi = obj.chi;
            d.rho = obj.volume;
            fC = FeasibleIndexComputer(d);
            obj.feasibleIndex = fC.index;
        end        
        
        function plotComponent(obj)
            ind = obj.feasibleIndex;
            x = obj.chi(ind);
            y = obj.volume(ind);
            val = squeeze(obj.Ptensor(:,:,4));
            z = val(ind);
            ncolors = 50;
            tri = delaunay(x,y);
            f = figure;
            tricontour(tri,x,y,z,ncolors)
            colorbar
            hold on
            plot(x,y,'+');
            xlabel('$\frac{m1}{m2}$','Interpreter','latex');
            ylabel('\rho');                        
        end
        
    end
    
end