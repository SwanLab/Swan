classdef SmoothingExponentComputerOptimal < SmoothingExponentComputer
    
   properties (Access = private)
       alpha
       beta
       gamma   
       m1
       m2
       qMin
       qMax
   end
    
   methods (Access = public)
       
       function obj = SmoothingExponentComputerOptimal(cParams)
           obj.init(cParams)
       end       
       
   end

   methods (Access = protected)
       
       function computeExponent(obj) 
            obj.alpha = 0.2505;
            obj.beta  = 0.5483;
            obj.gamma = 1;         
            obj.qMin = 2;
            obj.qMax = 32;           
            q = obj.computeQ();          
            obj.value = q;
       end
       
       function computeExponent2(obj) 
            obj.alpha = 6;
            obj.beta  = 20;
            obj.gamma = 4;  
            obj.qMin = 2;
            obj.qMax = 512;            
            a = obj.alpha;
            b = obj.beta;
            c = obj.gamma;           
            x = max(obj.m1,obj.m2);
            q = min(512,c*(1/(1-x^b))^a);  
            obj.value = q;
       end
       
       function computeExponent3(obj)
           fN = 'OptimalSuperEllipseExponent';
           pD = 'Topology Optimization/Vademecums';
           file2SaveName = [pD,'/',fN,'.mat'];
           s = load(file2SaveName,'d');  
           d = s.d;
           mx = d.mx;
           my = d.my;
           q  = d.q;
           if (obj.m1 > max(mx(:)) || obj.m2 > max(my(:)))
               q = 32;
           elseif (obj.m1 < min(mx(:)) || obj.m2 < min(my(:)))
               q = 2;                      
           else
             q = interp2(mx,my,q,obj.m1,obj.m2);
           end
           obj.value = q;
       end
       
   end   
   
   methods (Access = private)
       
       function init(obj,cParams)
          obj.m1 = cParams.m1;
          obj.m2 = cParams.m2;
       end
       
       function q = computeQ(obj)
            a = obj.alpha;
            b = obj.beta;
            c = obj.gamma;       
            qmin = obj.qMin;
            qmax = obj.qMax;
            x = max(obj.m1,obj.m2);
            q = max(qmin,min(qmax,c*(1/(1-x^b))^a));             
       end
       
   end   
 
end