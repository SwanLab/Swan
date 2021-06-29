classdef OldVademecumReader < handle
    
   properties (Access = public)
      rhoV
      xiV
      phiV
      qV
   end
   
   properties (Access = private)
      path
      fileName      
   end
    
   methods (Access = public)
       
       function obj = OldVademecumReader()
           obj.init();
           obj.read();
       end
       
   end
   
   methods (Access = private)
       
       function init(obj)
           obj.path = 'Topology Optimization/Vademecums/';
           obj.fileName = 'OptimalSuperEllipseExponentData';           
       end
       
       function read(obj)
           nPhi = 10;
            for i = 1:nPhi
                fName = [obj.path,obj.fileName,num2str(i),'.mat'];
                d = load(fName);
                q(:,i) = d.x.q(i,:);
            end      
           nPoints = size(q,1);            
           obj.qV = q;
           obj.rhoV = repmat(d.x.rho',1,nPhi);
           obj.xiV  = repmat(d.x.txi',1,nPhi);d.x.txi;
           obj.phiV = repmat(d.x.psi,nPoints,1);
       end
   end
    
    
    
    
end