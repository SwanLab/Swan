function [P,Pder ]= CoordinateMatrixPolynomialTRI(COOR,ORDER_POLY)

if nargin == 0
    load('tmp1.mat')
elseif nargin  ==1
    ORDER_POLY = 1; 
end


if nargout == 1
    COMPUTE_DERIVATIVE = 0 ;
else
    COMPUTE_DERIVATIVE = 1;
end

[nnode,ndim] = size(COOR) ;


if   ORDER_POLY == 1 
    
   
     P = [ones(nnode,1),COOR] ; 
     Pder = cell(2,1) ; 
     Pder{1} = zeros(nnode,3) ; 
      Pder{1}(:,2) = 1 ; 
        
  Pder{2} = zeros(nnode,3) ; 
      Pder{2}(:,3) = 1 ; 
    
else
    
    error('Option not implemented yet')
    
    
end


nnode  = size(COOR,1) ;
% 
% if ndim  == 1
%     
%     iacum = 0;
%     for ix=0:n(1)
%      %   for iy=0:n(2)
%             iacum = iacum+1;
%             P(:,iacum) = (COOR(:,1).^ix)   ; ;
%       %  end
%     end
%     Pder = {} ;
%     if COMPUTE_DERIVATIVE == 1
%         % For the derivative with respect to  x
%         iacum = 0;
%         Pderx = zeros(nnode,nmon) ;
%         for ix=0:n(1)
%           %  for iy=0:n(2)
%                 iacum = iacum+1;
%                 if ix >=1
%                     Pderx(:,iacum) = ix*(COOR(:,1).^(ix-1))  ;
%                 end
%            % end
%         end
%         
% %         % For the derivative with respect to  y
% %         iacum = 0;
% %         Pdery = zeros(nnode,nmon) ;
% %         for ix=0:n(1)
% %             for iy=0:n(2)
% %                 iacum = iacum+1;
% %                 if iy >=1
% %                     Pdery(:,iacum) = (COOR(:,1).^(ix)).*(COOR(:,2).^(iy-1))*iy  ;
% %                 end
% %             end
% %         end
%         
%         Pder = {Pderx} ;
%     end
%     
% elseif ndim == 2
%     iacum = 0;
%     for ix=0:n(1)
%         for iy=0:n(2)
%             iacum = iacum+1;
%             P(:,iacum) = (COOR(:,1).^ix).*(COOR(:,2).^iy)  ; ;
%         end
%     end
%     Pder = {} ;
%     if COMPUTE_DERIVATIVE == 1
%         % For the derivative with respect to  x
%         iacum = 0;
%         Pderx = zeros(nnode,nmon) ;
%         for ix=0:n(1)
%             for iy=0:n(2)
%                 iacum = iacum+1;
%                 if ix >=1
%                     Pderx(:,iacum) = ix*(COOR(:,1).^(ix-1)).*(COOR(:,2).^iy)  ;
%                 end
%             end
%         end
%         
%         % For the derivative with respect to  y
%         iacum = 0;
%         Pdery = zeros(nnode,nmon) ;
%         for ix=0:n(1)
%             for iy=0:n(2)
%                 iacum = iacum+1;
%                 if iy >=1
%                     Pdery(:,iacum) = (COOR(:,1).^(ix)).*(COOR(:,2).^(iy-1))*iy  ;
%                 end
%             end
%         end
%         
%         Pder = {Pderx,Pdery} ;
%     end
% elseif ndim ==3
%     % 3D -----------------------------------
%     
%     
%     iacum = 0;
%     for ix=0:n(1)
%         for iy=0:n(2)
%             for iz = 0:n(3)
%                 iacum = iacum+1;
%                 P(:,iacum) = (COOR(:,1).^ix).*(COOR(:,2).^iy).*(COOR(:,3).^iz)   ;
%             end
%         end
%     end
%     Pder = {} ;
%     if COMPUTE_DERIVATIVE == 1
%         % For the derivative with respect to  x
%         iacum = 0;
%         Pderx = zeros(nnode,nmon) ;
%         for ix=0:n(1)
%             for iy=0:n(2)
%                 for iz = 0:n(3)
%                     iacum = iacum+1;
%                     if ix >=1
%                         Pderx(:,iacum) = ix*(COOR(:,1).^(ix-1)).*(COOR(:,2).^iy).*(COOR(:,3).^iz)  ;
%                     end
%                 end
%             end
%         end
%         
%         % For the derivative with respect to  y
%         iacum = 0;
%         Pdery = zeros(nnode,nmon) ;
%         for ix=0:n(1)
%             for iy=0:n(2)
%                 for iz = 0:n(3)
%                     iacum = iacum+1;
%                     if iy >=1
%                         Pdery(:,iacum) = (COOR(:,1).^(ix)).*((COOR(:,2).^(iy-1))*iy).*(COOR(:,3).^iz)   ;
%                     end
%                 end
%             end
%         end
%         
%         % For the derivative with respect to  z
%         iacum = 0;
%         Pderz = zeros(nnode,nmon) ;
%         for ix=0:n(1)
%             for iy=0:n(2)
%                 for iz = 0:n(3)
%                     iacum = iacum+1;
%                     if iz >=1
%                         Pderz(:,iacum) = (COOR(:,1).^(ix)).*((COOR(:,3).^(iz-1))*iz).*(COOR(:,2).^iy)   ;
%                     end
%                 end
%             end
%         end
%         
%         Pder = {Pderx,Pdery,Pderz} ;
%     end
%     
% else
%     error('Option not implemented')
%     
% end
% 
% 
% 
