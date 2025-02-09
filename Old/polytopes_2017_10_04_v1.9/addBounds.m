function [A,b,Aeq,beq]=addBounds(A0,b0,Aeq0,beq0,lb,ub)
%addBounds - given a sextuplet {A, b, Aeq, beq, lb, ub} of matrices and vectors 
%of the kind used by the Optimization Toolbox to express linear
%constraints, the routine will convert the sextuplet to a quadruplet {A, b,
%Aeq, beq}. 
%
%      [A,b,Aeq,beq]=addBounds(A,b,Aeq,beq,lb,ub)
%
%In other words, it will re-write the bounds given by vectors lb, ub as equalities
%and inequalities C*x<=c, D*x=d and append these as new rows to the matrices
%A, b, Aeq, beq.
%
%
%EXAMPLE:
%
%   >> [A,b,Aeq,beq]=addBounds([1 1 1],1,[],[],[0;0;1],[1;1;1])
% 
%         A =
% 
%              1     1     1
%             -1     0     0
%              0    -1     0
%              1     0     0
%              0     1     0
% 
% 
%         b =
% 
%              1
%              0
%              0
%              1
%              1
% 
% 
%         Aeq =
% 
%              0     0     1
% 
% 
%         beq =
% 
%              1




  %%%%begin parsing
  
   if nargin<5, error 'At least 5 arguments required'; end
   
   if size(A0,1)~=length(b0)
       error 'Incompatible inequality matrix data sizes: size(A,1) ~= length(b)'
   end
       
   if size(Aeq0,1)~=length(beq0)
       error 'Incompatible equality matrix data sizes: size(Aeq,1) ~= length(beq)'
   end



     if ~exist('lb','var'), lb=[]; end
      if ~exist('ub','var'), ub=[]; end      
    
      if ~isempty(lb)
         N=length(lb);
      elseif ~isempty(ub);
         N=length(ub);
      else
          error 'Either lb or ub must be nonempty.'
      end
    
      if isempty(lb)
        lb=-inf(N,1); 
      end

      if isempty(ub)
        ub=+inf(N,1); 
      end
      
      if length(ub)~=length(lb)
         error 'Arguments lb and ub must be either [] or the same lengths' 
      end
      
     
      lb=lb(:); ub=ub(:);

      %%%end parsing
      
      
      lisinf=~isfinite(lb);
      uisinf=~isfinite(ub); 
      idx_eq=(ub==lb)& ~lisinf & ~uisinf;

      Au=[eye(N),ub];
      Al=[-eye(N),-lb];      
      Aeq=Au(idx_eq,:);
      Au(idx_eq|uisinf,:)=[];
      Al(idx_eq|lisinf,:)=[];
      
      A=[Al;Au]; b=A(:,end); A(:,end)=[];
                 beq=Aeq(:,end); Aeq(:,end)=[];
      
                 
      A=[A0;A];
      b=[b0;b];
      Aeq=[Aeq0;Aeq];
      beq=[beq0;beq];
      
      if any(lb==inf | ub==-inf)
         A(end+1,end)=0; b(end+1)=-1; 
      end
      
                 