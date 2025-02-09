function [A,b,Aeq,beq,lb,ub,report] = separateBounds( A,b,Aeq,beq,tol )
%separateBounds - given a quadruplet {A, b, Aeq, beq} of linear constraint 
%matrices and vectors corresponding to the region A*x<=b, Aeq*x=beq, the routine
%will convert to a sextuplet {A, b, Aeq, beq, lb, ub} of the kind used
%in the Optimization Toolbox.
%
%      [A,b,Aeq,beq,lb,ub,report] = separateBounds( A,b,Aeq,beq )
%
%In other words, the code will look through the rows of A and Aeq for
%constraints corresponding to simple upper and lower bounds. It will then
%separate these constraints from the rest, expressing them instead using
%vectors lb, ub. The "report" output argument is a structure containing
%some relevnt stats,
%
%    report.infeasibleBounds: indices, i, of any lb(i)>ub(i)
%    report.inequalitiesRemoved: indices of rows removed from A,b
%    report.equalitiesRemoved: indices of rows removed from Aeq,beq
%
%EXAMPLE:
%
%     A =  [1     1     1
%          -1     0     0
%           0    -1     0
%           1     0     0
%           0     1     0] ;
%
%      b=[1 0 0 1 1].';
%
%     Aeq=[0 0 1]; 
%     beq=1;
%
%     >>[A,b,Aeq,beq,lb,ub] = separateBounds( A,b,Aeq,beq )
% 
%         A =
% 
%              1     1     1
% 
% 
%         b =
% 
%              1
% 
% 
%         Aeq =
% 
%            Empty matrix: 0-by-3
% 
% 
%         beq =
% 
%            Empty matrix: 1-by-0
% 
% 
%         lb =
% 
%              0
%              0
%              1
% 
% 
%         ub =
% 
%              1
%              1
%              1
%
%By default, the code will consider a constraint A(i,:)<=b(i) to correspond
%to a pure bound if a row A(i,:) or similarly Aeq(i,:) contains no more
%than one non-zero element. In certain situations, however, the matrices A
%and Aeq are the output of floating point calculations and such rows will
%contain non-zero elements representing numerical noise. In such instances,
%one can call the code with a tolerance parameter
%
%      [A,b,Aeq,beq,lb,ub,report] = separateBounds( A,b,Aeq,beq,tol )
%
%where 0<=tol<=1. The criterion for deciding whether a constraint is a pure
%bound is then,
%
%        max(abs(A(:,i)))>=tol*sum(abs(A(i,:)))
%
%The default beavior corresponds to tol=1.


  %%%%begin parsing

          if ~exist('Aeq','var'), Aeq=[]; end
          if ~exist('beq','var'),  beq=[]; end      
          if ~exist('tol','var')||isempty(tol), 
              tol=1; 
          elseif tol<0 || tol>1
              error 'Tolerance parameter must satisfy 0<=tol<=1';
          end 

          Na=nan; Nb=nan;
          if xor(isempty(A),isempty(b)) 
            error 'If A is empty then b must also be empty and vice versa'
          elseif ~isempty(A)
              Na=size(A,2);
          elseif isempty(A)
              A=[];b=[];
          end

          if xor(isempty(Aeq),isempty(beq))
              error 'If Aeq is empty then beq must also be empty and vice versa'
          elseif ~isempty(Aeq)
              Nb=size(Aeq,2);
          elseif isempty(Aeq)
                Aeq=[]; beq=[];
          end

          if ~isnan(Na) && ~isnan(Nb) 

              if Na~=Nb
               error 'If both A and Aeq are both non-empty, they must have same number of columns'
              end


          elseif isnan(Na) && isnan(Nb)

              lb=[]; ub=[];
              return
 

          end

          N=max(Na,Nb);
          

          b=b(:);  beq=beq(:); %henceforth we are sure b-data are column vectors
  
  
           if size(A,1)~=length(b)
               error 'Incompatible inequality matrix data sizes: size(A,1) ~= length(b)'
           end

           if size(Aeq,1)~=length(beq)
               error 'Incompatible equality matrix data sizes: size(Aeq,1) ~= length(beq)'
           end
          
          
  %%%%end parsing
  
  
  
   [ls1,lv1,us1,uv1,A1,b1,rows1] = extract(A,b,tol);
   [ls2,lv2,us2,uv2,A2,b2,rows2] = extract(Aeq,beq,tol);
   [ls3,lv3,us3,uv3] = deal(us2,uv2,ls2,lv2); 
  
   A=A1; 
   b=b1;
   Aeq=A2;
   beq=b2;  
   lb=accumarray([ls1;ls2;ls3]  ,  [lv1;lv2;lv3]  , [N,1], @max , -inf);
   ub=accumarray([us1;us2;us3]  ,  [uv1;uv2;uv3]  , [N,1], @min , +inf);
  
   report.infeasibleBounds=find(ub<lb).';
   report.inequalitiesRemoved=rows1.';
   report.equalitiesRemoved=rows2.';
   
   if nargout<7 &&  ~isempty(report.infeasibleBounds)
      
       warnstr=[sprintf('Infeasible combinations of bounds ( LB(i)>UB(i) ) were detected.\n\n'),...
               sprintf( 'To disable this warning, call with 7 output argument syntax:\n\n'),...
               sprintf( '          [A,b,Aeq,beq,lb,ub,report] = separateBounds(...)')];
       
       warning('separateBounds:infeasible', warnstr);
       
   end
   
   
 function [lsubs,lvals,usubs,uvals,A,b,rows] = extract(A,b,tol)

      lsubs=[]; usubs=[];
      lvals=[];  uvals=[];
      rows=[];
     
     if isempty(A), return; end
     
    b=b(:);
  
    absA=abs(A);
  
    [maxA,idx]=max(absA,[],2);
    sumA=sum(absA,2);
    
    rows =find(  maxA>=tol*sumA & sumA>0 );
  
   
    subs=idx(rows);
    
      ind=sub2ind(size(A),rows, subs);
      num=b(rows);
      den=A(ind);
    
    vals=num./den;
    
    lidx=(den<0);
    uidx=(den>0);
    
    lsubs=subs(lidx);
    lvals=vals(lidx);
   
    
    usubs=subs(uidx);
    uvals=vals(uidx);
    

   zrows=(sumA==0);
   rowsInfeas=find(zrows & b<0);
   rowsFeas=find(zrows & b>=0);
   
   if ~isempty(rowsInfeas) %row of all zeros
   
       N=size(A,2);
       lsubs=(1:N).';
       lvals=inf(N,1);
       usubs=lsubs;
       uvals=-lvals;
       rows=[rows;rowsInfeas];
   end
   
   if ~isempty(rowsFeas) %row of all zeros
       rows=[rows;rowsFeas];
   end   
   
    A(rows,:)=[];
    b(rows)=[];
   
   
   
