function A = adjacency(options,stX)     %求邻接矩阵的函数；            
 n=size(stX,1);
 p=2:(options.NN+1);
    if size(stX,1)<500 % block size: 500
      step=n;
      else
	    step=500;
         end
       idy=zeros(n*options.NN,1);
         DI=zeros(n*options.NN,1); t=0; s=1;
        for i1=1:step:n    
                t=t+1;i2=i1+step-1;
           if (i2>n) 
               i2=n;
           end
             Xblock=stX(i1:i2,:);  
             dt=feval(options.GraphDistanceFunction,Xblock,stX);
             [Z,I]=sort(dt,2); Z=Z(:,p)'; I=I(:,p)'; [g1,g2]=size(I);
             idy(s:s+g1*g2-1)=I(:);DI(s:s+g1*g2-1)=Z(:);s=s+g1*g2;
        end 
       I=repmat((1:n),[options.NN 1]); I=I(:);
       t=mean(DI(DI~=0)); 
     A=sparse(I,idy,ones(prod(size(idy)),1),n,n);
end
function D = euclidean(A,B)   
    if (size(A,2) ~= size(B,2))
    error('A and B must be of same dimensionality.');
    end
    if (size(A,2) == 1) 
    A = [A, zeros(size(A,1),1)]; B = [B, zeros(size(B,1),1)];
    end
    a = dot(A,A,2);b = dot(B,B,2);ab=A*B';
    D = real(sqrt(bsxfun(@plus,a,b')-2*ab));
end