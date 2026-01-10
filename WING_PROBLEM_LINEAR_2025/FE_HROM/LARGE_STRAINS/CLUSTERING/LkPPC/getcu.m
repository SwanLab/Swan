function [cluster]=getcu(a)   
cluster=[];m=length(a); B=1;i=1;
for j=1:m
    if all(j*ones(size(B))-B) || j==1 
        test=j;Y=[];
            while (~isempty(find(a(test,:)==1)))  
                [Y,a,N]=f2(a,test,Y);  
                test=unique(N);
               
            end
            cluster{i}=unique(Y);
            i=i+1;B=[B Y];B=unique(B);
    end
end
end
function [Y,a,N]=f2(a,test,Y)
N=[];  
for i=1:length(test)
    [~,n]=find(a(test(i),:)==1);
    if ~isempty(n)   
        for j=1:length(n)
            a(test(i),n(j))=0;
            a(n(j),test(i))=0;
        end
    end
    N=unique([N n]);
end
b=unique([N test]);
Y=unique([Y b]);
end