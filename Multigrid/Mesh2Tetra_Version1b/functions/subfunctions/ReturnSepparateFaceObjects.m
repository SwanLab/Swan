function [nObjects Object]=ReturnSepparateFaceObjects(F)
F_new=zeros(size(F));  nF=size(F,1);
N=zeros(1,numel(F));

Object=cell(1,100);
nObjects=0;
while(nF>0)
    nN=1; N(nN)=F(1,1);
    nF_new=0;
    while(nN>0)
        row=find((F(1:nF,1)==N(1))|(F(1:nF,2)==N(1))|(F(1:nF,3)==N(1)),1,'first');
        if(~isempty(row))
            nF_new=nF_new+1;
            if(F(row,1)==N(1))
                nN=nN+1; N(nN)=F(row,2);
                nN=nN+1; N(nN)=F(row,3);
            elseif(F(row,2)==N(1))
                nN=nN+1; N(nN)=F(row,1);
                nN=nN+1; N(nN)=F(row,3);
            else
                nN=nN+1; N(nN)=F(row,1);
                nN=nN+1; N(nN)=F(row,2);
            end
            
            F_new(nF_new,:)=F(row,:);
            F(row,:)=F(nF,:); nF=nF-1;    
        else
            N(1)=N(nN); nN=nN-1;
        end
    end 
    nObjects=nObjects+1;
    Object{nObjects}=F_new(1:nF_new,:);
end
