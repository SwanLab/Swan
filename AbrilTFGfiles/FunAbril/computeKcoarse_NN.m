function K=computeKcoarse_NN(K_NN,r,nf)
    K_aux1=K_NN.computeOutputValues(r);
    K_aux2=zeros(nf);
    idx=1;
    for n=1:nf
        for m=n:nf
            K_aux2(n,m)=K_aux1(idx);
            idx=idx+1;
        end
    end
    K=K_aux2+triu(K_aux2,1).';
end