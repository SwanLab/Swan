function K=computeKcoarse_NN(K_NN,r)
    K_aux1=K_NN.computeOutputValues(r);
    K_aux2=zeros(8);
    idx=1;
    for n=1:8
        for m=n:8
            K_aux2(n,m)=K_aux1(idx);
            idx=idx+1;
        end
    end
    K=K_aux2+triu(K_aux2,1).';
end