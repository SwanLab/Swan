function [ ak] = rec2(alpha,beta,a0,k)

    function[res] = aux(k,accu1,accu2)
        if k == 0
            res = accu1.*a0 + accu2;
        else
            alphak = alpha(k);
            betak = beta(k);
            accu2 = accu2 + accu1.*betak;
            accu1 = accu1.*alphak;
            res = aux(k-1,accu1,accu2);
        end
    end
ak = aux(k,1,0);
end

