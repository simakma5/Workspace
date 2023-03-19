function Pxx = PSD_ARMA_first(un,B,fs,wlen,type1,type2)
%
% use as ARMA_first(....,wlen,'MA','NF')
 
Ts = 1/fs;
mod = exp(-B*pi*Ts);
if type1=='AR' %AR
    b=1;
    if type2=='NF' % NF
        a = [1 -mod];
    else %VF
        a = [1 mod];
    end
else %MA
    a=1;
    if type2=='NF' % NF
        b = [1 mod];
    else %VF
        b = [1 -mod];
    end
end

nc1 = filter(b,a,un);
Pxx = pwelch(nc1,wlen,wlen/2,wlen,fs,'twosided','psd');
Pxx_dB = 10*log10(Pxx);
disp(Pxx_dB(1:10))
end




