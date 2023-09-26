NTc = 100; %pocet vzorku na cip kodu
Np = NTc*1023; %poèet vzorku na periodu kodu
T = 1; % pocet period - integracni doba v ms

fv = Np * 1000;
t1 = 25;
% generovani signalu
s = sGPSGen(1,Np,T,0,0) + (___)*sGPSGen(1,Np,T,0,0+___/100/1023);
rp = sGPSGen(1,Np,T,0,0.5);	%maximum korelace v polovine Rcross
Rcross = korelace(s, rp);

d=___;
for k = Np/2-NTc:Np/2
    if Rcross(k)<=Rcross(k+d)
        te=k-Np/2+d/2;
    end
end
_____ % velikost chyby mereni pseudovzdalenosti
