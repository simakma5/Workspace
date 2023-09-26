NTc = 100; %pocet vzorku na cip kodu
Np = NTc*1023; %po�et vzorku na periodu kodu
T = 1; % pocet period - integracni doba v ms
fv = Np * 1000;

s0 = sGPSGen(1,Np,T,0,0);
rp = sGPSGen(1,Np,T,0,0.5);
for t1 = 1:150
	% generovani signalu
    s1 = (0.1)*sGPSGen(1,Np,T,0,0+t1/100/1023);
	
	RcrossP = korelace(s0 + s1, rp);
    RcrossM = korelace(s0 - s1, rp);
        
	d1 = NTc;
	d2 = NTc/10;
	for k = Np/2-NTc:Np/2
        if RcrossP(k)<=RcrossP(k+d1)
            te1P=k-Np/2+d1/2;
        end
        if RcrossP(k)<=RcrossP(k+d2)
            te2P=k-Np/2+d2/2;
        end
        if RcrossM(k)<=RcrossM(k+d1)
            te1M=k-Np/2+d1/2;
        end
        if RcrossM(k)<=RcrossM(k+d2)
            te2M=k-Np/2+d2/2;
        end
	end
	teg1P(t1) = te1P/fv*3e8; % velikost chyby mereni pseudovzdalenosti
    teg2P(t1) = te2P/fv*3e8;
    teg1M(t1) = te1M/fv*3e8;
    teg2M(t1) = te2M/fv*3e8;
    t1
end
t1 = 1:150;
t1 = t1/100;
figure(2)
plot(t1,teg1P,t1,teg2P,t1,teg1M,t1,teg2M);
xlabel('delay [Tchip]');
ylabel('error [m]');

x = [3 9 25 50 75 95 100 105 125 145 150];
zobraz = [x/100; teg1P(x); teg2P(x); teg2M(x); teg1M(x)] 