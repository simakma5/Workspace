clear all
NTc = 1; %pocet vzorku na cip kodu
Np = NTc*1023; %pocet vzorku na periodu kodu
T = 1; % pocet period - integracni doba v ms

fv = Np * 1000;

% generovani signalu
PRN=1;
fd_rx=150;
tau_rx=0.5;
secnd_path_delay=160/1023;
% s = sGPSGen(PRN,Np,T,fd_rx,tau_rx);
s = sGPSGen(PRN,Np,T,fd_rx,tau_rx) + 0.5*sGPSGen(1,Np,T,fd_rx,tau_rx+secnd_path_delay);
% load('s_rec10101.mat');
% s = s_rec;
% krok_fd = 100;
krok_fd = 1000/(2*T); % [Hz] = 1/(2*doba integrace)
for fd = -10000:krok_fd:10000
    rp = sGPSGen(1,Np,T,fd,0);
    Rcross = korelace(s, rp);
    Rcross3D((fd+10000)/krok_fd+1,:) = Rcross;
end

% accumulation through T periods
Rcross3D_accum=zeros(2*10000/krok_fd+1,1023);
for fd = -10000:krok_fd:10000
    for td = 1:1023
      for taccum=0:T-1
          Rcross3D_accum((fd+10000)/krok_fd+1,td) = Rcross3D_accum((fd+10000)/krok_fd+1,td) + Rcross3D((fd+10000)/krok_fd+1,td+1023*(taccum));
      end
    end
end

figure(1)
surf(abs(Rcross3D_accum))
figure(2)
surf(abs(Rcross3D))
figure(3)
zoom on;
plot(-10000:krok_fd:10000,abs(Rcross3D_accum(:,round(tau_rx*Np)+1)));