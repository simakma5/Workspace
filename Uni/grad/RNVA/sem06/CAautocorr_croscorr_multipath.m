Np = 10*1023;
T = 1;

% sum
CN0dB = 45;
fv = Np * 1000;
Pn = fv/(10^(CN0dB/10));
n = (Pn^0.5)*(randn(1,Np*T)+1j*randn(1,Np*T))/2^0.5;

% generovani signalu
s = sGPSGen(1,Np,T,0,0) + 0.3*sGPSGen(1,Np,T,0,0+1.6/1023);
% s = sGPSGen(1,Np,T,0,0) + 0.3*sGPSGen(1,Np,T,0,0+1.6/1023) + sGPSGen(5,Np,T,0,0.8);
% s = sGPSGen(1,Np,T,0,0) + n;
%s = sGPSGen(1,Np,T,0,0) + 0.3*sGPSGen(1,Np,T,0,0+1.6/1023) + sGPSGen(5,Np,T,0,0.8) + n;
% s = sGPSGen(1,Np,T,0,1.5/1023);
PRNno1=1;

s_spekt = fft(s);
figure(1)
s_spekt(2000:length(s_spekt)-2000)=zeros(1,length(s_spekt)-3999);
plot(fftshift(abs(s_spekt)))
s = real(ifft(s_spekt));

figure(2)
    

for k = 1:9
    % generovani repliky
    rp = sGPSGen(k,Np,T,0,0.5);
    PRNno2= k;
    Rcross = xcorr([s, s], rp); %* (1/(Np*T);
    Rcross = (Rcross(2*Np:3*Np-1));
    
    zoom on;
    xaxis = (-(Np-1)/2:(Np-1)/2)*1023/Np;
    subplot(3,3,k);
    plot(xaxis, Rcross);
    %plot(Rcross);
    title(strcat('PRN ', num2str(PRNno1), ' and ', num2str(PRNno2)));
    xlabel('Shift [number of chips]');
    
end
figure(2)
    
