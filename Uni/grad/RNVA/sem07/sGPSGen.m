function s = sGPSGen(PRNno, Npc, Np, fd, tau);
% PRNno - Cislo PRN kodu
% Npc - pocet vzorku na periodu kodu
% Np - pocet period
% fd - doppleruv kmitocet v Hz
% tau - zpozdeni kodu v ms
PRN = CAcode(PRNno);
k = 1:Npc*Np;
sc = PRN.PRN11(floor(mod((k/Npc-tau)*1023,1023))+1);
fn = exp(2*pi*1j*k*fd/1000/Npc);
s = sc.*fn;