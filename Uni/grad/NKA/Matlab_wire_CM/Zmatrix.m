function [Z,delta] = Zmatrix(f,L,N,a,d)

% f.. frekvence, L.. delka dipolu, N.. segmentu, a.. polomer
% d.. vzdalenost dipolu pro pripad vypoctu vzajemne vazby
% dvou kolinearnich dipolu
% vraci [Z] a delta (delka jednoho segmentu)

c=3e+8;
lambda  = c/f;                 % wavelength
mi      = 4*pi*1e-7;         % permeability of vacuum
epsilon = 8.85e-12;          % permittivity of vacuum

% N number of segments
delta   = L/(N+1);  		     % length of one segment
k	     = 2*pi/lambda;        % wavenumber
omega   = k*3e+8;            % angular frequency
alpha   = 0.5*delta;         % one half of the segment's length

psi = zeros( 1, N+1);        % numerical integration
for m = 1:(N+1)
    x = (m-1)*delta;
%  psi(m) = quadl( 'green', -alpha, alpha, 1e-5, [], x, a, k);
    psi(m) = quad( 'green', -alpha, alpha,  1e-5, [], x, a, k, d);
end

psi = psi/delta;
%SelfTerm pro d=0 (samotny dipol)
if d==0;
% d2a = delta/(2*a);
% psiself = log(d2a+sqrt(1+d2a^2))/(2*pi)-(j*k*delta)/(4*pi);
% psi(1) = psiself/delta;
% psiself/delta


d4a = delta/(4*a);
F = @(x)log(d4a+((d4a)^2+(sin(x/2)).^2).^0.5);
Q = (1/(4*pi^2))*quad(F,0,2*pi);
T = (-i*k*delta/(4*pi))+(log(2)/(2*pi))+Q;
psi(1) = T/delta;


end;

Z = zeros( N, N);            % computing impedance matrix
multi = 1j*omega*mi;
divid = 1j*omega*epsilon;
% for m = 1:N
%     for n = m:N
%         dist   = abs(m-n);
%         hlp    = 2*psi(1+dist) - psi(1+abs(dist-1)) - psi(1+abs(dist+1));
%         Z(m,n) = multi*(delta^2)*psi(1+dist) + hlp/divid;
%         Z(n,m) = Z(m,n);
%     end
% end
% 
for m = 1:N
    for n = 1:N
        dist   = abs(m-n);
        hlp    = 2*psi(1+dist) - psi(1+abs(dist-1)) - psi(1+abs(dist+1));
        Z(m,n) = multi*(delta^2)*psi(1+dist) + hlp/divid;
    end
end




