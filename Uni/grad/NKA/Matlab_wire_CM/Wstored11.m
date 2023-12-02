function [We11,Wm11,Pr11] = Wstored11(J,N,delta,f,a)

%%
% a..polomer vodice
% d..vzdalenost mezi vodici

om=2*pi*f;
c=3e8;
epsilon = 8.85e-12;
k=om/c;

Jd=diff(J)./(delta); Jd(N)=Jd(N-1);
Jdconj=diff(conj(J))./(delta); Jdconj(N)=Jdconj(N-1);

a = a*(1-0.40976*(a/delta));

Wetosum = zeros(N,N);
Wmtosum = zeros(N,N);
WradG2tosum = zeros(N,N);
Prtosum = zeros(N,N);

for m = 1:N;
    for n = 1:N;
        dist = sqrt((abs(m-n)*delta)^2+a^2);
        if m==n;
        %%%% SELFTERMY
        sinkaovera = k*delta^2-(k^3*a^2*delta^2)/6;
        coskaovera = 2*a-2*sqrt(delta^2+a^2)-2*delta*log(a)+2*delta*log(delta+sqrt(delta^2+a^2));
        sinkr      = k*(((2*a^3)/3)+((delta^2+a^2)^(3/2)/3)-a^2*sqrt(delta^2+a^2)-...
                     a^2*delta*log(a)+a^2*delta*log(delta+sqrt(delta^2+a^2)));
        Prtosum(m,n)=(k^2*(J(m)*conj(J(n)))-(Jd(m)*Jdconj(n)))*...
                    sinkaovera;
        Wetosum(m,n)=((Jd(m)*Jdconj(n)))*...
                    coskaovera;
        Wmtosum(m,n)=(J(m)*conj(J(n)))*...
                    coskaovera;
        WradG2tosum(m,n)=(k^2*(J(m)*conj(J(n)))-(Jd(m)*Jdconj(n)))*...
                   sinkr;
        else
        Wetosum(m,n)=delta^2*((Jd(m)*Jdconj(n)))*...
                    (cos(k*dist)/dist);
        Wmtosum(m,n)=delta^2*(J(m)*conj(J(n)))*...
                    (cos(k*dist)/dist);
        WradG2tosum(m,n)=delta^2*((k^2*(J(m)*conj(J(n)))-(Jd(m)*Jdconj(n)))*...
                    (sin(k*dist)));
        Prtosum(m,n)=delta^2*(k^2*(J(m)*conj(J(n)))-(Jd(m)*Jdconj(n)))*...
                    (sin(k*dist)/dist);
        end;
    end;
end;

WradG2=sum(sum(WradG2tosum));
mpk=(16*pi*om^2*epsilon)^-1;
We11=real((mpk*(sum(sum(Wetosum))-(k/2)*WradG2)));
Wm11=real((mpk*(k^2*sum(sum(Wmtosum))-(k/2)*WradG2)));
Pr11=real(((1/(8*pi*om*epsilon))*sum(sum(Prtosum))));
