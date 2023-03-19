% Code dispersion.m	Courtesy: Debesh Bhatta
% the software describes pulse propagation in free space in one-dimension
% we study the effect of numerical dispersion on pulse shape by changing the value of alpha.
% provision for Gaussian and rectangular pulse shapes
% instability of the solution process may be checked by using alpha>1.

clc
clear
close all
clf

%specify the number of time steps over which the simulations may be performed
tend=600; 		%value of nmax
alpha=1; 		%specify alpha
mu0=4*pi*1e-7;   eps0=8.854e-12;
c=1/sqrt(mu0*eps0);
dx=1.5e-3;	 	%specify space resolution dx in meter
dt = alpha*dx/c;  	% compute time step dt
 
% compute the number of nodes. The upper limit is determined by tendXalpha
spc_end=300;
spc_indx=0:spc_end+1; 

e=zeros(length(spc_indx),1); h=e; ep=e; hp=h;     	%initialize previous values of e and h

hndl=plot(e);   			% get the handle to the plot
set(hndl,'EraseMode','xor') 
tic
for n = 1:tend      			%simulate over the time interval
    for i=1:length(spc_indx)-1   	%simulate over the nodes
        h(i)=hp(i)-dt*(e(i+1)-e(i))/(mu0*dx);	 % curl e = -del b/ del t;
    end
                    
    h(i+1)=hp(i)+((alpha-1)/(alpha+1))*(h(i)-hp(i+1)); % ABC for right boundary
    
    for i=2:length(spc_indx)-1  	 %simulate over the nodes
        e(i)=ep(i)-dt*(h(i)-h(i-1))/(eps0*dx);  % curl b = mu*epsilon*del e/ del t;
    end
   
	e(i+1)=ep(i)+((alpha-1)/(alpha+1))*(e(i)-ep(i+1)); 	% ABC for right boundary
    e(1)=ep(2)+((alpha-1)/(alpha+1))*(e(2)-ep(1)); 	% ABC for left boundary
        
    %%%%%%%%%%% Excitation of grid
    pwidth=500e-12; 		%specify pulsewidth in sec
    if n*dt <= pwidth 		%select this and next statement for rectangular pulse
        e(1) = 1;  			
%         if n*dt < pwidth*1.1     %select this and the next statement for Gaussian pulse
%             e(1)= exp(-((n*dt-pwidth/2)/(pwidth/4))^2);
%         end
    else
       e(1)=0;
    end
%     e1(1)=e(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hp=h;ep=e;  	%store the previous values in appropriate vectors         
    drawnow;           	%refresh data and plot the current condition of the pulse
    set(hndl,'YData',(e));

%specify the time instants at which the propagating pulse is to be displayed
     if n==50 || n==150 || n==250 || n==350 || n== 480 || n==600 
         figure;
         plot(e);
        xlabel('node number');
    	ylabel('amplitude');
     end
end
toc
