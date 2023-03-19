% Function COHEGRAM: Signal coheregram (time dependancy of short-time MSC coherence)
%
% Usage: [C,time,freq]=cohegram(x,y[,flen,fstep,wlen, wstep,NFFT,fs]);
%
% Mandatory input parameters:
%   x      - 1st input signal
%   y      - 2nd input signal
%  
% Optional input parameters:                                       <Default values>
%   flen   - length of the frame for short-time analysis                     <1024>
%   fstep  - step of the frame for short-time analysis                     <flen/2>
%   wlen   - length of inner frame for MSC coherence estimation            <flen/4>
%   wstep  - step of inner frame for MSC coherence estimation              <wlen/2>
%   NFFT   - number of points for FFT computation        <minimum from (256, wlen)>
%   fs     - sampling frequency                                                 <2>
%
% Return value:
%   C      - matrix with time-dependent MSC
%   time   - vector of times related to the short-time frame beginnings 
%   freq   - vector of frequencies related to given FFT resolution 

%                              Made by PP
%                            CVUT FEL K13131
%                         Last change 8 Mar 2017

function [C,vtime,freq]=cohegram(x,y,flen,fstep,wlen,wstep,NFFT,fs);

if nargin<3,
  flen=1024;
end;
if nargin<4,
  fstep=flen/2;
end;
if nargin<5,
  wlen=flen/4;
end;
if nargin<6,
  wstep=wlen/2;
end;
if nargin<7,
  NFFT=max([256 wlen]);
end;
if nargin<8,
  fs=2;
end;


if ( length(x) ~= length(y) ),
  error('COHEGRAM: Vectors x and y must be same length!');
end
 
slen=  length(x);
wnum = floor((slen-flen)/fstep + 1) ;
  
C=[];

for i=1:wnum,

  framex=x((i-1)*fstep+1:(i-1)*fstep+flen);
  framey=y((i-1)*fstep+1:(i-1)*fstep+flen);

  [Cxy, freq] = mscohere ( framex,framey,wlen,wlen-wstep,NFFT,fs ) ;
  C=[C Cxy];

end;

vtime = (0:wnum-1)*fstep./fs ;

if nargout<1,
  pcolor(vtime,freq,C);
  xlabel('Time [s]')
  ylabel('Frequency [Hz]')
  shading flat ;
end
