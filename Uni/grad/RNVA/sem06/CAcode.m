% Function CAcode generates C/A code for a given number PRN.
% It generates one period of the code.
% Code is provided in three formats
% 1. sequences containing 1, 0 - a vector length 1023 for each satellite
% 2. sequences containing 1, -1 (mapping from ad 1: 1->1, 0->-1)
%                             - a vector length 1023 for each satellite
%
% Syntax:
% PRNcode = CAcode(PRNnumber)
%
% PRNcode - structure containing vectors PRN10, PRN11, PRNhex
% PRN10(1, 1023) - Gold code in representation ad 1
% PRN11(1, 1023) - Gold code in representation ad 2
%


function PRN = CAcode(PRNno)
% Jiri Fajt, October 2001, Josef Spacek 2008

if (PRNno<1) || (PRNno>37)
    msgbox('Wrong PRN number requested.', 'Error', 'error');
   exit;
end

G1SR = ones(1,10);       % G1 shift register initial
G2SR = ones(1,10);       % G2 shift refister initial
G1mask = [0 0 1 0 0 0 0 0 0 1]';  % Shift register masks
G2mask = [0 1 1 0 0 1 0 1 1 1]';

G2jmask = zeros(10, 37);
G2jmask([2, 6], 1) = 1;          % PRN = 1
G2jmask([3, 7], 2) = 1;          % PRN = 2
G2jmask([4, 8], 3) = 1;          % PRN = 3
G2jmask([5, 9], 4) = 1;          % PRN = 4
G2jmask([1, 9], 5) = 1;          % PRN = 5
G2jmask([2, 10], 6) = 1;          % PRN = 6
G2jmask([1, 8], 7) = 1;          % PRN = 7
G2jmask([2, 9], 8) = 1;          % PRN = 8
G2jmask([3, 10], 9) = 1;          % PRN = 9
G2jmask([2, 3], 10) = 1;          % PRN = 10
G2jmask([3, 4], 11) = 1;          % PRN = 11
G2jmask([5, 6], 12) = 1;          % PRN = 12
G2jmask([6, 7], 13) = 1;          % PRN = 13
G2jmask([7, 8], 14) = 1;          % PRN = 14
G2jmask([8, 9], 15) = 1;          % PRN = 15
G2jmask([9, 10], 16) = 1;          % PRN = 16
G2jmask([1, 4], 17) = 1;          % PRN = 17
G2jmask([2, 5], 18) = 1;          % PRN = 18
G2jmask([3, 6], 19) = 1;          % PRN = 19
G2jmask([4, 7], 20) = 1;          % PRN = 20
G2jmask([5, 8], 21) = 1;          % PRN = 21
G2jmask([6, 9], 22) = 1;          % PRN = 22
G2jmask([1, 3], 23) = 1;          % PRN = 23
G2jmask([4, 6], 24) = 1;          % PRN = 24
G2jmask([5, 7], 25) = 1;          % PRN = 25
G2jmask([6, 8], 26) = 1;          % PRN = 26
G2jmask([7, 9], 27) = 1;          % PRN = 27
G2jmask([8, 10], 28) = 1;          % PRN = 28
G2jmask([1, 6], 29) = 1;          % PRN = 29
G2jmask([2, 7], 30) = 1;          % PRN = 30
G2jmask([3, 8], 31) = 1;          % PRN = 31
G2jmask([4, 9], 32) = 1;          % PRN = 32
G2jmask([5, 10], 33) = 1;          % PRN = 33
G2jmask([4, 10], 34) = 1;          % PRN = 34
G2jmask([1, 7], 35) = 1;          % PRN = 35
G2jmask([2, 8], 36) = 1;          % PRN = 36
G2jmask([4, 10], 37) = 1;          % PRN = 37

% generation G1 and G2 codes
for i = 1:1023
    G1(i) = G1SR(10);
    G2(i) = G2SR(10);
    G2j = mod(G2SR * G2jmask(:,PRNno), 2);
    
    G1SR = [mod(G1SR * G1mask, 2), G1SR(1:9)];
    G2SR = [mod(G2SR * G2mask, 2), G2SR(1:9)];
    
    PRN10(i) = xor(G1(i), G2j); 
end

% conversion to 1, -1 sequence

PRN11 = 2 * PRN10 - 1;

PRN.PRN10 = PRN10;
PRN.PRN11 = PRN11;

%%%%%% end of CAcode.m %%%%%%