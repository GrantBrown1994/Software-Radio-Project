%
% QPSK transmitter
%
clear all
close all
%
% parameters
%
Tb=1e-4 ;               % (nominal) symbol rate
fb=1/Tb;    
epsilon=0;           % fractional error in symbol rate
L=100;                  % oversampling factor
Ts=Tb/L;                % sampling period
fs=1/Ts;
fc=1e5;                 % (nominal) carrier frequency
phic=0;                 % carrier phase offset
Dfc=0;                 % carrier frequency offset (unknown to the receiver)
sigmav2=0.0001;           % variance of channel noise (unknown to the receiver)
sigmav=sqrt(sigmav2);
alpha=0.5;              % transmit filter roll-off factor 
K=8;                    % length of the transmit filter in unit of symbol intervals
c=1;                    % multipath channel
c=[1 zeros(1,91) 0.4];
c=[1 zeros(1,60) 2 zeros(1,123) 0.25]; 
c=[1 zeros(1,67) 0.75 zeros(1,145) 0.4]; 
%c=[1 zeros(1,75) 0.6 zeros(1,103) 0.2]; 
s=sqrt(2)*CycPilot(31); % get one period of the cyclic preamble; length = 31+1 = 32 
cp=sqrt(2)*CycPilot(31);        % get one period of the cyclic preamble; length = 31+1 = 32 
preamble=[cp; cp; cp; cp];  % constrauct the preamble
%
% source
%
info_bits=file2bin('testtext.txt')';      % convert source to bits
%
% generate transmit signal
%     
s=bits2QPSK(info_bits);          % convert bits to symbols
s=[preamble; s];                 % add preamble
pT=sr_cos_p(L*K,L,alpha);        % transmitter filter
xBB=conv(pT,expander(s,L));      % baseband pulse-shaping 
xBB=resample(xBB,epsilon);       % change the sampling rate fs to (1+epsilon)fs
t=[0:length(xBB)-1]'*Ts;         % time
xRF=real(xBB.*exp(j*(2*pi*(fc+Dfc)*t+phic)));       % modulation
xRF=conv(c,xRF);                            % pass the RF signal through a multipath channel
%
% Add adjacent channels
%
user1=1.5*bits2QPSK(sign(randn(2*length(s),1)));
user1=conv(pT,expander(user1,L)); 
t=[0:length(user1)-1]'*Ts;
user1=real(user1.*exp(-j*2*pi*(fc+3*fb)*t));
user2=0.5*bits2QPSK(sign(randn(2*length(s),1)));
user2=conv(pT,expander(user2,L));
user2=real(user2.*exp(-j*2*pi*(fc-3*fb)*t));
xRF=xRF+[user1+user2; zeros(length(xRF)-length(user1),1)]; 
%
% Add noise  
%
xRF=xRF+sigmav*randn(size(xRF));
spec_analysis(xRF,fs)

N = 32;
starting=N+N+2; ending=length(s)-starting;
i=1;
for n=starting:ending
    ryy_curr = 0;
    for k=0:N
        ryy_curr = ryy_curr + s(n-k-1)*conj(s(n-N-1-k));
    end
    ryy(i) = ryy_curr;
    i = i+1;
end
plot(abs(ryy))