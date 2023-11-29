clear all;
close all;
load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF1.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine Spectral Content of xRF %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_analysis(xRF,1/Ts)

phic=0;                 % carrier phase offset
Dfc=0;                 % carrier frequency offset (unknown to the receiver)
L=100;

%%%%%%%%%%%%%%%%%%%%%%
%    DEMODULATION    %
%%%%%%%%%%%%%%%%%%%%%%
t=[0:length(xRF)-1]'*Ts;         % Set the time indices
xbbRF=2*exp(-i*(2*pi*(fc+Dfc)*t-phic)).*xRF;

%%%%%%%%%%%%%%%%%%%%%%
% RECEIVE FILTERING  %
%%%%%%%%%%%%%%%%%%%%%%
pR=pT;    
xBB=conv(xbbRF,conj(pT));
xBBd=xBB(1:L:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine Spectral Content of y %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_analysis(xBB,1/Ts)

%%%%%%%%%%%%%%%%%%%%%%
% DECIMATION         %
%%%%%%%%%%%%%%%%%%%%%%

figure(2)
hold on
plot(xBBd,'b')
plot(xBBd,'r.')
axis('square')
xlabel('real part')
ylabel('imaginary part')
hold off

%%%%%%%%%%%%%%%%%%%%%%
% Extract Payload   %
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Detection of s[n] (pilot)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 32;
starting=N+N+2; ending=length(xBBd)-starting;
i=1;
for n=starting:ending
    ryy_curr = 0;
    for k=0:N
        ryy_curr = ryy_curr + xBBd(n-k-1)*conj(xBBd(n-N-1-k));
    end
    ryy(i) = ryy_curr;
    i = i+1;
end
plot(abs(ryy));

figure(6);
for i=1:length(xBBd)- N
    xBBd_cp = xBBd(i:i+N-1);
    rxx_curr = 0;
    for k=1:N
        rxx_curr = rxx_curr + xBBd_cp(k)*conj(cp(k));
    end
    rxx(i) = rxx_curr;
end
plot(abs(rxx));

payload=xBBd(106+31:end);
figure(3)
hold on
plot(payload,'b')
plot(payload,'xr')
axis('square')
xlabel('real part')
ylabel('imaginary part')
hold off

info_bits = QPSK2bits(payload);
data = bin2file(info_bits , 'Part1_Output.txt');

