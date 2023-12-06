clear all;
close all;
load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF1.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine Spectral Content of xRF %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'CTFT of xRF')
spec_analysis(xRF,1/Ts)
title('CTFT of xRF')
fontsize(16,"points")

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
figure('Name', 'CTFT of xBB')
spec_analysis(xBB,1/Ts)
title('CTFT of xBB')
fontsize(16,"points")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Timing Phase %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=450;
p_t = zeros(4*L, 1);
j=1;
for tau=[0:4*L]
    p_t(j)=mean(sum(abs(xBB(500+tau:L:500+L*n+tau)).^2));
    j=j+1;
end
tau=[0:4*L];
figure('Name', 'Ensamble Power of xBB')
plot(tau/Tb, p_t)
title('Ensamble Power of xBB')
fontsize(16,"points")


%%%%%%%%%%%%%%%%%%%%%%
% DECIMATION         %
%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Constellation of xBBd')
hold on
plot(xBBd,'b')
plot(xBBd,'r.')
axis('square')
xlabel('real part')
ylabel('imaginary part')
hold off
title('Constellation of xBBd')
fontsize(16,"points")

%%%%%%%%%%%%%%%%%%%%%%
% Extract Payload   %
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Detection of s[n] (pilot)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 32;
for k=1:length(xBBd)-2*N
    ryy(k)=xBBd(k:k+N-1)'*xBBd(k+N:k+2*N-1);
end
figure('Name', 'Auto-Correlation Graph for Pilot Detection')
plot(abs(ryy));
title('Auto-Correlation Graph for Pilot Detection')
fontsize(16,"points")


for i=1:length(xBBd)- N
    xBBd_cp = xBBd(i:i+N-1);
    rxx_curr = 0;
    for k=1:N
        rxx_curr = rxx_curr + xBBd_cp(k)*conj(cp(k));
    end
    rxx(i) = rxx_curr;
end
figure('Name', 'Cross-Correlation Graph for Pilot Detection')
plot(abs(rxx));
title('Cross-Correlation Graph for Pilot Detection')
fontsize(16,"points")
[M, I] = maxk(abs(rxx), 4);
I = max(I);

payload=xBBd(I+32:end);
figure('Name', 'Payload Constellation')
hold on
plot(payload,'b')
plot(payload,'xr')
axis('square')
xlabel('real part')
ylabel('imaginary part')
hold off
title('Payload Constellation')
fontsize(16,"points")

info_bits = QPSK2bits(payload);
data = bin2file(info_bits , 'Part1_Output.txt');

