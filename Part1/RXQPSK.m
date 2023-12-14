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
for tau=[11*L:17*L]
     p_t(tau)=mean(abs(xBB(500+tau:L:500+tau+200*L)).^2);
end
figure('Name', 'Ensamble Power of xBB')
plot(p_t)
title('Ensamble Power of xBB')
fontsize(16,"points")
[M, packet_start] = max(p_t);
packet_start = packet_start - roundn(packet_start,2);
if packet_start <= 0
    packet_start = packet_start + 100;
end
xBBd=xBB(packet_start:L:end);

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
N=32;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Incorrect Timing Phase         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xBBd=xBB(40:L:end);
figure('Name', 'Incorrect Timing Phase Constellation')
plot(xBBd);
hold on;
plot(xBBd, 'xr');
hold off;
title('Decimated Constellation before Equalization')
fontsize(16,"points")

