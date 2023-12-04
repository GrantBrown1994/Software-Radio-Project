close all;
clear all;
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF1.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF2.mat');
load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF3.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF4.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF5.mat');

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
%%%%%%%%%%%%%%%%%%%%%%
% Extract Payload   %
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Detection of s[n] (pilot)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 32;
for k=1:length(xBBd)-4*N
    ryy(k)=xBBd(k:k+2*N-1)'*xBBd(k+2*N:k+4*N-1);
end
figure('Name', 'Auto-Correlation Graph for Pilot Detection')
plot(abs(ryy));
title('Auto-Correlation Graph for Pilot Detection')
fontsize(16,"points")

figure('Name', 'Auto-Correlation Graph for Pilot Detection')
plot(abs(ryy));
title('Auto-Correlation Graph for Pilot Detection')
fontsize(16,"points")

epsilon=5;
preamble_started=false;
end_point=0;
starting_point=0;
for n=1+N:length(ryy)
    if preamble_started == false && abs(abs(ryy(n))- abs(ryy(n-N))) <= epsilon
        preamble_started = true;
        starting_point = n-N;
    elseif preamble_started == true && abs(abs(ryy(n))- abs(ryy(n-N))) >= epsilon
        end_point = n-1;
        break
    end
end

preamble=xBBd(1:end_point);
pilot=preamble(length(preamble)-2*N: length(preamble)-N - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust Equalizer Weights       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=zeros(N,1);
e=zeros(N,1);
mu=0.00003;
for k=1:100000
        e(k)= cp(mod(k-1,N)+1) - (w'*pilot);
        w = w + 2*mu*conj(e(k))*pilot;
        pilot=circshift(pilot,-1);
end
figure('Name', 'Error of Equalization')
plot(abs(e));
title('Error of Equalization')
fontsize(16,"points")

figure('Name', 'Equalizer Weights')
plot(abs(w));
title('Equalizer Weights')
fontsize(16,"points")

% Y = zeros(32,32);
% for k=1:32
%     Y(:,k) = pilot;
%     pilot=circshift(pilot,-1);
% end
% Y_H = Y';
% I = epsilon * eye(32,32);
% Y_w = Y*Y_H + I;
% Y_w = inv(Y_w);
% Y_s = Y*conj(cp);
% w = Y_w*Y_s;

[m,i] = max(w);
w=circshift(w,(length(w)/2) -i);
figure('Name', 'Centered Equalizer Weights')
plot(abs(w));
title('Centered Equalizer Weights')
fontsize(16,"points")

xBBe = conv(xBBd,conj(flip(w)));

figure('Name', 'Decimated vs Equalized Constellations')
subplot(2,1,1);
plot(xBBd);
hold on;
plot(xBBd, 'xr');
hold off;
title('Decimated Constellation before Equalization')
fontsize(16,"points")
subplot(2,1,2)
plot(xBBe);
hold on;
plot(xBBe, 'xr');
hold off;
title('Decimated Constellation after Equalization')
fontsize(16,"points")
hold off

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
[M, I] = maxk(abs(ryy), 4);
I = max(I);

payload=xBBd(I+32:end);
xBBe_payload = conv(xBBd,conj(flip(w)));
xBBe_payload = xBBe_payload(I+32:end);
figure('Name', 'Payload vs Equalized Payload Constellation')
subplot(2,1,1);
plot(payload);
hold on;
plot(payload, 'xr');
hold off;
title('Payload Constellation before Equalization')
fontsize(16,"points")
subplot(2,1,2)
plot(xBBe_payload);
hold on;
plot(xBBe_payload, 'xr');
hold off;
title('Payload Constellation after Equalization')
fontsize(16,"points")
hold off

info_bits = QPSK2bits(xBBe_payload);
data = bin2file(info_bits , 'Part2_Output.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Incorrect Timing Phase         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xBBd=xBB(4:L:end);
xBBe = conv(xBBd,conj(flip(w)));
figure('Name', 'Incorrect Timing Phase Constellation')
subplot(2,1,1);
plot(xBBd);
hold on;
plot(xBBd, 'xr');
hold off;
title('Decimated Constellation before Equalization')
fontsize(16,"points")
subplot(2,1,2)
plot(xBBe);
hold on;
plot(xBBe, 'xr');
hold off;
title('Decimated Constellation after Equalization')
fontsize(16,"points")
hold off

