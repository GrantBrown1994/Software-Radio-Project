clear all;
close all;
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF1.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF2.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF3.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF4.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF5.mat');
load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF6.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF7.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF8.mat');

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
%% Find Timing Phase %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=450;
p_t = zeros(4*L, 1);
for tau=[0:4*L]
    p_t(tau+1)=mean(sum(abs(xBB(500+tau:L:500+L*n+tau)).^2));
end
tau=[0:4*L];
figure('Name', 'Ensamble Power of xBB')
plot(tau/Tb, p_t)
title('Ensamble Power of xBB')
fontsize(16,"points")
p_t_timing_phase = p_t(1:L/1.5);
[M, I] = maxk(abs(p_t_timing_phase),4);


xBBd=xBB(I:L/2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine Spectral Content of y %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'CTFT of xBB')
spec_analysis(xBB,1/Ts)
title('CTFT of xBB')
fontsize(16,"points")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Carrier Aquisition       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=32;
N1 = 25;
N2 = N1 + 2*N;
J_coarse = xBBd(N1:N1+(2*N) -1)'*xBBd(N1+(2*N):N1+(4*N)-1);
deltaFC_coarse = (1/(2*pi*N*Tb))*angle(J_coarse);

% %%%%%%%%%%%%%%%%%%%%%%
% %    9.47 Code       %
% %%%%%%%%%%%%%%%%%%%%%%
% N=32;
% N1 = 50;
% N2 = N1 + 2*N;
% J_coarse = 0;
% for n = N1:N2
%     J_coarse = J_coarse + xBBd(n+2*N)*conj(xBBd(n));
% end
% 
% deltaFC_coarse = (1/(2*pi*N*Tb))*angle(J_coarse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Remove Carrier Offset    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=(0:length(xBBd)-1)'*Tb;         % Set the time indices
xBBd=exp(-i*(2*pi*deltaFC_coarse*t)).*xBBd;


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
ryy = zeros(length(xBBd)-4*N,1);
for k=1:length(xBBd)-4*N
    ryy(k)=xBBd(k:k+2*N-1)'*xBBd(k+2*N:k+4*N-1);
end
figure('Name', 'Auto-Correlation for Pilot Detection')
plot(abs(ryy));
title('Auto-Correlation for Pilot Detection')
fontsize(16,"points")

epsilon=10;
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
y_pilot=preamble(length(preamble)-4*N: length(preamble)-2*N - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust Equalizer Weights       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=zeros(2*N,1);
e=zeros(2*N,1);
mu=0.9;
N=32;
for k=1:100000
        e(k)= cp(mod(k-1,N)+1) - (w'*y_pilot);
        w = w + 2*mu*conj(e(k))*y_pilot/(y_pilot'*y_pilot);
        y_pilot=circshift(y_pilot,-2);
end
figure('Name', 'Error of Equalization')
plot(abs(e));
title('Error of Equalization')
fontsize(16,"points")

figure('Name', 'Equalizer Weights')
plot(abs(w));
title('Equalizer Weights')
fontsize(16,"points")

[m,i] = max(w);
shift_length = length(w)/2 -i;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cross Correlation with Pilot    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(shift_length, 2) == 0
    xBBe = xBBe(1:2:end);
else 
    xBBe = xBBe(2:2:end);
end
ryy = zeros(length(xBBe)- N,1);
for k=1:length(xBBe)- N
    ryy(k) = xBBe(k:k+N-1)'*cp(1:32);
end
figure('Name', 'Cross-Correlation Graph for Pilot Detection')
plot(abs(ryy));
title('Cross-Correlation Graph for Pilot Detection')
fontsize(16,"points")
[M, I] = maxk(abs(ryy), 4);
I = max(I);

payload=xBBd(I+32:end);
xBBe_payload = conv(xBBd,conj(flip(w)));
if mod(shift_length, 2) == 0
    xBBe_payload = xBBe_payload(1:2:end);
else 
    xBBe_payload = xBBe_payload(2:2:end);
end
xBBe_payload = xBBe_payload(I+N:end);

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

phi=zeros(size(xBBe_payload)); s1=zeros(size(xBBe_payload)); mu=0.1;
for n=1:length(xBBe_payload)-1
    s1(n)=xBBe_payload(n)*exp(-1i*phi(n));
    s2=sign(real(s1(n)))+1i*sign(imag(s1(n)));
    s12=s1(n)*s2'; e=imag(s12)/real(s12);
    phi(n+1)=phi(n)+mu*e;
end

figure('Name', 'Constellation After Decision Directed')
plot(s1);
hold on;
plot(s1, 'xr');
hold off;
title('Constellation After Decision Directed')
fontsize(16,"points")

info_bits = QPSK2bits(s1);
data = bin2file(info_bits , 'Part4_Output_Fractional.txt');