close all;
clear all;
load('../CD/xRF1.mat');
%load('../CD/xRF2.mat');
%load('../CD/xRF3.mat');
%load('../CD/xRF4.mat');
%load('../CD/xRF5.mat');
%load('../CD/xRF6.mat');
%load('../CD/xRF7.mat');
%load('../CD/xRF8.mat');

%load('../CD/xRF2ans.mat');
%load('../CD/xRF3ans.mat');
%load('../CD/xRF4ans.mat');
%load('../CD/xRF5ans.mat');
%load('../CD/xRF6ans.mat');

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
for tau=[11*L:17*L]
     p_t(tau)=mean(abs(xBB(500+tau:L:500+tau+100*L)).^2);
end
figure('Name', 'Ensamble Power of xBB')
plot(p_t)
title('Ensamble Power of xBB')
fontsize(16,"points")
[M, packet_start] = max(p_t);
packet_start = packet_start - floor(packet_start/100)*100;
if packet_start <= 0
    packet_start = packet_start + 100;
end
timing_phase=0;
xBBd=xBB(packet_start+timing_phase:L:end);

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
N1 = 50;
N2 = N1 + 4*N;
J_coarse = 0;
for n = N1:N2
    J_coarse = J_coarse + xBBd(n+N)*conj(xBBd(n));
end

deltaFC_coarse = (1/(2*pi*N*Tb))*angle(J_coarse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Remove Carrier Offset    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=[0:length(xBBd)-1]'*Tb;         % Set the time indices
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
ryy = zeros(length(xBBd),1);
% for k=1:length(xBBd)-4*N
%     ryy(k)=xBBd(k:k+2*N-1)'*xBBd(k+2*N:k+4*N-1);
% end
for n=2*N+1:length(xBBd)
    ryy(n) = xBBd(n-1:-1:n-N)'*xBBd(n-N-1:-1:n-2*N);
end
figure('Name', 'Auto-Correlation for Pilot Detection')
plot(abs(ryy));
title('Auto-Correlation for Pilot Detection')
fontsize(16,"points")

ryy_start=2*N+1;
ryy_deriv = abs(conv(ryy, [1 -1]));
plot(ryy_deriv)
epsilon = 4.25;
flat_top_length=0;
for k=ryy_start:length(ryy_deriv)
    if ryy_deriv(k) < epsilon && flat_top_length == 0
        flat_top_start=k;
        flat_top_length = 1;
    elseif  ryy_deriv(k) < epsilon && flat_top_length > 0
        flat_top_length = flat_top_length + 1;
    elseif ryy_deriv(k) > epsilon && flat_top_length <= N
        flat_top_length = 0;
    elseif ryy_deriv(k) > epsilon && flat_top_length >= N
        flat_top_end=k;
        break
    end
end

y_pilot=xBBd(flat_top_start+20+N-1:-1:flat_top_start+20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust Equalizer Weights       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=zeros(N,1);
e=zeros(N,1);
mu=0.005;
for k=1:100000
        e(k)= cp(mod(k-1,N)+1) - (w'*y_pilot);
        w = w + 2*mu*conj(e(k))*y_pilot/(y_pilot'*y_pilot);
        y_pilot=circshift(y_pilot,-1);
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
shift_length = length(w)/2 -i;
w=circshift(w,(length(w)/2) -i);
figure('Name', 'Centered Equalizer Weights')
plot(abs(w));
title('Centered Equalizer Weights')
fontsize(16,"points")

xBBe = conv(xBBd,conj(w));

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
for i=1:length(xBBe)- N
    rxx(i)=xBBe(i:i+N-1)'*cp(1:N);
end
figure('Name', 'Cross-Correlation Graph for Pilot Detection')
plot(abs(rxx));
title('Cross-Correlation Graph for Pilot Detection')
fontsize(16,"points")
[M, I] = maxk(abs(rxx), 4);
I = max(I);

payload=xBBd(I+32:end);
xBBe_payload = conv(xBBd,conj(w));
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

phi=zeros(size(xBBe_payload)); s1=zeros(size(xBBe_payload)); mu=0.1;
for n=1:length(xBBe_payload)-1
    s1(n)=xBBe_payload(n)*exp(-1i*phi(n));
    s2=sign(real(s1(n)))+1i*sign(imag(s1(n)));
    s12=s1(n)*s2'; e=imag(s12)/real(s12);
    phi(n+1)=phi(n)+mu*e;
end
xBBe_phi=s1;

figure('Name', 'Constellation After Decision Directed')
plot(xBBe_phi);
hold on;
plot(xBBe_phi, 'xr');
hold off;
title('Constellation After Decision Directed')
fontsize(16,"points")

info_bits = QPSK2bits(xBBe_phi);
data = bin2file(info_bits , 'Part4_Output.txt');

