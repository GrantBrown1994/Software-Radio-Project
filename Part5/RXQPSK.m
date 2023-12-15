clear all;
close all;
%load('../CD/xRF1.mat');
load('../CD/xRF9.mat');

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
xBB=conv(xbbRF,pT);

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
    p_t(tau)=mean(abs(xBB(500+tau:L:500+tau+100*L)).^2);
end
figure('Name', 'Ensamble Power of xBB')
plot(p_t)
title('Ensamble Power of xBB')
fontsize(16,"points")
[M, packet_start] = max(p_t);
packet_start = packet_start - floor(packet_start/100)*100;
if packet_start < 0
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


N=32;
ryy = zeros(length(xBBd)- N,1);
for k=1:length(xBBd)- N
    ryy(k) = xBBd(k:k+N-1)'*cp(1:32);
end
figure('Name', 'Cross-Correlation Graph for Pilot Detection')
plot(abs(ryy));
title('Cross-Correlation Graph for Pilot Detection')
fontsize(16,"points")
[M, I] = maxk(abs(ryy), 4);
I = (max(I)+N-1)*L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  TIMING RECOVERY: Decision Directed   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=0.01;dtau=6;
Ly=length(xBB);
kk=1;
yp=0;ym=0;
start=1;
tau=0.1*ones(1,floor((Ly-start)/L));
for k=I:L:length(tau)*L
    tauTb=round(tau(kk)*L);
    tauTB_array(kk)=tauTb;
    if k+tauTb+dtau > length(xBB)
        break
    end
    sk(kk)=slicer(xBB(k+tauTb),4);
    error = real(sk(kk)-xBB(k+tauTb));
    tau(kk+1)=tau(kk)+mu*real((sk(kk)-xBB(k+tauTb))*(xBB(k+tauTb+dtau)-xBB(k+tauTb-dtau))');
    kk=kk+1;
end
figure, axes('position',[0.1 0.25 0.8 0.5]), plot(tau(1:kk-1),'k')
xlabel('Iteration Number, n'), ylabel('\tau[n]')
xBBd = sk;

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

info_bits = QPSK2bits(xBBd);
data = bin2file(info_bits , 'Part5_Output.txt');

