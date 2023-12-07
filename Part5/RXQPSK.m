clear all;
close all;
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF1.mat');
load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF9.mat');

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  TIMING RECOVERY: Decision Directed   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=0.05;dtau=6;
Ly=length(xBB);
kk=1;
yp=0;ym=0;
start=1;
tau=0.3*ones(1,floor((Ly-start)/L));
x=tau;
for k=start:L:length(tau)*L-1
    tauTb=round(tau(kk)*L); %Moves sample offset to sample value in discrete time
    sk=slicer(xBB(k+tauTb),4); %moves it to closest constellation point
    tau(kk+1)=tau(kk)+mu*real(conj(sk-xBB(k+tauTb))*(xBB(k+tauTb+dtau)-xBB(k+tauTb-dtau)));
    kk=kk+1;
end
tau_final=tau(end);
I=round(abs(tau_final*L));
figure, axes('position',[0.1 0.25 0.8 0.5]), plot(tau(1:kk-1),'k')
xlabel('Iteration Number, n'), ylabel('\tau[n]')

xBBd=xBB(I:L:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine Spectral Content of y %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'CTFT of xBB')
spec_analysis(xBB,1/Ts)
title('CTFT of xBB')
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

