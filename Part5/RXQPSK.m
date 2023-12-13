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
j=1
for tau=[10*L:15*L]
    %p_t(tau+1)=mean(abs(xBB(L+tau:L:length(xBB)-L)).^2);
    p_t(j)=mean(abs(xBB(500+tau:L:500+tau+10*L)).^2);
    j=j+1;
end
figure('Name', 'Ensamble Power of xBB')
plot(p_t)
title('Ensamble Power of xBB')
fontsize(16,"points")
[M, packet_start] = max(p_t);
packet_start = packet_start - roundn(packet_start,2);
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
%I = (max(I)+N)*L;
%I=(max(I)+N)*L -payload_start;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  TIMING RECOVERY: Decision Directed   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=0.01;dtau=6;
Ly=length(xBB);
kk=1;
yp=0;ym=0;
start=1;
tau=0.001*ones(1,floor((Ly-start)/L));
% for k=I+payload_start:L:length(tau)*L-1
%     tauTb=round(tau(kk)*L); %Moves sample offset to sample value in discrete time
%     sk(kk)=slicer(xBB(k+tauTb),4); %moves it to closest constellation point
%     tau(kk+1)=tau(kk)+mu*real(conj(sk(kk)-xBB(k+tauTb))*(xBB(k+tauTb+dtau)-xBB(k+tauTb-dtau)));
%     kk=kk+1;
% end
test_xBBd = xBBd(137:end);
for k=I:L:length(tau)*L+L
    tauTb=round(tau(kk)*L);
    tauTB_array(kk)=tauTb;
    if k+tauTb+dtau > length(xBB)
        break
    end
    test_xBB(kk)=xBB(k+tauTb);
    sk(kk)=slicer(xBB(k+tauTb),4);
    error = real(sk(kk)-xBB(k+tauTb));
    tau(kk+1)=tau(kk)+mu*real((sk(kk)-xBB(k+tauTb))*(xBB(k+tauTb+dtau)-xBB(k+tauTb-dtau))');
    kk=kk+1;
end
figure, axes('position',[0.1 0.25 0.8 0.5]), plot(tau(1:kk-1),'k')
xlabel('Iteration Number, n'), ylabel('\tau[n]')
plot(xBBd(137:end));

for k=1:length(test_xBBd)
    test_slicer(k) = slicer(test_xBBd(k), 4);
end
test_slicer=test_slicer.';
xBBd=sk;
sk=sk.';

index=1;
for k=1:length(sk)
    if sk(k) == test_slicer(k)
        continue
    else
        mismitch_index(index)=k;
        index = index + 1;
    end 
end

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

info_bits = QPSK2bits(xBBd);
data = bin2file(info_bits , 'Part5_Output.txt');

info_bits = QPSK2bits(test_xBBd);
data = bin2file(info_bits , 'Part5_Output_Test.txt');

