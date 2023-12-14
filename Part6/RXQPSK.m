close all;
clear all;
load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF3.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF8.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF2.mat');
%load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF10.mat');

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
epsilon = 5;
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

for k=1:length(xBBd)
    s(k) = slicer(xBBd(k),4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust Equalizer Weights       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N1=31;          % Equalizer order
c=[1];
Delta=round((N+(length(c)+2*length(pT))/L)/2);
w=zeros(N,1); nn=1;
Psi_inv=10000*eye(N);
lambda=0.95;
for n=N:length(xBBd)-Delta
    tdl=xBBd(n:-1:n-N1);
    u=Psi_inv*tdl;
    k=u/(lambda+tdl'*u);
    x(n)=w'*tdl;
    e=s(n-Delta)-x(n);
    w=w+k*e';
    Psi_inv=(1/lambda)*(Psi_inv-k*(tdl'*Psi_inv));
    xi(nn)=abs(e)^2;
    nn=nn+1;
    if rem(n,100)==0 
        figure(1),plot(x(n-99:n),'.'),pause(0.1) 
    end
end
figure(2),semilogy(xi)

% figure('Name', 'Error of Equalization')
% plot(abs(e));
% title('Error of Equalization')
% fontsize(16,"points")

figure('Name', 'Equalizer Weights')
plot(abs(w));
title('Equalizer Weights')
fontsize(16,"points")

[m,i] = max(w);
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

for i=1:length(xBBe)- N
    xBBe_cp = xBBe(i:i+N-1);
    rxx_curr = 0;
    for k=1:N
        rxx_curr = rxx_curr + xBBe_cp(k)*conj(cp(k));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  TIMING RECOVERY: Decision Directed   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=0.01;dtau=6;
Ly=length(xBBe_payload);
kk=1;
yp=0;ym=0;
start=1;
tau=0.1*ones(1,floor((Ly-start)));
for k=I:length(tau)
    tauTb=round(tau(kk));
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


info_bits = QPSK2bits(xBBe_payload);
data = bin2file(info_bits , 'Part6_Output.txt');

