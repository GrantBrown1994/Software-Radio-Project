clear all;
close all;
load('../CD/xRF1.mat');
%load('../CD/xRF2.mat');
%load('../CD/xRF3.mat');
%load('../CD/xRF4.mat');
%load('../CD/xRF5.mat');

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
timing_phase=0;
xBBd=xBB(1+timing_phase:L/2:end);

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
N2 = N1 + 4*N;
J_coarse = xBBd(N1:N1+(2*N) -1)'*xBBd(N1+(2*N):N1+(4*N)-1);
deltaFC_coarse = ((1/(2*pi*2*N*Tb))*angle(J_coarse));

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
ryy = zeros(length(xBBd),1);
% for k=1:length(xBBd)-4*N
%     ryy(k)=xBBd(k:k+2*N-1)'*xBBd(k+2*N:k+4*N-1);
% end
for n=4*N+1:length(xBBd)
    %ryy(n) = 
    ryy_curr = 0;
    for k=0:2*N-1
        ryy_curr = ryy_curr + xBBd(n-k-1)*conj(xBBd(n-2*N-k));
    end
    ryy(n) = ryy_curr;
end
figure('Name', 'Auto-Correlation for Pilot Detection')
plot(abs(ryy));
title('Auto-Correlation for Pilot Detection')
fontsize(16,"points")

ryy_start=4*N+1;
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

if mod(flat_top_start, 2) == 1
    flat_top_start = flat_top_start - 1;
end
y_pilot=xBBd(flat_top_start+2*N-1:-1:flat_top_start);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust Equalizer Weights       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=zeros(2*N,1);
e=zeros(2*N,1);
mu=0.1;
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

xBBd_payload = xBBd(I+32:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust Equalizer Weights       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N1=31;          % Equalizer order
c=[1];
nn=1;
Psi_inv=10000*eye(N);
lambda=0.75;
Delta=round((N1/2+(length(c)+length(pT))/L)/2);
w=w(1:2:end);
for n=N:2:length(xBBd_payload)-N
    tdl=xBBd_payload(n:-1:n-N+1);
    u=Psi_inv*tdl;
    k=u/(lambda+tdl'*u);
    x(n)=w'*tdl;
    s(n) = slicer(x(n),4);
    e=s(n)-x(n);
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
shift_length = length(w)/2 -i;
w=circshift(w,(length(w)/2) -i);
figure('Name', 'Centered Equalizer Weights')
plot(abs(w));
title('Centered Equalizer Weights')
fontsize(16,"points")

xBBe_payload = conv(xBBd_payload,conj(w));
if mod(shift_length, 2) == 0
    xBBe_payload = xBBe_payload(1:2:end);
else 
    xBBe_payload = xBBe_payload(2:2:end);
end


figure('Name', 'Decimated vs Equalized Constellations')
subplot(2,1,1);
plot(xBBd_payload);
hold on;
plot(xBBd_payload, 'xr');
hold off;
title('Decimated Constellation before Equalization')
fontsize(16,"points")
subplot(2,1,2)
plot(xBBe_payload);
hold on;
plot(xBBe_payload, 'xr');
hold off;
title('Decimated Constellation after Equalization')
fontsize(16,"points")
hold off

info_bits = QPSK2bits(xBBe_payload);
data = bin2file(info_bits , 'Part6_Output_Fractional.txt');
