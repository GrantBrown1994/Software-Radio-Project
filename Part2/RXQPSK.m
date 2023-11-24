close all;
clear all;
% load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF1.mat');
load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF2.mat');
% load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF3.mat');
% load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF4.mat');
% load('/Users/grantbrown/Library/Mobile Documents/com~apple~CloudDocs/Documents_UofU/Software Radio/CD/xRF5.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine Spectral Content of xRF %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_analysis(xRF,1/Ts)

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
spec_analysis(xBB,1/Ts)

%%%%%%%%%%%%%%%%%%%%%%
% DECIMATION         %
%%%%%%%%%%%%%%%%%%%%%%

figure(2)
hold on
plot(xBBd,'b')
plot(xBBd(1:L:end),'r.')
axis('square')
xlabel('real part')
ylabel('imaginary part')
hold off

%%%%%%%%%%%%%%%%%%%%%%
% Extract Payload   %
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Detection of s[n] (pilot)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 32;
starting=N+N+2; ending=length(xBBd)-starting;
i=1;
for n=starting:ending
    ryy_curr = 0;
    for k=0:N
        ryy_curr = ryy_curr + xBBd(n-k-1)*conj(xBBd(n-N-1-k));
    end
    ryy(i) = ryy_curr;
    i = i+1;
end
plot(abs(ryy));

epsilon=0.05;
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

payload=xBBd(end_point:end);
figure(3)
hold on
plot(payload,'b')
plot(payload,'r.')
axis('square')
xlabel('real part')
ylabel('imaginary part')
hold off


preamble=xBBd(1:end_point+1);

pilot=preamble(length(preamble)-2*N: length(preamble)-(1*N)-1);
xBBd_pilot=xBBd(length(preamble)-2*N: length(preamble)-(1*N)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust Equalizer Weights       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=zeros(N,1);
e=zeros(N,1);
mu=0.005;
for k=1:100000
        e(k)= pilot(mod(k-1,N)+1) - (w'*xBBd_pilot);
        w = w + 2*mu*conj(e(k))*xBBd_pilot;
        xBBd_pilot=circshift(xBBd_pilot,-1);
end
plot(abs(e));
plot(abs(w));

[m,i] = max(w);
w=circshift(w,(length(w)/2) -i);
plot(abs(w));

xBBe = conv(xBBd,conj(flip(w)));

figure(5);
subplot(2,1,1);
plot(xBBd);
hold on;
plot(xBBd, 'xr');
hold off;
title('Eye diagram before Equalization')
subplot(2,1,2)
plot(xBBe);
hold on;
plot(xBBe, 'xr');
hold off;
title('Eye diagram after Equalization')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Incorrect Timing Phase         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xBBd=xBB(4:L:end);
xBBe = conv(xBBd,conj(flip(w)));
figure(6);
subplot(2,1,1);
plot(xBBd);
hold on;
plot(xBBd, 'xr');
hold off;
title('Eye diagram before Equalization')
subplot(2,1,2)
plot(xBBe);
hold on;
plot(xBBe, 'xr');
hold off;
title('Eye diagram after Equalization')
hold off

