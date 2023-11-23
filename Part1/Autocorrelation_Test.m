%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Construction of s[n]       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=31;           % Equalizer order
N=N+1;          % Equalizer length
pilot=CycPilot(N-1);
s=pilot;
for k=1:3
    s=[s; s];
end
s=[s;sign(randn(1000,1)) + i*sign(randn(1000,1))];

N = 32;
starting=N+N+2; ending=length(s)-starting;
i=1;
for n=starting:ending
    ryy_curr = 0;
    for k=0:N
        ryy_curr = ryy_curr + s(n-k-1)*conj(s(n-N-1-k));
    end
    ryy(i) = ryy_curr;
    i = i+1;
end
plot(abs(ryy))