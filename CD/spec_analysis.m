%
% Spectrum analysis: [X,F]=s_analysis(x,fs)
% This function calculates the magnitude spectrum of x, sampled at a rate fs.
% 
function [XX,F]=spec_analysis(x,fs)
N = length(x);                      % the length of input
Ts=1/fs;                            % sampling interval
[X,F] = freqz(Ts*x,1,N,'whole',fs); % perform Fourier transform
F=F-fs/2;                           % adjust the frequency range to (-fs/2,+fs/2)
X=fftshift(abs(X));                 % take amplitude and shift the samples of
                                    %  the spectrum to the correct position
if nargout==0                                    
    t=[0:N-1]*Ts;                       % time scale
    plot(F,X)           % plot the spectrum
    xlabel('frequency (Hz)'); 
    ylabel('magnitude')                 % label the axes
else
    XX=X;
end
