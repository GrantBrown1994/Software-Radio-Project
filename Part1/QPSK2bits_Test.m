%Test
clear all;
close all;
N=1000;
bits = randi([0 1],N,1);
sym = bits2QPSK(bits);
data = QPSK2bits(sym);
if bits == data
    disp('Arrays are Equal')
else
    disp('Code is wrong')
end
