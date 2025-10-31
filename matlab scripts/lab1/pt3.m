clear; close all; clc;

% Low and high operating modes (transfer function coefficients)
% G_low(z)  = (0.5 z + 0.15) / (z^2 - 0.4 z + 0.05)
% G_high(z) = (0.29 z + 0.27) / (z^2 - 1.25 z + 0.8)

% Denominator and numerator coefficient vectors (MATLAB convention)
sys_lowA  = [1, -0.4, 0.05];
sys_lowB  = [0.5, 0.15];

sys_highA = [1, -1.25, 0.8];
sys_highB = [0.29, 0.27];

N = 64e3;
fs = 8e3;

frequencies = (0:N-1)*fs/N;

[h_low, ~]  = freqz(sys_lowB,  sys_lowA,  N);
[h_high, ~] = freqz(sys_highB, sys_highA, N);

figure;
subplot(2, 1, 1);
plot(frequencies, 20*log10(abs(h_low)), 'b');
title('Frequency Response of G\_low');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
subplot(2, 1, 2);
plot(frequencies, 20*log10(abs(h_high)), 'r');
title('Frequency Response of G\_high');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

