clear; close all; clc;


sys1A = [1, -1.5, 0.7];
sys1B = [0.5, 0.3];

sys2A = [1, -0.8, 0.3];
sys2B = [0.4, 0.2];

N = 64e3;
fs = 8e3;

frequencies = (0:N-1)*fs/N;

[h1, ~] = freqz(sys1B, sys1A, N);
[h2, ~] = freqz(sys2B, sys2A, N);

figure;
subplot(2, 1, 1);
plot(frequencies, 20*log10(abs(h1)), 'b');
title('Frequency Response of System 1');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
subplot(2, 1, 2);
plot(frequencies, 20*log10(abs(h2)), 'r');
title('Frequency Response of System 2');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

