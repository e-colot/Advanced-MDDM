clear; close all; clc;

L0 = 5e-6;
R0 = 10;
C0 = 3e-9;

Ac = [0 1/C0; -1/L0 -R0/L0];
Bc = [0; 1/L0];
Cc = [1 0];
Dc = 0;

sysContinuous = ss(Ac, Bc, Cc, Dc);

poles = eig(Ac);
disp('Poles of the continuous-time system:');
disp(poles);

fs = 4e6;

sysDiscreteTustin = c2d(sysContinuous, 1/fs, 'tustin');
sysDiscreteZOH = c2d(sysContinuous, 1/fs, 'zoh');

Ad = sysDiscreteZOH.A;
Bd = sysDiscreteZOH.B;
Cd = sysDiscreteZOH.C;
Dd = sysDiscreteZOH.D;

% step response comparison
tc = 0:1/(20*fs):0.00001;
td = 0:1/fs:0.00001;

[continuousStep, tc] = step(sysContinuous, tc);
[discreteStepTustin, td] = step(sysDiscreteTustin, td);
[discreteStepZOH, td] = step(sysDiscreteZOH, td);

figure;
plot(tc*1e6, continuousStep, 'r-', 'LineWidth', 1.2);
hold on;
plot(td*1e6, discreteStepTustin, 'b--', 'LineWidth', 1.2);
plot(td*1e6, discreteStepZOH, 'k-.', 'LineWidth', 1.2);
title('Step Response Comparison');
xlabel('Time (μs)');
ylabel('Amplitude');
legend('Continuous-Time', 'Discrete-Time (Tustin)', 'Discrete-Time (ZOH)');
grid on;

disp('Discrete-time system matrices (ZOH):');
disp('Ad = '); disp(Ad);
disp('Bd = '); disp(Bd);
disp('Cd = '); disp(Cd);
disp('Dd = '); disp(Dd);

disp('Poles of the discrete-time system (ZOH):');
disp(eig(Ad));
