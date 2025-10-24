clear; close all; clc;


%% 1. Parameters

N = 64e3;
fs = 8e3;
g = 0.99;

f0 = fs/N;

nb = 1;
na = 2;

sampleForLLS = 200; % arbitrarily chosen as 50 * (nb + na)
if sampleForLLS < max(na, nb) + 1
    error('sampleForLLS must be greater than max(na, nb) + 1');
end

%% 2. Simulate system

u = randn(N,1); 

sys1A = [1, -1.5, 0.7];
sys1B = [0.5, 0.3];

y1 = filter(sys1B, sys1A, u(1:N/2));

sys2A = [1, -0.8, 0.3];
sys2B = [0.4, 0.2];

y2 = filter(sys2B, sys2A, u(N/2+1:end));

y = [y1; y2] + 0.05 * randn(N,1);


%% 3. Recursive estimation

% initialization with an LLS estimate over a few of the first samples
% y_lls = H * theta
H = zeros(sampleForLLS - max(na, nb), na + nb + 1);
for uCol = 0:nb
    H(:, uCol + 1) = u((max(na, nb) + 1 - uCol:sampleForLLS - uCol));
end
for yCol = 1:na
    H(:, nb + 1 + yCol) = -y((max(na, nb) + 1 - yCol:sampleForLLS - yCol));
end
y_lls = y((max(na, nb) + 1:sampleForLLS));

theta = H\y_lls;
P = inv(H'*H);

theta_sav = zeros(na + nb + 1, N - max(na, nb) + 2);
P_diag_sav = zeros(na + nb + 1, N - max(na, nb) + 2);

for itr = sampleForLLS+1:N

    % saving previous values
    theta_sav(:, itr - max(na, nb) + 1) = theta;
    P_diag_sav(:, itr - max(na, nb) + 1) = diag(P);

    K_itr = [u(itr - (0:nb))', -y(itr - (1:na))'];

    theta = theta - (P * (K_itr') * (K_itr * theta - y(itr))) / (g + K_itr * P * (K_itr'));
    P = 1/g * (P - (P * (K_itr') * K_itr * P) / (g + K_itr * P * (K_itr')));

end

%% 4. Plot results

timeVec = (0:N-1)/fs;
figure;
subplot(311);
plot(timeVec, theta_sav(1, :), 'LineWidth', 1.5);
hold on;
for idx = 2:nb+na+1
    plot(timeVec, theta_sav(idx, :), 'LineWidth', 1.5);
end
xlabel('Time (s)');
ylabel('Estimate value');
title('Recursive estimates of the parameters');
% build legend entries: b0..bnb then a1..ana
labels = cell(nb + na + 1, 1);
for i = 0:nb
    labels{i+1} = sprintf('b%d', i);
end
for j = 1:na
    labels{nb + 1 + j} = sprintf('a%d', j);
end
legend(labels, 'Interpreter', 'none', 'Location', 'best');

subplot(312);
plot(timeVec, P_diag_sav(1, :), 'LineWidth', 1.5);
hold on;
for idx = 2:nb+na+1
    plot(timeVec, P_diag_sav(idx, :), 'LineWidth', 1.5);
end
xlabel('Time (s)');
ylabel('Variance estimate');
title('Evolution of the variance estimates');
% build legend entries: b0..bnb then a1..ana
labels = cell(nb + na + 1, 1);
for i = 0:nb
    labels{i+1} = sprintf('b%d', i);
end
for j = 1:na
    labels{nb + 1 + j} = sprintf('a%d', j);
end
legend(labels, 'Interpreter', 'none', 'Location', 'best');


subplot(313);
plot(timeVec(2:end), db(sum((theta_sav(:, 2:end) - theta_sav(:, 1:end-1)).^2, 1)), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Norm of parameter change');
title('Convergence of the parameter estimates');

