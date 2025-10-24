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

%% 2. Load data

%load();
% u = input signal (N x 1)
% y = output signal (N x 1)


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

    K_itr = [u(itr - (0:nb)), -y(itr - (1:na))];

    theta = theta - (P * (K_itr') * (K_itr * theta - y(itr))) / (g + K_itr * P * (K_itr'));
    P = 1/g * (P - (P * (K_itr') * K_itr * P) / (g + K_itr * P * (K_itr')));

end


