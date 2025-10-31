clear; close all; clc;


%% 1. Parameters

N = 64e3;
fs = 8e3;
g_list = [0.8 0.9 0.95 0.99 0.995 0.999];

f0 = fs/N; % frequency resolution

nb = 1;
na = 2;

n_start = max(na, nb); % first index where estimation can start
n_LLS = 500; % arbitrarily chosen as 50 * size(theta)

if n_LLS < n_start + 1
    error('n_LLS must be greater than max(na, nb) + 1');
end

%% 2. Load data

%load();
% u = input signal (N x 1)
% y = output signal (N x 1)
load("results/test4.mat");
u = Su;
y = Sy;


%% 3. Recursive estimation

% initialization with an LLS estimate over a few of the first samples
% y_lls = H * theta
H = zeros(n_LLS - n_start, na + nb + 1);
for uCol = 0:nb
    H(:, uCol + 1) = u((n_start + 1 - uCol:n_LLS - uCol)); % shifting input
end
for yCol = 1:na
    H(:, nb + 1 + yCol) = -y((n_start + 1 - yCol:n_LLS - yCol)); % shifting output
end
y_lls = y((n_start + 1:n_LLS));

theta = H\y_lls;
P = inv(H'*H);

figure;

for gIndex = 1:length(g_list)
    g = g_list(gIndex);

    theta_sav = zeros(na + nb + 1, N - n_start + 2);
    P_diag_sav = zeros(na + nb + 1, N - n_start + 2);

    for itr = n_LLS+1:N

        % saving previous values
        theta_sav(:, itr - n_start + 1) = theta;
        P_diag_sav(:, itr - n_start + 1) = diag(P);

        K_itr = [u(itr - (0:nb))', -y(itr - (1:na))'];

        theta = theta - (P * (K_itr') * (K_itr * theta - y(itr))) / (g + K_itr * P * (K_itr'));
        P = 1/g * (P - (P * (K_itr') * K_itr * P) / (g + K_itr * P * (K_itr')));

    end

    %% 4. Plot results with shaded variance areas

    % time vector corresponding to columns of theta_sav / P_diag_sav
    timeTheta = ((0:size(theta_sav,2)-1) + n_start - 1) / fs;

    subplot(2,3,gIndex);
    colors = lines(nb + na + 1);
    labels = cell(nb + na + 1, 1);
    for i = 0:nb
        labels{i+1} = sprintf('b%d', i);
    end
    for j = 1:na
        labels{nb + 1 + j} = sprintf('a%d', j);
    end

    % plot parameter estimates with shaded std dev
    hold on;
    hLine = gobjects(nb + na + 1, 1);
    for idx = 1:(nb + na + 1)
        yVal = theta_sav(idx, :);
        stdv = sqrt(P_diag_sav(idx, :));
        upper = yVal + stdv;
        lower = yVal - stdv;
        t = timeTheta;
        % shaded area (exclude from legend)
        fill([t, fliplr(t)], [upper, fliplr(lower)], colors(idx,:), ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        % plot mean line and save handle for legend
        hLine(idx) = plot(t, yVal, 'Color', colors(idx,:), 'LineWidth', 1.5);
    end
    xlabel('Time (s)');
    ylabel('Estimate value');
    ylim([-2 1.5]);
    title('RLS, g = ' + string(g));
    legend(hLine, labels, 'Interpreter', 'none');
    hold off;

end
