%% loads data from simulink

uin = u.Data;
yin = y.Data;
xin = x.Data;

% load system matrices
A = Ad; B = Bd; C = Cd; D = Dd; 

% load initial conditions
x0 = [0; 0];

% load covariances
Rn = Rnsim; Rv = Rvsim;


%% State estimation

% xPred computed using x(t+1) = A x(t) + B u(t)
xPred = [x0  zeros(size(xin, 2), size(xin, 1)-1)];
% xKalman computed using the Kalman filter
xKalman = [x0 zeros(size(xin, 2), size(xin, 1)-1)];

% Kalman gain
K = zeros(size(xin, 2), size(xin, 1)); % 2D matrix to store Kalman gain at each time step
% covariance
P = zeros(size(xin, 2), size(xin, 2), size(xin, 1)); % 3D tensor to store covariance at each time step
x1Std = zeros(size(xin, 1), 1);
x2Std = zeros(size(xin, 1), 1);

% high values as high uncertainty at first
P(:, :, 1) = [100 0; 0 100];

for itr = 1:size(xin, 1)-1

    % prediction only
    xPred(:, itr+1) = A * xPred(:, itr) + B * uin(itr, :);

    % Kalman
    Q = A * P(:, :, itr) * A' + Rv; % prediction covariance
    K(:, itr+1) = Q * C' / (C * Q * C' + Rn); % Kalman gain
    P(:, :, itr+1) = (eye(size(A, 1)) - K(:, itr+1) * C) * Q; % updated covariance

    x1Std(itr+1) = sqrt(P(1, 1, itr+1));
    x2Std(itr+1) = sqrt(P(2, 2, itr+1));

    prediction = A * xKalman(:, itr) + B * uin(itr, :);
    measurement = (yin(itr+1) - C*prediction - D*uin(itr+1, :));

    xKalman(:, itr+1) = prediction + K(:, itr+1) * measurement;

end

%% Plots

t = 1:size(xin, 1);

figure('Position', [100, 100, 1200, 800]);
subplot(2, 2, 1);
plot(xin(:,1), 'k', LineWidth=1.5); hold on;
plot(xPred(1, :), 'b--', LineWidth=1.5);
plot(xKalman(1, :), 'r-.', LineWidth=1.5);
fill([t, fliplr(t)], [xKalman(1, :)+x1Std', fliplr(xKalman(1, :)-x1Std')], 'r', ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
title('State Estimation - V_C');
xlabel('Time step');
ylabel('V_C');
legend('True V_C', 'Predicted V_C', 'Kalman Estimated V_C');
xlim([1, size(xin, 1)]);
grid on;

subplot(2, 2, 2);
plot(xin(:, 2), 'k', LineWidth=1.5); hold on;
plot(xPred(2, :), 'b--', LineWidth=1.5);
plot(xKalman(2, :), 'r-.', LineWidth=1.5);
fill([t, fliplr(t)], [xKalman(2, :)+x2Std', fliplr(xKalman(2, :)-x2Std')], 'r', ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
title('State Estimation - I');
xlabel('Time step');
ylabel('I');
legend('True I', 'Predicted I', 'Kalman Estimated I');
xlim([1, size(xin, 1)]);
grid on;

subplot(223);
plot(db(abs((xin(:,1)'-xPred(1,:))./rms(xin(:,1)'))), 'b', LineWidth=1.5); hold on;
plot(db(abs((xin(:,1)'-xKalman(1, :))./rms(xin(:,1)'))), 'r', LineWidth=1.5);
plot(db(x1Std./rms(xin(:,1)')), 'k--', LineWidth=1.5);
title('Normalized State Error - V_C');
xlabel('Time step');
ylabel('V_C error (dB)');
legend('Prediction Error', 'Kalman Estimation Error', 'Normalized Std Dev');
xlim([1, size(xin, 1)]);
grid on;

subplot(224);
plot(db(abs((xin(:,2)'-xPred(2,:))./rms(xin(:,2)'))), 'b', LineWidth=1.5); hold on;
plot(db(abs((xin(:,2)'-xKalman(2, :))./rms(xin(:,2)'))), 'r', LineWidth=1.5);
plot(db(x2Std./rms(xin(:,2)')), 'k--', LineWidth=1.5);
title('Normalized State Error - I');
xlabel('Time step');
ylabel('I error (dB)');
legend('Prediction Error', 'Kalman Estimation Error', 'Normalized Std Dev');
xlim([1, size(xin, 1)]);
grid on;
