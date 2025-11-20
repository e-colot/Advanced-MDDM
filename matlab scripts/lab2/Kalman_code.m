%% KALMAN FILTER

% Amaury Arico
% MA2-IRELE
% Bruafec VUB ULB
% 20251118

%% -----------------------------------------------------

close all;
% Initialisation of the algorithm
% case = 0 is the state space matrix computation
% case = 1 is the kalman filter application in discrete time

run_case = 1;


%% Task 2.1 - Continuous-time state-space model
if (run_case == 0)     
    % Parameters
    L0 = 5e-6;     % Inductance [H]
    R0 = 10;       % Resistance [Ohm]
    C0 = 3e-9;     % Capacitance [F]
    
    % Continuous-time matrices
    Ac = [0, 1/C0;
         -1/L0, -R0/L0];
    
    Bc = [0; 1/L0];
    Cc = [1 0];
    Dc = 0;
    
    % Display matrices
    disp('Continuous-time matrices:')
    Ac, Bc, Cc, Dc;
    
    % Create state-space system
    sys_c = ss(Ac, Bc, Cc, Dc);
    
    %% Task 2.2 - Discretization
    
    % Choose an appropriate sampling time
    % Rule of thumb: Ts << (1 / (10 * natural frequency))
    wn = 1 / sqrt(L0 * C0);   % Natural frequency (rad/s)
    Ts = 2*pi / (10* wn);       % Sampling time (choose high sampling rate)
    
    disp(['Chosen sampling time Ts = ', num2str(Ts), ' s'])
    
    % Discretize with zero-order hold
    sys_d = c2d(sys_c, Ts, 'zoh');
    sys_d_tust = c2d(sys_c, Ts, 'tustin');
    
    % Extract discrete matrices
    [Ad, Bd, Cd, Dd] = ssdata(sys_d);
    
    disp('Discrete-time matrices:')
    Ad, Bd, Cd, Dd;
    
    %% Compare step responses
    figure;
    step(sys_c, 'b', sys_d_tust, 'r--', sys_d, 'k--');
    legend('Continuous-time', 'Discrete-time (Tustin)', 'Discrete-time (ZOH)');
    title('Step Response Comparison');
    grid on;
    
    %% Stability check
    eig_Ac = eig(Ac);
    eig_Ad = eig(Ad);
    
    disp('Eigenvalues of Ac (continuous):'), disp(eig_Ac)
    disp('Eigenvalues of Ad (discrete):'), disp(eig_Ad)
    
    if all(real(eig_Ac) < 0)
        disp('The continuous-time system is stable.');
    else
        disp('The continuous-time system is unstable.');
    end
    
    if all(abs(eig_Ad) < 1)
        disp('The discrete-time system is stable.');
    else
        disp('The discrete-time system is unstable.');
    end
    
    %% Observability check
    O = obsv(Ac, Cc);
    if rank(O) == size(Ac,1)
        disp('The system is observable.');
    else
        disp('The system is NOT observable.');
    end
    
    %% Rv and Rn
    
    Rv = [0.025 0 ; 0 1e-6];
    Rn = 0.0025;
    
    stop_time = 199*Ts;
end

%%  Kalman filtering loop

if (run_case == 1)
    time_axis = (1:200)';
    y_meas = y.data;
    x_meas = x.data;
    u_exc = u.data;
    Rv = [0.025 0 ; 0 1e-6];
    Rn = 0.0025;

    figure;
    subplot(121);
    plot(time_axis, y_meas);
    title("Voltage Output measurements");
    xlabel('Time index');
    ylabel('Voltage (V)');
    subplot(122);
    plot(time_axis, u_exc);
    title("Signal input");
    xlabel('Time index');
    ylabel('Current (V)');
    
    %% Observer - 1 : Time prediction / 2 : Kalman
    
    % Initialization
    P_0 = 100*eye(2);
    x_0 = [0 ; 0];
    x_est = zeros(2,200);
    x_kalman = zeros(2,200);
    P = zeros(2,2,200);
    Q = zeros(2,2,200);
    K = zeros(2,200);
    P(:,:,1) = P_0;
    A = Ad;
    C = Cd;
    B = Bd;
    D = Dd;
    %A=2*Ad;                            %Case 2
    x_0 = [5 ; 5];               %Case 3
    %Rn = 100*Rn;                       %Case 4
    %Rv = 100*Rv;                       %Case 5

    x_est(:,1) = x_0;
    x_kalman(:,1) = x_0;
    
    % Loop of Kalman

    for i= 1:1:199
        x_est(:,i+1) = A*x_est(:,i)+B*u_exc(i);
        Q(:,:,i+1) = A*P(:,:,i)*A'+Rv;
        K(:,i+1) = Q(:,:,i+1)*C'/(C*Q(:,:,i+1)*C'+Rn);
        P(:,:,i+1)=(eye(2)-K(:,i+1)*C)*Q(:,:,i+1);
        time_pred = A*x_kalman(:,i)+B*u_exc(i);
        x_kalman(:,i+1)=time_pred + K(:,i+1)*(y_meas(i+1)-C*time_pred-D*u_exc(i));
    end
    kal_std = sqrt([P(1,1,:);P(2,2,:)]);
    kal_std = squeeze(kal_std);
    pred_std = sqrt([Q(1,1,:);Q(2,2,:)]);
    pred_std = squeeze(pred_std);
    norm_er_kal = (x_kalman'-x_meas)./(rms(x_meas));
    norm_er_pred = (x_est'-x_meas)./(rms(x_meas));
    norm_std_kal = kal_std'./(rms(x_meas));

    uncertainty_kal_low = x_kalman-kal_std;
    uncertainty_kal_high = x_kalman+kal_std;

    uncertainty_pred_low = x_est- pred_std;
    uncertainty_pred_high = x_est+ pred_std;

    t = 1:200;
    figure;
    subplot(121);
    fill([t(2:end) fliplr(t(2:end))],[uncertainty_kal_low(1,2:end) fliplr(uncertainty_kal_high(1,2:end))],[0.4 0.6 1],'EdgeColor','none','FaceAlpha',0.5);
    hold on
    fill([t(5:end) fliplr(t(5:end))],[uncertainty_pred_low(1,5:end) fliplr(uncertainty_pred_high(1,5:end))],[1 0.6 1],'EdgeColor','none','FaceAlpha',0.5);
    plot(time_axis, x_meas(:,1),'black--','LineWidth',1');
    plot(time_axis, x_est(1,:)','r');
    plot(time_axis, x_kalman(1,:)','b');
    legend("Kal Uncert","Pred Uncert","True Volt", "Predict Volt", "Kalman Volt");
    title("Voltage State measurements");
    xlabel('Time index');
    ylabel('Voltage (V)');

    subplot(122);
    fill([t(4:end) fliplr(t(4:end))],[uncertainty_kal_low(2,4:end) fliplr(uncertainty_kal_high(2,4:end))],[0.4 0.6 1],'EdgeColor','none','FaceAlpha',0.5);
    hold on
    fill([t(5:end) fliplr(t(5:end))],[uncertainty_pred_low(2,5:end) fliplr(uncertainty_pred_high(2,5:end))],[1 0.6 1],'EdgeColor','none','FaceAlpha',0.5);
    plot(time_axis, x_meas(:,2),'black--','LineWidth',1');
    hold on
    plot(time_axis, x_est(2,:)','r');
    hold on
    plot(time_axis, x_kalman(2,:)','b');
    legend("Kal Uncert","Pred Uncert","True Cur", "Predict Cur", "Kalman Cur");
    xlabel('Time index');
    ylabel('Current (A)')
    title("Current State measurements");
    

    figure;
    subplot(121);
    plot(time_axis, db(norm_er_pred(:,1)));
    hold on
    plot(time_axis, db(norm_er_kal(:,1)));
    plot(time_axis, db(norm_std_kal(:,1)), 'black--');
    legend("Pred error", "Kalman error","Kal. std deviation");
    title("Normalized error on the Voltage Vc");
    xlabel('Time index');
    ylabel('Norm. Error (dB)')
    ylim([-80 20]);
    subplot(122);
    plot(time_axis, db(norm_er_pred(:,2)));
    hold on
    plot(time_axis, db(norm_er_kal(:,2)));
    plot(time_axis, db(norm_std_kal(:,2)), 'black--');
    legend("Pred error", "Kalman error","Kal. std deviation");
    title("Normalized error on the Current");
    xlabel('Time index');
    ylabel('Norm. Error (dB)')
    ylim([-80 20]);


    
end
