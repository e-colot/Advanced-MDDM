close all;
clear;

%% Task 1.1 - SIGNAL GENERATION
% Parameters
fs = 8000;        % Sampling rate (Hz)
N = 64000;        % Number of samples
sigma = 1;        % Standard deviation

% Generate white Gaussian noise signal
u0 = sigma * randn(N, 1);

% Time vector (optional)
t = (0:N-1)'/fs;

% % Plot the signal (optional)
% figure;
% plot(t(1:1000), u0(1:1000)); % Plot only first 1000 samples for clarity
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('White Gaussian Noise Signal u_0(t)');
% grid on;

%% Load Output signal

system = 1;

switch system
    case 1
        load("Low_LTI.mat");
    case 2
        load("High_LTI.mat");
    case 3
        load("Slow_Shift.mat");
    case 4
        load("Fast_Shift.mat");
end
Output = [Su_Time, Su, Sy_Time, Sy];


%% Task 1.2 - Regression Matrix Total

y_0_full = Output(3:N,4);
y_1_full = [Output(2:N-1,4)];
y_2_full = [Output(1:N-2,4)];
u_0_full = Output(3:N,2);
u_1_full = [Output(2:N-1,2)];
K_tot = [-y_1_full, -y_2_full, u_0_full, u_1_full];

disp(['K(N) = ', num2str(K_tot(N-2,:))]);

Theta_True = K_tot \  y_0_full;

disp(['True parameters :', num2str(Theta_True')]);

P_true = inv(K_tot'*K_tot);
Var_LLS = [P_true(1,1); P_true(2,2); P_true(3,3); P_true(4,4)];

disp(['True variance :', num2str(Var_LLS')]);

%% Task 1.3.1 - Regression Matrix

N_iteration1 = 100;
y_0 = Output(3:N_iteration1,4);
y_1 = [Output(2:N_iteration1-1,4)];
y_2 = [Output(1:N_iteration1-2,4)];
u_0 = Output(3:N_iteration1,2);
u_1 = [Output(2:N_iteration1-1,2)];
K = [-y_1, -y_2, u_0, u_1];

disp(['K(N) = ', num2str(K(N_iteration1-2,:))]);

Theta_LLS = K \ y_0;

disp(['Parameters after 1st iteration :', num2str(Theta_LLS')]);

%% Task 1.3.2 \ Task 1.4 \ Task 1.5 - Recursive parameters

g_matrix = [0.7, 0.8, 0.9, 0.95, 0.99, 0.9999];
size_forg = size(g_matrix,2);
Theta = zeros(4,N-N_iteration1,size_forg);
Var_para = zeros(4,N-N_iteration1,size_forg);
for forg = 1 : 1 : size_forg
    
    g = g_matrix(forg);
    Theta(:,1,forg) = Theta_LLS;
    P_prev = inv(K'*K);
    Var_para(:,1,forg) = [P_prev(1,1); P_prev(2,2); P_prev(3,3); P_prev(4,4)];
    
    for t = 1 : 1 : N-N_iteration1
        K_Nrow = K_tot(N_iteration1-2+t,:);
        y_Nrow = Output(N_iteration1+t,4);
        Theta(:,t+1,forg) = Theta(:,t,forg)-((P_prev*K_Nrow')*(K_Nrow*Theta(:,t,forg)-y_Nrow))/(g+K_Nrow*P_prev*K_Nrow');
        P =(1/g)*(P_prev-((P_prev*K_Nrow'*K_Nrow*P_prev)/(g+(K_Nrow*P_prev*K_Nrow'))));
        P_prev = P;
        Var_para(:,t+1,forg) = [P_prev(1,1); P_prev(2,2); P_prev(3,3); P_prev(4,4)];
    end
end

disp(['Parameters after recursivity:', num2str(Theta(:,end,size_forg)')]);

%% Task 1.6 - Visualization

visu = 1;

if visu == 1

    % Uncertainty calculation
    
    Std_para = sqrt(Var_para);
    lim_inf = Theta - Std_para;
    lim_sup = Theta + Std_para;
    sample = 1 : N - N_iteration1 + 1;
    
    colors = {'magenta','red','blue','green'};
    param_names = {'a1','a2','b0','b1'};
    
    figure;
    for k = 1:size_forg
        subplot(2, 3, k); % Adjust layout automatically for up to 6 forgetting factors
        hold on;
        
        % Plot each parameter with its uncertainty band
        for i = 1:4
            x = [sample, fliplr(sample)];
            y = [lim_sup(i,:,k), fliplr(lim_inf(i,:,k))];
            fill(x, y, colors{i}, 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility', 'off');
            plot(sample, Theta(i,:,k), 'Color', colors{i}, 'LineWidth', 1.5);
        end
        
        ylim ([-6 6]);
        xlabel('Time');
        ylabel('\theta estimation');
        title(['Parameter evolution (g = ', num2str(g_matrix(k)), ')']);
        grid on;
        legend(param_names, 'Location', 'best');
        hold off;
    end
end

%% True FRF for LTI system

if system == 1 | system == 2

    Nfft = 3200;
    
    uAperiodic = reshape(Output(1:N,2), Nfft, 20);
    yAperiodic = reshape(Output(1:N,4), Nfft, 20);

    uDFT = fft(uAperiodic);
    yDFT = fft(yAperiodic);
    
    freqaxis = (0:Nfft/2-1)*fs/Nfft;
    
    FRF = mean(yDFT.*conj(uDFT),2) ./ mean(uDFT.*conj(uDFT),2);
    FRF = FRF(1:Nfft/2);
    
    figure;
    plot(freqaxis,db(FRF),'Linewidth',1);
    xlabel('Frequency axis (Hz)');
    ylabel('Magnitude')

end