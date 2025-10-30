clear; close all; clc;

%% generation

testSize = 1000;
noisePow = 0.2;

inFFT = [randn(1, 10), zeros(1, testSize-10)];
in = real(ifft(inFFT));
in = in / rms(in);
noise = sqrt(noisePow) * randn(1, testSize);

%% hyperparameters
A = 1;
l_vals = [1 10 50 200];

% squared exponential function
tStar = 1:0.5:testSize;
t = 1:testSize;

figure;

for i = 1:length(l_vals)
    l = l_vals(i);

    Ctranspose = A * exp(-(tStar(:) - t(:).').^2/l^2);
    B = A * exp(-(t(:) - t(:).').^2/l^2) + noisePow*eye(length(t));
    y = (in+noise)';
    
    f_est = Ctranspose * (B\y);
    
    subplot(2,2,i);
    plot(t, y, '*');
    hold on;
    plot(tStar, f_est, LineWidth=2.5);
    plot(t, in, LineWidth=1.5);
    legend("Noisy signal", "estimated fit", "True signal");
    title("l = " + num2str(l));

end

