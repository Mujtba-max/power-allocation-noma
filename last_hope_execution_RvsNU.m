alpha_rng_length = 10; % number of alpha values to compute.
P = 16;
nvar = 1.9905e-08; % Noise Variance ?????
epsilon = 1e-5; % For convergence test.
inner_radius = 500; 
minR_ratio = 0.01;
numIter = 2000;
num_reals = 1000;
alpha_rng = [1, alpha_rng_length];
seed = 1;
fontSize = 20;

clear RR RR_max WRR WRR_max convv;
figure;
for NC=4:6
    i = NC-3;
    subplot(1,3,i);
    t = sprintf('Number of cells = %d', NC);
    title(t, 'FontSize', fontSize);
    xlabel('Number of users', 'FontSize', fontSize);
    ylabel('Sum rate (bits/s/Hz)', 'FontSize', fontSize);
    xlim([2,20]);
    ylim([0,25]);
    hold on; grid on;
    for alpha_idx = alpha_rng
        for NU = 2:20
%             clear H in D; % clear the variables H, in, and D.
%             fileName = sprintf('channels_for_NU/Channels%dx%dpower%d.mat', NC, NU, P);
%             load(fileName,'H', 'in', 'D');
%             current = sprintf('NC = %d, NU = %d, alpha = %d', NC, NU, alpha_idx);
%             disp(current)
%             last_hope

            fileName = sprintf('WMMSE_for_NU/WMMSE_%dx%dpower%dalpha%dabs.mat', NC, NU, P, alpha_idx);
            load(fileName, 'R_sums', 'Rmax_sums', 'WR_sums', 'WRmax_sums', 'tdma_rates');
            
            RR(NU) = mean(R_sums);
            RR_max(NU) = mean(Rmax_sums);
            WRR(NU) = mean(WR_sums);
            WRR_max(NU) = mean(WRmax_sums);
            tdma(NU) = mean(tdma_rates);

            % PP(:, :, P) = mean(Powers, 3);
    %         convv(NU) = mean(conv); 
        end
        if (alpha_idx == 1)
            plot(2:20, RR(2:20), 'bo-', 'linewidth',2);
        elseif (alpha_idx == 10)
            plot(2:20, RR(2:20), 'm<-', 'linewidth',2);
        end
    end
% 
    plot(2:20, RR_max(2:20), 'ks-.', 'linewidth',2);
    plot(2:20, tdma(2:20), 'gd--', 'linewidth',2);
end

legend('distance-based-alpha WMMSE', 'uniformly-distributed-alpha WMMSE', 'uniform power allocation', 'OMA', 'FontSize', 20);

% plot(2:20, convv);

% NC = 2;
% 
% for NU = 2:8
%     
%     file_name = sprintf('channels_for_NU_last/Channels%dx%d.mat', NC, NU);
%     load(file_name);
%     NU
%     sum(in, 1)
% end






