NC = 6; % #Cells which equals #BSs
P = 16;
nvar = 1.9905e-08; % Noise Variance ?????
epsilon = 1e-5; % For convergence test.
inner_radius = 500; 
minR_ratio = 0.01;
numIter = 2000;
num_reals = 1000;
alpha_rng = [1, 10];
seed = 1;
fontSize = 15;

% for NU=2:20
%     clear H  in D;
%     fileName = sprintf('channels_for_NU_fixed/Channels%dx%dpower%d.mat', NC, NU, P);
%     load(fileName,'H', 'in', 'D'); 
%     last_hope
% end
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
            file_name = sprintf('WMMSE_for_powers_fixed/WMMSE_%dx%dpower%dalpha%d.mat', NC, NU,P, alpha_idx);
            load(file_name, 'Powers', 'conv', 'R_sums', 'Rmax_sums', 'WR_sums', 'WRmax_sums', 'tdma_rates');
            RR(NU) = mean(R_sums);
            RR_max(NU) = mean(Rmax_sums);
            WRR(NU) = mean(WR_sums);
            WRR_max(NU) = mean(WRmax_sums);
            tdma(NU) = mean(tdma_rates);
            % PP(:, :, P) = mean(Powers, 3);
    %         convv(NU) = mean(conv); 
        end
        if (alpha_idx == 1)
            plot(2:20, RR(2:20), 'm<-', 'linewidth',1);
        elseif (alpha_idx == 10)
            plot(2:20, RR(2:20), 'bo-', 'linewidth',1);
        end
    end

    plot(2:20, RR_max(2:20), 'ks-.', 'linewidth',1);
    plot(2:20, tdma(2:20), 'gd--', 'linewidth',1);
end
legend('distance-based alpha', 'uniformly distributed alpha', 'uniform power allocation', 'OMA', 'FontSize', fontSize);
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


    
    
    
    
    
    