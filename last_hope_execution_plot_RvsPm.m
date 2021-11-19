NC = 4; % #Cells which equals #BSs
NU = 10; % #USers in each cell.
nvar = 1.9905e-08; % Noise Variance ?????
epsilon = 1e-5; % For convergence test.
inner_radius = 500; 
minR_ratio = 0.01;
numIter = 2000;
num_reals = 1000;
alpha_rng = [1, 10];
P_max_idx = 20;
seed = 1;

% 
% for P=0:P_max_idx
%     clear H  in D;
%     fileName = sprintf('channels_for_powers_fixed/Channels%dx%dpower%d.mat', NC, NU, P);
%     load(fileName,'H', 'in', 'D'); 
%     last_hope
% end

clear RR RR_max tdma WRR WRR_max convv;
RR = zeros(10, P_max_idx+1);
RR_max = zeros(10, P_max_idx+1);
WRR = zeros(10, P_max_idx+1);
WRR_max = zeros(10, P_max_idx+1);
tdma = zeros(10, P_max_idx+1);
fdma = zeros(10, P_max_idx+1);
convv = zeros(10, P_max_idx+1); 
for alpha_idx = alpha_rng
    for P = 1:P_max_idx+1
        file_name = sprintf('WMMSE_for_powers_fixed/WMMSE_%dx%dpower%dalpha%d.mat', NC, NU,P-1, alpha_idx);
        load(file_name, 'Powers', 'conv', 'R_sums', 'Rmax_sums', 'WR_sums', 'WRmax_sums', 'tdma_rates');

        RR(alpha_idx,P) = mean(R_sums);
        RR_max(alpha_idx, P) = mean(Rmax_sums);
        WRR(alpha_idx, P) = mean(WR_sums);
        WRR_max(alpha_idx, P) = mean(WRmax_sums);
        tdma(alpha_idx, P) = mean(tdma_rates);
        if (mod(alpha_idx,2) == 1)
            file_name = sprintf('OMA_for_powers/OMA_%dx%dpower%dalpha%d.mat', NC, NU,P-1, alpha_idx);
            load(file_name,'tdma_rates');
            fdma(alpha_idx, P) = mean(tdma_rates);
        end
        % PP(:, :, P) = mean(Powers, 3);
%         convv(alpha_idx, P) = mean(conv); 
    end

%     plot(0:20, RR);
%     hold on;
    % plot(0:20, convv);

end

figure; hold on; grid on;
plot(0:20, RR(10, :),'m<-', 'linewidth',2);
plot(0:20, RR(1, :), 'bo-', 'linewidth',2);
% plot(0:20, RR(2, :), 'bo--', 'linewidth',2);
% plot(0:20, RR(3, :), 'r*-', 'linewidth',2);
% plot(0:20, RR(4, :), 'r*--', 'linewidth',2);
% plot(0:20, RR(5, :), 'g+-', 'linewidth',2);
% plot(0:20, RR(6, :), 'g+--', 'linewidth',2);
% plot(0:20, RR(7, :), 'cv-', 'linewidth',2);
% plot(0:20, RR(8, :), 'cv--', 'linewidth',2);
% plot(0:20, RR(9, :), 'm<-', 'linewidth',2);


plot(0:20, RR_max(1, :), 'ks-.', 'linewidth',2);
plot(0:20, fdma(1, :), 'gd--', 'linewidth',2);

xlabel('Maximum power for each BS (dBW)', 'FontSize', 15);
ylabel('Sum rate (bits/s/Hz)', 'FontSize', 15);

legend('uniformly-distributed-alpha WMMSE', 'distance-based-alpha WMMSE', 'uniform power allocation', 'OMA', 'FontSize', 15);



% lgnd = cell(7, 3);
% for alpha_idx = 1:2:7
%     lgnd(alpha_idx, 1) = {sprintf('WMMSE NOMA -a=%d', alpha_idx)};
% %     lgnd(alpha_idx, 2) = {sprintf('Uniform NOMA-a=%d', alpha_idx)};
% %     lgnd(alpha_idx, 3) = {sprintf('OMA-a=%d', alpha_idx)};
% end
% 
% lgnd(1, 2) = {sprintf('Uniform NOMA-a=%d', alpha_idx)};
% lgnd(1, 3) = {sprintf('OMA-a=%d', alpha_idx)};
% % , lgnd{7, 1}, lgnd{7, 2}, lgnd{7, 3}
% % , lgnd{5, 1}, lgnd{5, 2}, lgnd{5, 3}
% legend(lgnd{1, 1}, lgnd{3, 1}, lgnd{5, 1}, lgnd{7, 1}, lgnd{1, 2}, lgnd{1, 3});
