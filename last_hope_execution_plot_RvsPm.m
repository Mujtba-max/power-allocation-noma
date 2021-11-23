alpha_rng_length = 10; % number of alpha values to compute.
NC = 4; % #Cells which equals #BSs
NU = 10; % #USers in each cell.
nvar = 1.9905e-08; % Noise Variance ?????
epsilon = 1e-5; % For convergence test.
inner_radius = 500; 
minR_ratio = 0.01;
numIter = 2000;
num_reals = 1000;
alpha_rng = [1, alpha_rng_length];
P_max_idx = 20;
seed = 1;

clear RR RR_max tdma WRR WRR_max convv;

%%% I replaced every 10 with alpha_rng_length down here ! Did I get it wrong?
RR = zeros(alpha_rng_length, P_max_idx+1); % The average rate for the WMMSE NOMA
RR_max = zeros(alpha_rng_length, P_max_idx+1); % The average rate for the fixed NOMA, and so on.
WRR = zeros(alpha_rng_length, P_max_idx+1);
WRR_max = zeros(alpha_rng_length, P_max_idx+1);
tdma = zeros(alpha_rng_length, P_max_idx+1);
fdma = zeros(alpha_rng_length, P_max_idx+1);
convv = zeros(alpha_rng_length, P_max_idx+1); 
for alpha_idx = alpha_rng
    for P = 1:P_max_idx+1
        file_name = sprintf('WMMSE_for_powers/WMMSE_%dx%dpower%dalpha%d.mat', NC, NU,P-1, alpha_idx);
        
        % check if the file exists, if not, execute lase_hope and save the parameters into the file.
        if isfile(file_name)
            % File exists.
            load(file_name, 'Powers', 'conv', 'R_sums', 'Rmax_sums', 'WR_sums', 'WRmax_sums', 'tdma_rates');
        else
           clear H  in D;
           fileName = sprintf('channels_for_powers/Channels%dx%dpower%d.mat', NC, NU, P-1);
           load(fileName,'H', 'in', 'D'); 
           last_hope
        end

        % here we are taking the mean of the rates?
        RR(alpha_idx,P) = mean(R_sums);
        RR_max(alpha_idx, P) = mean(Rmax_sums);
        WRR(alpha_idx, P) = mean(WR_sums);
        WRR_max(alpha_idx, P) = mean(WRmax_sums);
        tdma(alpha_idx, P) = mean(tdma_rates);
        if (mod(alpha_idx,2) == 1)
            file_name = sprintf('OMA_for_powers/OMA_%dx%dpower%dalpha%d.mat', NC, NU,P-1, alpha_idx);
            load(file_name,'tdma_rates'); % ???
            fdma(alpha_idx, P) = mean(tdma_rates);
        end
        
            %%%%% Ali from the future: What are these comments below???

        % PP(:, :, P) = mean(Powers, 3);
%         convv(alpha_idx, P) = mean(conv); 
    end
        
            %%%%% Ali from the future: What are these comments below???

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


function [] = execute_last_hope()
    clear H in D; % clear the variables H, in, and D.
    fileName = sprintf('channels_for_NU_fixed/Channels%dx%dpower%d.mat', NC, NU, P);
    load(fileName,'H', 'in', 'D'); 
    last_hope
end

              %%%%% Ali from the future: What are these comments below???

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
