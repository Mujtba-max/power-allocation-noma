% Calculating TDMA rates
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

close all;
for P=0:P_max_idx
    clear H  in D;
    fileName = sprintf('channels_for_powers_fixed/Channels%dx%dpower%d.mat', NC, NU, P);
    load(fileName,'H', 'in', 'D'); 
    last_hope_tdma;
end

clear RR RR_max tdma WRR WRR_max convv;

tdma = zeros(10, P_max_idx+1);
for alpha_idx = alpha_rng
    if (mod(alpha_idx,2)) == 1
        for P = 1:P_max_idx+1
            file_name = sprintf('OMA_for_powers/OMA_%dx%dpower%dalpha%d.mat', NC, NU,P-1, alpha_idx);
            load(file_name,'tdma_rates');
            tdma(alpha_idx, P) = mean(tdma_rates);
            % PP(:, :, P) = mean(Powers, 3);
    %         convv(alpha_idx, P) = mean(conv); 
        end
    end

%     plot(0:20, RR);
%     hold on;
    % plot(0:20, convv);

end

figure; hold on; grid on;

plot(0:20, tdma(1, :), 'gd--', 'linewidth',2);