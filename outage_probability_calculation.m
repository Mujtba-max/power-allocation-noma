% Calculating Outage Probabilities 
NC = 4; % #Cells which equals #BSs
NU = 10; % #USers in each cell.
nvar = 1.9905e-08; % Noise Variance ?????
epsilon = 1e-5; % For convergence test.
inner_radius = 500; 
minR_ratio = 0.01;
numIter = 2000;
num_reals = 1000;
alpha_rng = 1:10;
P_max_idx = 20;
seed = 1;
% P = 16;

outage = zeros(12,P_max_idx+1);
outage_uniform = zeros(1, P_max_idx+1);

for P = 0:P_max_idx
    if mod(P, 4) == 0
        disp(P);
    end
    clear H  in D;
    fileName = sprintf('channels_for_powers_fixed/Channels%dx%dpower%d.mat', NC, NU, P);
    load(fileName,'H', 'in', 'D'); 

    P_linear = 10^(P/10);
    P_uniform = sqrt(P_linear/NU) *ones(NU, NC);
    
    h = zeros(NU,NC,NC,num_reals);
    d = zeros(NU,NC,NC,num_reals);

    for c=1:NC
        num_cell_reals = sum(in(:,c),1);
        cell_indecies = find(in(:,c)==1);
        if(num_cell_reals < num_reals)
            i = 0;
            f = H(:,c,:,cell_indecies);
            dis = D(:,c,:,cell_indecies);
            while (i*num_cell_reals < num_reals)
                start_point = i*num_cell_reals+1;
                end_point = min((i+1)*num_cell_reals, num_reals);
                h(:,c,:,start_point:end_point) = f(:,:,:,1:end_point-start_point+1);
                d(:,c,:,start_point:end_point) = dis(:,:,:,1:end_point-start_point+1);
                i = i + 1;
            end
        else
            h(:,c,:,:) = H(:, c, :, cell_indecies(1:num_reals)); d(:,c,:,:) = D(:, c, :, cell_indecies(1:num_reals));
        end
    end

    for alpha_idx=[1,10]
        file_name = sprintf('WMMSE_for_powers_fixed/WMMSE_%dx%dpower%dalpha%d.mat', NC, NU,P, alpha_idx);
        load(file_name, 'Powers');

        for s=1:num_reals
            HH = abs(h(:,:,:,s));
            dd = d(:,:,:,s);
            
            R = rate(NC, NU, HH, Powers(:,:,s), nvar);
            outage(alpha_idx, P+1) = outage(alpha_idx, P+1) + sum(R<0.005, 'all');
            if alpha_idx == alpha_rng(1)
                R_uniform = rate(NC, NU, HH, P_uniform, nvar);
                outage_uniform(P+1) = outage_uniform(P+1) + sum(R_uniform<0.5, 'all');
            end
        end
        outage(alpha_idx, P+1) = outage(alpha_idx, P+1)/(NC*NU*num_reals);
        outage_uniform(P+1) = outage_uniform(P+1)/(NC*NU*num_reals);
    end
end


figure; hold on; grid on;
plot(0:20, outage(1, :), 'bo-', 'linewidth',2);
% plot(0:20, outage(2, :), 'bo--', 'linewidth',2);
% plot(0:20, outage(3, :), 'r*-', 'linewidth',2);
% plot(0:20, outage(4, :), 'r*--', 'linewidth',2);
% plot(0:20, outage(5, :), 'g+-', 'linewidth',2);
% plot(0:20, outage(6, :), 'g+--', 'linewidth',2);
% plot(0:20, outage(7, :), 'cv-', 'linewidth',2);
% plot(0:20, outage(8, :), 'cv--', 'linewidth',2);
% plot(0:20, outage(9, :), 'm<-', 'linewidth',2);
plot(0:20, outage(10, :),'m<--', 'linewidth',2);

% plot(0:20, outage(11, :), 'r*-', 'linewidth',2);
% plot(0:20, outage(12, :), 'r*--', 'linewidth',2);

plot(0:20, outage_uniform,'ks-', 'linewidth',2);




    

