% Calculating Outage Probabilities 
NC = 3; % #Cells which equals #BSs
NU = 4; % #USers in each cell.
nvar = 1.9905e-08; % Noise Variance ?????
epsilon = 1e-5; % For convergence test.
inner_radius = 500; 
minR_ratio = 0.01;
numIter = 2000;
num_reals = 1000;
alpha_rng = 1:10;
P_max_idx = 20;
P = 20;
seed = 1;

alpha_idx = 5;

outage = 0;

clear H  in D;
fileName = sprintf('channels_for_powers/Channels%dx%d power%d.mat', NC, NU, P);
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

last_hope;

for s=1:num_reals
    HH = abs(h(:,:,:,s));
    dd = d(:,:,:,s);

    R = rate(NC, NU, HH, Powers(:,:,s), nvar);
    outage = outage + sum(R<0.005, 'all');
end
outage = outage/(NC*NU*num_reals);
