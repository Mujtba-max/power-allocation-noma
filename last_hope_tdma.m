% Calculating CPU time with number of users
Pmax = 10^(P/10); % Maximum power for each cell.
Pm = Pmax*ones(1,NC);

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


for alpha_idx = alpha_rng
    if mod(alpha_idx,2) == 1
        for s = 1:num_reals
            HH = abs(h(:,:,:,s));
            dd = d(:,:,:,s);

            % Next, define weights
            alpha = alpha_computation(HH, dd, NU, NC, alpha_idx);
            tdma_rate = alphaTdmaRates(NC,NU,HH,1/(NU)*ones(NU,NC),Pmax,nvar);
            tdma_sum = sum(tdma_rate, 'all');
            tdma_rates(s) = tdma_sum;
            file_name = sprintf('OMA_for_powers/OMA_%dx%dpower%dalpha%d.mat', NC, NU, P, alpha_idx);
            save(file_name, 'tdma_rates');
        end
    end
end
    