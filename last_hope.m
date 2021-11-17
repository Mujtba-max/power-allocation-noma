% clear; close all;
% 
% NC = 5; % #Cells which equals #BSs
% NU = 3; % #USers in each cell.
% P = 16;
% nvar = 1.9905e-08; % Noise Variance ?????
% epsilon = 1e-5; % For convergence test.
% inner_radius = 500; 
% minR_ratio = 0.01;
% numIter = 2000;
% num_reals = 1000;
% alpha_rng = 2;
% seed = 1;

Pmax = 10^(P/10); % Maximum power for each cell.
Pm = Pmax*ones(1,NC);

%system model
% [h, ms, Cell] = generate_IBC_channel(NU, inner_radius, NC, minR_ratio, seed, 0);
% fileName = sprintf('Channels%dx%d.mat', NC, NU);
% load(fileName,'H', 'in', 'D'); 

% n = zeros(1,NC);
% for c=1:NC
%     a = find(in(:,c)==1);
%     n(c)= length(a);
% end
% 
% num_reals = min(n);

% 
% for c=1:NC
%     a = find(in(:,c)==1);
%     h(:,c,:,:) = H(:,c,:,a(1:num_reals));
%     d(:,c,:,:) = D(:,c,:,a(1:num_reals));
% end


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
R_vs_iter = zeros(num_reals, numIter, 10);
WR_vs_iter= zeros(num_reals, numIter, 10);
ppp = zeros(num_reals, numIter);

conv= zeros(10,num_reals);
 
for alpha_idx = alpha_rng
Powers = zeros(NU,NC,num_reals);
for s = 1:num_reals
    
%     disp(s)
    if mod(s, 500) == 0
        to_disp = sprintf('%d. %d', P, s);
        disp(to_disp)
    end
%     alpha = ones(NU,NC);
    HH = abs(h(:,:,:,s));
    dd = d(:,:,:,s);
    
    % Next, define weights
    alpha = alpha_computation(HH, dd, NU, NC, alpha_idx);
%     alpha = ones(NU, NC);
%     for c=1:NC
%         alpha_sum = 0;
%         for u=1:NU
% %             alpha(u,c) = u^2;
% %             alpha(u,c) = (HH(u,c,c)/(sum(HH(u,c,:))-HH(u,c,c)))^-2;
%             alpha(u,c) = (dd(u,c,c)/((sum(dd(u,c,:))-dd(u,c,c))))^2;
%             alpha_sum = alpha_sum+alpha(u,c);
%         end
%         alpha(:,c) = alpha(:,c)/alpha_sum;
%     end

    g = ones(NC,NU)'; % Receivers' gains.
    w = ones(NC,NU)'; % Weights in WMMSE.
    A = zeros(NU, NC);
    lambda = zeros(1, NC);

    % [P1max, P2max] = towUsersMaxPower(Pm, Ptol, h, NC);
%     v_init = zeros(NU,NC); vs = zeros(NU,NC);
%     v_init = sqrt(repmat(Pm,NU,1)./NU);
    v_init = sqrt(Pmax*rand(NU,NC));
    vs = v_init;
    % Calculating the corresponding values of g & w for all users.
    for c=1:NC
        for u=1:NU
            intra = HH(u, c, c)^2*sum(v_init(1:u-1, c).^2);
            inter = 0;
            for k=1:NC
                if k~=c
                    temp = HH(u, c, k)^2*sum(v_init(:, k).^2);
                    inter = inter + temp;
                end
            end
            g(u, c) = HH(u, c, c)*v_init(u, c)/(HH(u, c, c)^2 * v_init(u, c)^2 +inter+intra+nvar);
            w(u, c) = 1/((g(u, c)*HH(u, c, c)*v_init(u, c)-1)^2+g(u, c).^2*(inter+intra+nvar));
            A(u, c) = alpha(u, c) * w(u, c) * g(u, c);
        end
    end


    vnew = 0; vold = 0; iter=0;

    combination = zeros(1,2);
    

    while(iter <= numIter)
        iter = iter+1;
        R_vs_iter(s, iter, alpha_idx) = sum(rate(NC, NU, HH, vs, nvar), 'all');
        WR_vs_iter(s, iter, alpha_idx)= sum(Wrate(NC, NU, HH, vs, alpha, nvar), 'all');
        ppp(s, iter) = vs(1, 1);

        for c = 1:NC
            inter = 0;
            for k=1:NC
                if k~=c
                    temp = sum(g(:, k).*A(:, k) .*HH(:, k, c).^2);
                    inter = inter + temp;
                end
            end

            for j=1 

                for u=1:NU
                    vs(u,c) = A(u, c) * HH(u, c, c)/(sum(g(u+1:NU, c).*A(u+1:NU, c).*HH(u+1:NU,c,c).^2) + inter);
                end
                compPower = sum(vs(:,c).^2) - Pm(c);
                if compPower <= 0
                    combination(1) = combination(1) +1;
                    break;
                end
                
                %%%%%%%%%%%%%%%%

                l = bisection(Pmax, A, HH, g, inter, c, NU); % Finding lambda using bisection search.
                if numel(l(l>=0)) > 1
                    error("numel(lambda(c))> 1");
                end
                if numel(l(l>=0)) == 0
                    "there is no lambda that satisfies the KKT conditions for this iteration and this cell." %%???
                    error("numel(lambda(c))== 0");
                end
                lambda(c) = l(l>=0);

                for u=1:NU
                    vs(u,c) = A(u, c) * HH(u, c, c)/(sum(g(u+1:NU, c).*A(u+1:NU, c).*HH(u+1:NU,c,c).^2)+ inter + lambda(c));
                end

                compPower = sum(vs(:,c).^2) - Pm(c) ;       % <= 0

                if abs(compPower) > 1e-7
                    error("the solution of the quatric equation is wrong!");
                end

                combination(2) = combination(2) +1;
                break;


            end

        end

        vold = vnew;
        vnew = sum(log2(w),'all');
        if vnew-vold < epsilon
           conv(alpha_idx, s) = iter;
           break;
        end


        for c=1:NC
            for u=1:NU
                intra = HH(u, c, c)^2*sum(vs(1:u-1, c).^2);
                inter = 0;
                for k=1:NC
                    if k~=c
                        temp = HH(u, c, k)^2*sum(vs(:, k).^2);
                        inter = inter + temp;
                    end
                end
                g(u, c) = HH(u, c, c) * vs(u, c)/(HH(u, c, c)^2 * vs(u, c)^2 +inter+intra+nvar);
                w(u, c) = 1/((g(u, c)*HH(u, c, c)*vs(u, c)-1)^2+g(u, c).^2*(inter+intra+nvar));
                A(u, c) = alpha(u, c) * w(u, c) * g(u, c);
            end
        end
    end
    
     R_vs_iter(s, conv(alpha_idx, s)+1:end, alpha_idx) =  R_vs_iter(s, conv(alpha_idx, s), alpha_idx);
    WR_vs_iter(s, conv(alpha_idx, s)+1:end, alpha_idx) = WR_vs_iter(s, conv(alpha_idx, s), alpha_idx);

%     vs = vs'; g = g'; w = w'; alpha = alpha'; v_init = v_init';
    R = rate(NC, NU, HH, vs, nvar);
    Rmax = rate(NC, NU, HH, v_init, nvar);
    
    WR = Wrate(NC, NU, HH, vs, alpha, nvar);
    WRmax = Wrate(NC, NU, HH, v_init, alpha, nvar);
    
    tdma_rate = tdmaRates(NC,NU,HH,alpha,Pmax,nvar);

%     MSE    = MMSE(NC, NU, HH, vs, alpha, g, w, nvar);
%     MSE_max = MMSE(NC, NU, HH, v_init, alpha, g, w, nvar);

    R_sum    = sum(R   , 'all');
    Rmax_sum = sum(Rmax, 'all');
    
    WR_sum    = sum(WR   , 'all');
    WRmax_sum = sum(WRmax, 'all');
    
    tdma_sum = sum(tdma_rate, 'all');
    
    R_sums(s) = R_sum;
    Rmax_sums(s) = Rmax_sum;
    
    WR_sums(s) = WR_sum;
    WRmax_sums(s) = WRmax_sum;
    
    tdma_rates(s) = tdma_sum;
    
    Powers(:,:,s) = vs;
    
%     vs = vs'; g = g'; w = w'; alpha = alpha'; v_init = v_init';
    
end

% file_name = sprintf('WMMSE_for_powers_fixed/WMMSE_%dx%dpower%dalpha%d.mat', NC, NU, P, alpha_idx);
% save(file_name, 'Powers', 'conv', 'R_sums', 'Rmax_sums', 'WR_sums', 'WRmax_sums', 'tdma_rates');


end