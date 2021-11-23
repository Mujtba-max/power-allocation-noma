NC = 4; NU = 3; NS = 10^4;
inner_radius = 500; minR_ratio = 0.01; seed = 1;
% Pm = 10^1.6; 
Pmax = 16; % dBw
nvar_dBm = 10^-17.4; %1.9905e-08
BW = 5 * 10^6;
nvar = nvar_dBm * 10^-3 * BW; % linear scale
rng(seed);
    Pm = 10^(Pmax/10);
    tic;

for NC = 4
for NU = 10
%     Pm = 10^(Pmax/10);
    disp(Pmax);
    A = zeros(NC,NU,2,NS);
    R = inner_radius - minR_ratio*inner_radius;      % effective cell radius
    for c = 1:NC
        for u = 1:NU
            d = sum(rand(2,NS),1) .* R;              % user distributed in the cell uniformly
            d(d>R) = 2*R - d(d>R);
            A(c,u,1,:) = d + minR_ratio*inner_radius;               % real MS location
            A(c,u,2,:) = 2*pi*rand(1,NS);
        end               
    end
    
%     for c = 1:NC
%         for u = 1:NU
%             A(c,u,1,:) = min(inner_radius,inner_radius*rand(1,NS) + minR_ratio*inner_radius);
%             A(c,u,2,:) = 2*pi*rand(1,NS);
%         end
%     end
    A = permute(A, [4,2,3,1]);
    H = zeros(NU,NC,NC,NS);
    D = zeros(NU,NC,NC,NS);
    in = ones(NS, NC);

    for s=1:NS
        [h, distances, ms, Cell] = generate_IBC_channel_fixed(NU, inner_radius, NC, minR_ratio, seed, squeeze(A(s,:,1,:)), squeeze(A(s,:,2,:)), 0);
        H(:,:,:,s) = h;
        D(:,:,:,s) = distances;
        p = test(NC, NU, inner_radius, minR_ratio,h, Pm, nvar, seed);
        for c = 1:NC 
            if (p(c)~=0)
                in(s,c) = 0;
                H(:, c, :, s) = NaN;
            end
        end
    end
     % if there is a cell that doesn't satisfy the sufficient condition:
%      for c = 1:NC
%          if sum(in(:, c), 1) == 0
%             [h_c, d_c, users_clusters, n] = clustering_shit(NC, inner_radius, minR_ratio, seed, Pm, nvar);
%             h_c = permute(h_c, [3,2,1]);
%             i = 1;
%             for clus_indx = find(n(c, :) == NU)
%                 i_c = find(users_clusters(c,:) == clus_indx);
%                 h_res = h_c(i_c, c, :);
%                 d_res = d_c(i_c, c, :);
%                 H(:, c, :, i) = h_res;
%                 in(i, c) = 1;
%                 %check if it realy satisfy it or not:
%                 p = test(NC, NU, inner_radius, minR_ratio,H(:, :, :, i), Pm, nvar, seed);
%                 if (p(c)~=0)
%                     in(i,c) = 0;
%                     H(:, c, :, i) = NaN;
%                 end
%                 i = i + 1;                
%             end
%          end
%      end
    
%     file_name = sprintf('channels_for_NU_fixed/Channels%dx%dpower%d.mat', NC, NU, Pmax);
%     save(file_name, 'H', 'in','D', nvar);
    toc;
end
end







% n = zeros(1,NC);
% for c=1:NC
%     a = find(in2(:,c)==1);
%     n(c)= length(a);
% end
% 
% num_reals = min(n);
% h = zeros(NU,NC,NC,num_reals);
% d = zeros(NU,NC,NC,num_reals);
% 
% for c=1:NC
%     a = find(in2(:,c)==1);
%     h(:,c,:,:) = H2(:,c,:,a(1:num_reals));
%     d(:,c,:,:) = D(:,c,:,a(1:num_reals));
% end
% save('H5x3.mat', 'h', 'd', 'num_reals');
% save('H5x3V2.mat', 'h', 'd', 'num_reals');
%%













