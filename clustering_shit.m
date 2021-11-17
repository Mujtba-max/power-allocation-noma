% NC = 4; NU = 100; NS = 2;
% inner_radius = 500; minR_ratio = 0.04; seed = 1;
% Pm = 10^1.6; nv = 1.9905e-08;
function [h, d, users_clusters, n] = clustering_shit(NC, inner_radius, minR_ratio, seed, Pm, nv)
NS = 1; NU = 100;
A = zeros(NC,NU,2,NS);
for c = 1:NC
    for u = 1:NU
        A(c,u,1,:) = min(inner_radius,inner_radius*rand(1,NS) + minR_ratio*inner_radius);
        A(c,u,2,:) = 2*pi*rand(1,NS);
    end
end
A = permute(A, [4,2,3,1]);
H = zeros(NU,NC,NC,NS);

[h, d, ~, ~] = generate_IBC_channel_fixed(NU, inner_radius, NC, minR_ratio, seed, squeeze(A(1,:,1,:)), squeeze(A(1,:,2,:)), 0);
h = permute(h, [3,2,1]);
h_cell = zeros(NC,NU);
for c=1:NC
    h_cell(c,:) = h(c,c,:);
end
Pmax = Pm*ones(1,NC);
U = NU*ones(1,NC);
nvar = nv*ones(NC,NU);

users_clusters = zeros(NC,NU);
NU_per_cluster = 6;
for c=1:NC
    current_cluster = 1;
%     for u=1:NU
%         if(users_clusters(c,u)~=0)
%             break;
%         end
%         users_clusters(c,u)=current_cluster;
%         for i=u+1:NU
%             if (length(user_clusters(c,:))==NU_per_cluster)
%                 break;
%             end
%             if(suff(u,i) && users_clusters(c,i) ~= 0)
%                 users_clusters(c,i) = current_cluster;
%             end
%         end
%         current_cluster = current_cluster+1;
%     end
    ava_users = find(users_clusters(c,:)==0);
    while(true)
        current_user = ava_users(1);
        users_clusters(c,current_user) = current_cluster;
        if (length(ava_users) == 1)
            break;
        end
        i = 1; %Counter of users assigned for the current_cluster
        u = 2;
        while (u ~= length(ava_users) && i < NU_per_cluster)
            next_user = ava_users(u);
            if suff(current_user,next_user,Pmax,h,NC,c,nvar)
                current_user = next_user;
                users_clusters(c,next_user) = current_cluster;
                i = i+1;
            end
            u = u+1;
        end
        current_cluster = current_cluster+1;
        ava_users = find(users_clusters(c,:)==0);
    end
      
%     for u=current_user:length(ava_users)
%         if suff(current_user, u)
%             current_user = u;
%             users_clusters = current_cluster;
%         end
%     end
end




for k=1:NC
    m(k) = max(users_clusters(k,:));
end
ma = max(m);
n = zeros(NC,ma);
for j=1:NC
    for v=1:ma
        n(j,v) = length(find(users_clusters(j,:)==v));
    end
end

% for l=1:NC
%     disp(length(find(n(l,:)==1)));
% end

% desired_NU = 7;
% a = zeros(NC,desired_NU); 
% for c=1:NC
% [num, clus] = max(n(c,:));
% b = find(users_clusters==clus);
% a(c,:) = b(1:desired_NU);
% end
% HHH = zeros(desired_NU,NC,NC);
% 
% for c=1:NC
% HHH(:,c,:) = permute(h(:,c,a(c,:)), [3,2,1]);
% end
