function lambda = bisection(Pm, A, H, g, inter, c, NU)
lambda = 0;
start_point = 0;
end_point = 100;
while true
    x = sum_power(end_point,A, H, g, inter, c, NU);
    if x > Pm
        end_point = end_point*2;
    else
        break;
    end
end
while true
%     pause(1);
    lambda = (end_point+start_point)/2;
    s = sum_power(lambda,A, H, g, inter, c, NU);
    if abs(s-Pm) < 1e-7
        break;
    end
    if s>Pm
        start_point = lambda;
    else
        end_point = lambda;
    end
end
end

function s = sum_power(lambda,A, H, g, inter, c, NU)
vs = zeros(size(H,1),1);
for u=1:NU
    vs(u) = A(u, c) * H(u, c, c)/(sum(g(u+1:NU, c).*A(u+1:NU, c).*H(u+1:NU,c,c).^2) + inter + lambda);
%     intra_users = ones(NU,1);
%     intra_users(u) = 0;
%     vs(u) = A(u, c) * H(u, c, c)/(sum(intra_users.*g(:, c).*A(:, c).*H(:,c,c).^2)+ inter + lambda);
end
s = sum(vs.^2);
end
