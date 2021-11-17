
function [] = plot_conv(alpha_rng, conv, WR_vs_iter)
  "plotting convergence"
  
  NIter = zeros(1, 10);
  S = zeros(1, 10);
  for alpha_idx = alpha_rng
    % find the sample corresponds to the average number of iterations for each alpha:
    NIter(alpha_idx) = max(conv(alpha_idx, :));
    [~,S(alpha_idx)]=min(abs(conv(alpha_idx, :) - NIter(alpha_idx))); %S is the sample which takes the average number of iterations for the algorithm to converge.
  end
  maxNIter = max(NIter);
  markers = 'm<-bo-';
  figure; hold on; grid on;
    plot(0:maxNIter/2, WR_vs_iter(S(1), 1:maxNIter/2+1, 1), 'm<-', 'linewidth', 2)
    plot(0:maxNIter/2, WR_vs_iter(S(10), 1:maxNIter/2+1, 10), 'ro-','linewidth', 2)
%     plot(0:maxNIter/2, WR_vs_iter(S(3), 1:maxNIter/2+1, 3), 'b<-', 'linewidth', 1)
%     plot(0:maxNIter/2, WR_vs_iter(S(4), 1:maxNIter/2+1, 4), 'ko-','linewidth', 1)
%     plot(0:maxNIter/2, WR_vs_iter(S(5), 1:maxNIter/2+1, 5), 'g<-', 'linewidth', 1)
%     plot(0:maxNIter/2, WR_vs_iter(S(6), 1:maxNIter/2+1, 6), 'co-','linewidth', 1)
%     plot(0:maxNIter/2, WR_vs_iter(S(7), 1:maxNIter/2+1, 7), 'm<--', 'linewidth', 1)
%     plot(0:maxNIter/2, WR_vs_iter(S(8), 1:maxNIter/2+1, 8), 'ro--','linewidth', 1)
%     plot(0:maxNIter/2, WR_vs_iter(S(9), 1:maxNIter/2+1, 9), 'm<--', 'linewidth', 1)
%     plot(0:maxNIter/2, WR_vs_iter(S(10), 1:maxNIter/2+1, 10), 'ko--','linewidth', 1)
    



plot(0:20, RR_max(1, :), 'ks-.', 'linewidth',2);
plot(0:20, tdma(1, :), 'gd--', 'linewidth',2);
  
  xlabel("number of iterations", 'FontSize', 15);
  ylabel("Weighted sum rate", 'FontSize', 15);
  legend('distance-based alpha', 'uniform alpha', 'FontSize', 15);
end