for P=1:20
    clear H  in D;
    last_hope
end
NC = 3; NU = 2;
clear RR RR_max WRR WRR_max convv;
for P = 1:21
    file_name = sprintf('WMMSE_%dx%dpower%ddis.mat', NC, NU,P-1);

    load(file_name, 'Powers', 'conv', 'R_sums', 'Rmax_sums', 'WR_sums', 'WRmax_sums', 'tdma_rates');
    RR(P) = mean(R_sums);
    RR_max(P) = mean(Rmax_sums);
    WRR(P) = mean(WR_sums);
    WRR_max(P) = mean(WRmax_sums);
    tdma(P) = mean(tdma_rates);
    % PP(:, :, P) = mean(Powers, 3);
    convv(P) = mean(conv); 
end

plot(0:20, RR);
hold on;
plot(0:20, RR_max);
legend('WMMSE Rate', 'Equal Powers Rate');
figure;
plot(0:20, convv);