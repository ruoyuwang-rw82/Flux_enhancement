clear all

alpha = 0:0.01:0.4;

Jw_dcmd_F = zeros(1, length(alpha));
Jw_dcmd_D = zeros(1, length(alpha));
Jw_vmd = zeros(1, length(alpha));

Jq_dcmd_F = zeros(1, length(alpha));
Jq_dcmd_D = zeros(1, length(alpha));
Jq_vmd = zeros(1, length(alpha));

T_dcmd_F = zeros(2, length(alpha));
T_dcmd_D = zeros(2, length(alpha));
T_vmd = zeros(1, length(alpha));

for i=1:length(alpha)

    [Jw, Jq, Ts, Pv, B] = model0_janus(200e-6, 0.7, 0.2e-6/2, 0.2, 70, 20, 1e3, alpha(i), 1e-6, 'F'); % LMH
    Jw_dcmd_F(i) = Jw;
    Jq_dcmd_F(i) = Jq;
    T_dcmd_F(:,i) = Ts;

    [Jw2, Jq2, Ts2, Pv2, B2] = model0_janus(200e-6, 0.7, 0.2e-6/2, 0.2, 70, 20, 1e3, alpha(i), 1e-6, 'D'); % LMH
    Jw_dcmd_D(i) = Jw2;
    Jq_dcmd_D(i) = Jq2;
    T_dcmd_D(:,i) = Ts2;

    [Jw3, Jq3, Ts3, Pv3, B3] = modelvmd0_janus(200e-6, 0.7, 0.2e-6/2, 70, 11e3, 1e3, alpha(i), 1e-6); % LMH
    Jw_vmd(i) = Jw3;
    Jq_vmd(i) = Jq3;
    T_vmd(i) = Ts3;



end

figure

plot(alpha*100, Jw_dcmd_F,'-','Color', [253,185,18]/255,'LineWidth',6)
hold on
plot(alpha*100, Jw_dcmd_D,'--','Color', [253,185,18]/255,'LineWidth',6)
hold on
plot(alpha*100, Jw_vmd,'-','Color', [46,117,182]/255,'LineWidth',6)
hold on
xlim([0,40])
ylim([0,30])
pbaspect([1.5 1 1])
set(gca, 'FontSize',30,'linewidth', 3)
set(gca, 'YAxisLocation', 'left', 'TickDir', 'out');


% figure
% 
% plot(alpha*100, Jq_dcmd_F,'-','Color', [253,185,18]/255,'LineWidth',6)
% hold on
% plot(alpha*100, Jq_dcmd_D,'--','Color', [253,185,18]/255,'LineWidth',6)
% hold on
% plot(alpha*100, Jq_vmd,'-','Color', [46,117,182]/255,'LineWidth',6)
% hold on
% plot(alpha*100, ones(1,length(alpha))*Jq_vmd(1), 'k-.','LineWidth',3)
% hold on
% plot(alpha*100, ones(1,length(alpha))*Jq_dcmd_F(1), 'k-.','LineWidth',3)
% hold on
% 
% xlim([0,50])
% %ylim([10,60])
% pbaspect([1.5 1 1])
% set(gca, 'FontSize',30,'linewidth', 3)
% set(gca, 'YAxisLocation', 'left', 'TickDir', 'out');
% 
% figure
% plot(alpha*100,T_dcmd_F)
% 
% figure
% plot(alpha*100,T_dcmd_D)
% 
% figure
% plot(alpha*100,T_vmd)