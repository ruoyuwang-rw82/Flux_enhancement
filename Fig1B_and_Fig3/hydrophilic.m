clear all

alpha = 0:0.01:0.15;

Jw_dcmd_F = zeros(1, length(alpha));
Jw_dcmd_D = zeros(1, length(alpha));
Jw_vmd = zeros(1, length(alpha));

for i=1:length(alpha)

    [Jw, Jq, Ts, Pv, B] = model0_janus(200e-6, 0.7, 0.2e-6/2, 0.2, 70, 20, 1e3, alpha(i), 1e-6, 'F'); % LMH
    Jw_dcmd_F(i) = Jw;

    [Jw2, Jq2, Ts2, Pv2, B2] = model0_janus(200e-6, 0.7, 0.2e-6/2, 0.2, 70, 20, 1e3, alpha(i), 1e-6, 'D'); % LMH
    Jw_dcmd_D(i) = Jw2;

    [Jw3, Jq3, Ts3, Pv3, B3] = modelvmd0_janus(200e-6, 0.7, 0.2e-6/2, 70, 10e3, 1e3, alpha(i), 1e-6); % LMH
    Jw_vmd(i) = Jw3;



end

figure

plot(alpha*100, Jw_dcmd_F,'-','Color', [253,185,18]/255,'LineWidth',6)
hold on
plot(alpha*100, Jw_dcmd_D,'--','Color', [253,185,18]/255,'LineWidth',6)
hold on
plot(alpha*100, Jw_vmd,'-','Color', [46,117,182]/255,'LineWidth',6)
hold on

xlim([0,15])
ylim([-25,100])
pbaspect([1 1 1])
set(gca, 'FontSize',30,'linewidth', 3)
set(gca, 'YAxisLocation', 'left', 'TickDir', 'out');



%%%%%%%%%%%%%%%%%%%%%

% pressure Data
vmd_0 = [17509, 17475, 10e3];    
vmd_15 = [38827, 38799, 10e3];
dcmd_0 = [16144, 6520, 6512];    
dcmd_5 = [17742, 12882, 12878];




% figure 3D
dH = 40650;

Pvs = @(x) exp(23.5377-4016.3622./(x-38.6339));
T = 273:1:353;
Pv = Pvs(T);

Tr = 373;
P_pre = Pvs(Tr)*exp(-dH/8.314*(1./T-1/Tr));
P_adj = 1.03e5*exp(dH/8.314/373)*exp(-dH*(1-0.08)/8.314./T);
P_adj2 = 1.03e5*exp(dH/8.314/373)*exp(-dH*(1-0.15)/8.314./T);

figure
plot(T,P_pre/1e3,'-','Color', [46,117,182]/255,'LineWidth',6)
hold on
plot(T,P_adj2/1e3,'-','Color', [46,117,182,150]/255,'LineWidth',6)
hold on
%set(gca, 'YScale', 'log')
pbaspect([1 1 1])
set(gca, 'FontSize',30,'linewidth', 3)
set(gca, 'YAxisLocation', 'left', 'TickDir', 'out');
plot(296, 38.827,'.','MarkerSize',5)
hold on
plot(328, 17.509,'.','MarkerSize',5)
ylim([0,60])
xlim([280,340])