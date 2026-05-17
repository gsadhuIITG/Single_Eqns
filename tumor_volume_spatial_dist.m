clc; clear; close all;

%% ================= Parameters =================
Nt = 51;                 % time points t = 0,...,50
Nx = 1001;               % spatial grid
data_dir = '.';          % data directory

times_to_plot = [10 25 35 50];

%% ================= Preallocate =================
p = zeros(Nx, Nt);
m = zeros(Nx, Nt);
s = zeros(Nx, Nt);
x = zeros(Nx, 1);

%% ================= Load data =================
for t = 0:Nt-1
    % p and m
    data_pm = load(sprintf('%s/pm_t_%d.txt', data_dir, t));
    x         = data_pm(:,1);
    p(:,t+1) = data_pm(:,2);
    m(:,t+1) = data_pm(:,3);

    % s
    data_s   = load(sprintf('%s/s_t_%d.txt', data_dir, t));
    s(:,t+1) = data_s(:,2);
end

n = p + m;              % total density
tvec = 0:Nt-1;

%% ================= Total mass =================
dx = x(2) - x(1); %#ok<NASGU>
mass_p = trapz(x, p);
mass_m = trapz(x, m);
mass_s = trapz(x, s);
mass_n = mass_p + mass_m;

%% =========================================================
% Spatial profiles: p, m, p+m, s
% =========================================================
%% =========================================================
% Spatial profiles: subplots by time
% =========================================================
fig1 = figure('Position',[1 1 900 650]);

% --- Fixed colors for variables (consistent across time) ---
col_p = [0.85 0.33 0.10];   % p
col_m = [0.00 0.45 0.74];   % m
col_n = [0.00 0.00 0.00];   % p+m
col_s = [0.47 0.67 0.19];   % s

nT = length(times_to_plot);

for k = 1:nT
    tt = times_to_plot(k);

    subplot(2,2,k); hold on;

    plot(x, p(:,tt+1), '-',  'LineWidth',2,   'Color',col_p);
    plot(x, m(:,tt+1), '--', 'LineWidth',2,   'Color',col_m);
    plot(x, n(:,tt+1), ':',  'LineWidth',2.5, 'Color',col_n);
    plot(x, s(:,tt+1), '-.', 'LineWidth',1.8, 'Color',col_s);

    title(sprintf('t = %d', tt));
    xlabel('x');
    ylabel('Density');
    box on;
end

% --- Single shared legend ---
legend({'p','m','p+m','n'}, ...
       'Orientation','horizontal', ...
       'Location','southoutside');

set(findall(fig1,'-property','FontSize'),'FontSize',13);

%% =========================================================
% Temporal evolution of total mass
% Normoxic (light gray) vs Hypoxic (dark gray)
% =========================================================
fig2 = figure('Position',[1 1 620 518]); hold on;

% --- Plot curves ---
h_p = plot(tvec, mass_p, 'r-',  'LineWidth', 2);
h_m = plot(tvec, mass_m, 'b-',  'LineWidth', 2);
h_n = plot(tvec, mass_n, 'k-',  'LineWidth', 2);
h_s = plot(tvec, mass_s, 'g--', 'LineWidth', 2.5);

xlabel('t');
ylabel('Total mass');
set(gca,'FontSize',13);
box on;

% --- Axis limits ---
xlim([tvec(1) tvec(end)]);
yl = ylim;

% --- Stronger contrasting bands ---
band_width = 5;

normoxic_color = [0.95 0.95 0.95];  % very light gray
hypoxic_color  = [0.70 0.70 0.70];  % clearly darker gray

for t0 = tvec(1):band_width:tvec(end)
    if mod(floor(t0/band_width),2) == 0
        band_color = normoxic_color;   % Normoxic
    else
        band_color = hypoxic_color;    % Hypoxic
    end

    patch([t0 t0+band_width t0+band_width t0], ...
          [yl(1) yl(1) yl(2) yl(2)], ...
          band_color, ...
          'EdgeColor','none', ...
          'FaceAlpha',0.6);   % higher opacity
end

% --- Phase boundary markers (critical for clarity) ---
for t0 = tvec(1):band_width:tvec(end)
    xline(t0,'k:','LineWidth',0.8);
end

% --- Bring curves to front ---
uistack(findall(gca,'Type','line'),'top');

% --- Legend patches ---
h_norm = patch(NaN,NaN,normoxic_color,'FaceAlpha',0.6,'EdgeColor','none');
h_hypo = patch(NaN,NaN,hypoxic_color ,'FaceAlpha',0.6,'EdgeColor','none');

legend([h_p h_m h_n h_s h_norm h_hypo], ...
       {'p','m','p+m','n','Normoxic','Hypoxic'}, ...
       'Location','best');
box on;
