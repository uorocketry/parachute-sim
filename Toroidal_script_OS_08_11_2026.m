% Toroidal parachute: ONE parachute, TWO configurations (reefed / fully open)
% All quantities are in SI units (m, kg, s, N, kg/m^3).
%
% - Sizing from target speeds, with ISA air density that varies with altitude
% - 2D descent simulation (same model as parachute_sim.py) with labelled time plots
% - Full-scale gore template written directly to PDF (like paracalc.py)
%
% FOR A NEW ROCKET / NEW PARACHUTE: edit only section 1 (INPUTS).
% Lines marked  <<< EDIT  must be reviewed; lines marked  <<< EDIT (optional)  have defaults.


clear
close all
clc

%% =====================================================================
%  1. INPUTS  -  EDIT ONLY THIS SECTION FOR A NEW ROCKET / NEW PARACHUTE
%  =====================================================================

% ---------- Rocket and flight (CHANGE FOR EACH NEW ROCKET / FLIGHT) ----------
m               = 1.5;    % <<< EDIT  mass under parachute, dry mass (kg)
site_alt_asl    = 100;    % <<< EDIT  launch/landing site altitude (m ASL)
apogee_agl      = 432;   % <<< EDIT  expected apogee = reefed deployment altitude (m AGL)
open_deploy_agl = 200;    % <<< EDIT  reefing cut / fully open deployment altitude (m AGL)
vx_init         = 50.5;   % <<< EDIT  horizontal velocity at apogee (m/s)

% ---------- Descent speed targets and drag coefficients ----------
v_target_reefed = 20;     % <<< EDIT  reefed descent speed target (m/s)
v_target_open   = 7;      % <<< EDIT  fully open descent speed target (m/s, at landing)
Cd_reefed       = 1.4;    % <<< EDIT  PLACEHOLDER (not sourced): replace with measured value
Cd_open         = 2.2;    % <<< EDIT  Fruity Chutes specs for a comparable parachute + other research
reef_speed_ref_agl = open_deploy_agl;   % <<< EDIT (optional) altitude where v_target_reefed applies
                                        %     = open_deploy_agl (speed at reefing cut) or apogee_agl

% ---------- Number of gores ----------
num_gores = 10;           % <<< EDIT  NUMBER OF GORES

% ---------- Canopy proportions (design choice) ----------
b_over_A = 0.5;           % <<< EDIT (optional) max height / base radius (original design: 22.5/45)
C_over_A = 6.25/45;       % <<< EDIT (optional) spill hole radius / base radius (original design)

% ---------- Weather and test data ----------
dT_isa      = 0;          % <<< EDIT (optional) temperature offset from standard day (K), e.g. +15 on a hot day
vw_10       = 0;          % <<< EDIT (optional) wind velocity at 10 m (m/s)
vw_alpha    = 0.27;       % <<< EDIT (optional) Hellmann exponent
t_open_reef = 0.0;        % <<< EDIT (optional) opening time, reefed parachute at apogee (s), 0 = instant
t_open_full = 0.0;        % <<< EDIT (optional) opening time, reefed -> fully open (s), 0 = instant

% ---------- Gore template (PDF) ----------
seam          = 0.0254;   % <<< EDIT (optional) seam allowance: 1 inch (m)
pdf_margin_cm = 2;        % <<< EDIT (optional) page margin around the cut line (cm)
pdf_labels    = true;     % <<< EDIT (optional) title + 10 cm scale-check bar in the margins
pdf_centre    = true;     % <<< EDIT (optional) thin centre line (fold / grain reference)
pdf_folder    = pwd;      % <<< EDIT (optional) folder where the PDF is saved
%   The file name is built automatically: toroidal_gore_template_<open diameter in m>_<number of gores>.pdf

%% =====================================================================
%  END OF INPUTS - nothing to edit below this line
%  =====================================================================

%% ===== 2. Constants and atmosphere model =====
g = 9.81;                % m/s^2

% ISA density (troposphere, valid up to 11000 m ASL); z in m ASL
rho_fun  = @(z) 101325*(1 - 2.25577e-5*z).^5.25588 ./ (287.05*(288.15 - 0.0065*z + dT_isa));
wind_fun = @(h) vw_10*(max(h, 0)/10).^vw_alpha;     % h in m AGL

%% ===== 3. Configurations =====
% Cd is defined with respect to the EFFECTIVE projected area of each
% configuration: S_eff = pi*R_eff^2 (reefed skirt for reefed, full base for open).
% rho for sizing is evaluated at the altitude where the target speed applies.
cfg(1).name         = 'Reefed';
cfg(1).v_target     = v_target_reefed;       % m/s
cfg(1).z_design_agl = reef_speed_ref_agl;    % m AGL where v_target applies
cfg(1).Cd           = Cd_reefed;

cfg(2).name         = 'Fully open';
cfg(2).v_target     = v_target_open;         % m/s
cfg(2).z_design_agl = 0;                     % landing speed, at site altitude
cfg(2).Cd           = Cd_open;

%% ===== 4. Sizing from target speeds =====
for k = 1:numel(cfg)
    cfg(k).rho   = rho_fun(site_alt_asl + cfg(k).z_design_agl);   % kg/m^3
    cfg(k).CdS   = 2*m*g/(cfg(k).rho*cfg(k).v_target^2);          % required drag area (m^2)
    cfg(k).S_eff = cfg(k).CdS/cfg(k).Cd;                          % effective projected area (m^2)
    cfg(k).R_eff = sqrt(cfg(k).S_eff/pi);                         % effective radius (m)
end

%% ===== 5. Canopy geometry (m) - same parachute for both configurations =====
A = cfg(2).R_eff;        % base radius = fully open effective radius
b = b_over_A*A;          % maximum canopy height (m)
C = C_over_A*A;          % spill hole radius (m)

% Vertical distance between the skirt plane (radius A) and the spill hole edge
% (radius C). Assumes a circular cross-section (A = 2*b).
h = sqrt(b^2 - (b - C)^2);

D_open     = 2*A;                % fully open diameter (m)
D_spill    = 2*C;                % spill hole diameter (m)
D_reefed   = 2*cfg(1).R_eff;     % reefed skirt diameter (m)
reef_ratio = cfg(1).R_eff/A;     % D_reefed / D_open
L_reef     = 2*pi*cfg(1).R_eff;  % reefing line perimeter, no allowance (m)

%% ===== 6. 2D descent simulation (same model as parachute_sim.py) =====
dt        = 0.001;               % s
max_steps = 600000;

t     = zeros(1, max_steps);
pos_x = zeros(1, max_steps);
alt   = zeros(1, max_steps);     % altitude AGL (m)
vel_x = zeros(1, max_steps);
vel_y = zeros(1, max_steps);
acc_x = zeros(1, max_steps);
acc_y = zeros(1, max_steps);
drag  = zeros(1, max_steps);

alt(1)   = apogee_agl;
vel_x(1) = vx_init;
t_main   = NaN;                  % time when the reefing is cut
n = 1;

while alt(n) > 0 && n < max_steps
    % Openness of the reefed parachute (apogee) and of the full opening
    if t_open_reef > 0
        o_reef = min(1, 0.5*(t(n)/t_open_reef)^2);
    else
        o_reef = 1;
    end
    if alt(n) < open_deploy_agl
        if isnan(t_main)
            t_main = t(n);
        end
        if t_open_full > 0
            o_full = min(1, 0.5*((t(n) - t_main)/t_open_full)^2);
        else
            o_full = 1;
        end
    else
        o_full = 0;
    end
    CdS = o_reef*cfg(1).CdS + o_full*(cfg(2).CdS - cfg(1).CdS);   % drag area (m^2)

    % Relative airflow and drag force
    ux = -vel_x(n) + wind_fun(alt(n));
    uy = -vel_y(n);
    u  = hypot(ux, uy);                                   % airspeed
    fd = 0.5*rho_fun(site_alt_asl + alt(n))*CdS*u^2;
    if u == 0
        fdx = 0;
        fdy = 0;
    else
        fdx = fd*ux/u;
        fdy = fd*uy/u;
    end

    % Accelerations
    acc_x(n) = fdx/m;
    acc_y(n) = fdy/m - g;
    drag(n)  = fd;

    % Time step (semi-implicit Euler)
    vel_x(n+1) = vel_x(n) + acc_x(n)*dt;
    vel_y(n+1) = vel_y(n) + acc_y(n)*dt;
    pos_x(n+1) = pos_x(n) + vel_x(n+1)*dt;
    alt(n+1)   = alt(n)   + vel_y(n+1)*dt;
    t(n+1)     = t(n) + dt;
    n = n + 1;
end

idx   = 1:n-1;                   % samples with computed accelerations / drag
t     = t(idx);      pos_x = pos_x(idx);  alt   = alt(idx);
vel_x = vel_x(idx);  vel_y = vel_y(idx);
acc_x = acc_x(idx);  acc_y = acc_y(idx);  drag  = drag(idx);

i_cut = find(t >= t_main, 1);
v_tot = hypot(vel_x, vel_y);
a_tot = hypot(acc_x, acc_y);
[vz_peak, i_pk] = max(-vel_y(1:i_cut));

%% ===== 7. Results =====
fprintf('\n=== Air density ===\n');
fprintf('At site (%.0f m ASL)            : %.3f kg/m^3\n', site_alt_asl, rho_fun(site_alt_asl));
fprintf('At open deployment (%.0f m AGL)  : %.3f kg/m^3\n', open_deploy_agl, rho_fun(site_alt_asl + open_deploy_agl));
fprintf('At apogee (%.0f m AGL)         : %.3f kg/m^3\n', apogee_agl, rho_fun(site_alt_asl + apogee_agl));

fprintf('\n=== Configurations (sizing) ===\n');
fprintf('%-11s %8s %6s %10s %9s %10s %10s\n', ...
        'Config', 'rho', 'Cd', 'v_tgt(m/s)', 'CdS(m2)', 'S_eff(m2)', 'R_eff(m)');
for k = 1:numel(cfg)
    fprintf('%-11s %8.3f %6.2f %10.1f %9.3f %10.3f %10.3f\n', ...
            cfg(k).name, cfg(k).rho, cfg(k).Cd, cfg(k).v_target, ...
            cfg(k).CdS, cfg(k).S_eff, cfg(k).R_eff);
end

fprintf('\n=== Sensitivity to reefed Cd ===\n');
for Cd_try = [1.2 1.6 2.0 2.2]
    R_try = sqrt(cfg(1).CdS/(Cd_try*pi));
    fprintf('Cd_reefed = %.2f -> reefed radius = %.3f m, reefing ratio = %.3f\n', ...
            Cd_try, R_try, R_try/A);
end

fprintf('\n=== Descent simulation ===\n');
fprintf('Reefed phase duration                  : %.1f s\n', t_main);
fprintf('Peak reefed vertical speed             : %.1f m/s (at %.0f m AGL)\n', vz_peak, alt(i_pk));
fprintf('Vertical speed at reefing cut          : %.1f m/s\n', -vel_y(i_cut));
fprintf('Total speed at reefing cut             : %.1f m/s\n', v_tot(i_cut));
fprintf('Landing vertical speed                 : %.2f m/s\n', -vel_y(end));
fprintf('Peak total drag force                  : %.0f N\n', max(drag));
fprintf('Peak total acceleration                : %.1f g\n', max(a_tot)/g);
fprintf('Total descent time                     : %.1f s\n', t(end));
fprintf('Horizontal drift                       : %.0f m\n', pos_x(end));

%% ===== 8. Gore pattern =====
numpoints = 200;                         % points along the gore edge

angmax = pi - asin(h/b);
ell_a  = (A - C)/(1 + abs(cos(angmax))); % horizontal semi-axis of the profile ellipse (m)
ell_c  = A - 2*ell_a;
if ell_c < -1e-9
    error('Invalid geometry: c < 0. Check A, b, C and h.');
end

theta  = linspace(0, angmax, numpoints);
gore_x = (1/num_gores)*pi*(ell_a*(cos(theta) + 1) + ell_c);   % gore half-width (m)
fun    = @(s) sqrt((ell_a*sin(s)).^2 + (b*cos(s)).^2);
gore_y = arrayfun(@(s) integral(fun, 0, s), theta);           % arc length along the meridian (m)

% Finished gore outline (sewing line) and cut line (outward offset on all edges)
inner = polyshape([gore_x, -fliplr(gore_x)], [gore_y, fliplr(gore_y)]);
outer = polybuffer(inner, seam, 'JointType', 'miter', 'MiterLimit', 2);
Vi = inner.Vertices;
Vo = outer.Vertices;

fprintf('\nGore: base width = %.3f m, top width = %.3f m, length = %.3f m\n', ...
        2*gore_x(1), 2*gore_x(end), gore_y(end));

%% ===== 9. Full-scale PDF template (written directly, 1 cm = 1 cm on paper) =====
pdf_file = fullfile(pdf_folder, sprintf('toroidal_gore_template_%.3fm_%d.pdf', D_open, num_gores));
write_gore_pdf(pdf_file, Vi, Vo, pdf_margin_cm, num_gores, pdf_labels, pdf_centre, gore_y(end));

xspan = (max(Vo(:,1)) - min(Vo(:,1)))*100 + 2*pdf_margin_cm;
yspan = (max(Vo(:,2)) - min(Vo(:,2)))*100 + 2*pdf_margin_cm;
fprintf('PDF template saved: %s\n', pdf_file);
fprintf('Page size: %.1f cm x %.1f cm (%.1f in x %.1f in). Print at 100%%.\n', ...
        xspan, yspan, xspan/2.54, yspan/2.54);

%% ===== 10. 3D view =====
figure
hold on
ph = linspace(0, 2*pi, 100);
for r = linspace(C, A, 40)
    z = b*sqrt(max(0, 1 - ((r - ell_a - ell_c)/ell_a)^2));
    plot3(r*cos(ph), r*sin(ph), z*ones(size(ph)), 'b');
end
plot3(A*cos(ph), A*sin(ph), zeros(size(ph)), 'k', 'LineWidth', 2);
plot3(cfg(1).R_eff*cos(ph), cfg(1).R_eff*sin(ph), zeros(size(ph)), 'r', 'LineWidth', 2);
axis equal
grid on
view(30, 25)
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)')
title('Fully open canopy (black) and reefed skirt (red)')
hold off

%% ===== 11. Descent plots =====
% Labels written above the points: the maximum of each phase (reefed deployment at
% apogee, reefing cut) and the plateaus (constant value for >= 10 s).
figure
plot(-vel_y, alt, 'LineWidth', 1.5)
hold on
yline(open_deploy_agl, '--k');
text(vz_peak, alt(i_pk), sprintf('%.1f m/s', vz_peak), ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 9);
text(-vel_y(i_cut), alt(i_cut), sprintf('%.1f m/s', -vel_y(i_cut)), ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 9);
text(-vel_y(end), 0.5*open_deploy_agl, sprintf('%.1f m/s', -vel_y(end)), ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 9);
xlim([0 1.15*vz_peak])
grid on
xlabel('Vertical speed (m/s)')
ylabel('Altitude AGL (m)')
title('Vertical speed vs altitude')

plots = {
    alt,     'Altitude AGL (m)',              'Altitude',              'm',     '%.1f';
    -vel_y,  'Vertical velocity (m/s)',       'Vertical Velocity',     'm/s',   '%.1f';
    a_tot,   'Total acceleration (m/s^2)',    'Total Acceleration',    'm/s^2', '%.2f';
    acc_y,   'Vertical acceleration (m/s^2)', 'Vertical Acceleration', 'm/s^2', '%.2f';
    drag,    'Total drag force (N)',          'Total Drag Force',      'N',     '%.1f';
    v_tot,   'Total velocity (m/s)',          'Total Velocity',        'm/s',   '%.1f'};
for i = 1:size(plots, 1)
    figure
    plot(t, plots{i, 1}, 'LineWidth', 1.5)
    hold on
    xline(t_main, '--k');
    label_series(t, plots{i, 1}, i_cut, t_main, dt, plots{i, 4}, plots{i, 5});
    grid on
    xlabel('Time (s)')
    ylabel(plots{i, 2})
    title([plots{i, 3} ' vs time'])
end

%% ===== 12. KEY DESIGN VALUES =====
fprintf('\n=== KEY DESIGN VALUES (SI) ===\n');
fprintf('Number of gores (panels)                : %d\n',   num_gores);
fprintf('Fully open diameter (2A)                : %.3f m\n', D_open);
fprintf('Spill hole diameter (2C)                : %.3f m\n', D_spill);
fprintf('Vertical distance pulled inward (h)     : %.3f m\n', h);
fprintf('Maximum canopy height (b)               : %.3f m\n', b);
fprintf('Reefed skirt diameter                   : %.3f m\n', D_reefed);
fprintf('Reefing ratio (D_reefed / D_open)       : %.3f\n',   reef_ratio);
fprintf('Reefing line perimeter (no allowance)   : %.3f m\n', L_reef);

%% =====================================================================
%  LOCAL FUNCTIONS
%  =====================================================================

function label_series(t, data, i_cut, t_main, dt, unit, fmt)
% Writes values directly above the points on the current axes:
%   - the maximum of each phase (reefed deployment at apogee, reefing cut)
%   - the plateaus (constant value for at least 10 s, within 0.5 % of the data range)
% A phase maximum that lies on a plateau is not repeated.
    ax  = gca;
    hold(ax, 'on');
    rng = max(data) - min(data);
    if rng == 0
        rng = 1;
    end

    % Plateaus
    step = max(1, round(0.1/dt));
    P = find_plateaus(t(1:step:end), data(1:step:end), 0.005*rng, 10, t_main);
    for p = 1:size(P, 1)
        text(ax, 0.5*(P(p, 1) + P(p, 2)), P(p, 3), sprintf([fmt ' %s'], P(p, 3), unit), ...
             'FontSize', 9, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end

    % Maximum of each phase
    phases = {1:(i_cut - 1), i_cut:numel(t)};
    for p = 1:2
        ii = phases{p};
        [vm, k] = max(data(ii));
        k = ii(k);
        if any(t(k) >= P(:, 1) & t(k) <= P(:, 2))
            continue
        end
        frac = (t(k) - t(1))/(t(end) - t(1));
        if frac < 0.2
            ha = 'left';
        elseif frac > 0.8
            ha = 'right';
        else
            ha = 'center';
        end
        text(ax, t(k), vm, sprintf([fmt ' %s'], vm, unit), ...
             'FontSize', 9, 'HorizontalAlignment', ha, 'VerticalAlignment', 'bottom');
    end

    ylim(ax, [min(data) - 0.05*rng, max(data) + 0.12*rng]);
end

function P = find_plateaus(t, y, tol, min_dur, t_cut)
% Returns rows [t_start t_end mean_value] for stretches where y stays within a band
% of width tol for at least min_dur seconds. Plateaus never straddle t_cut.
    P = zeros(0, 3);
    ic = find(t >= t_cut, 1);
    if isempty(ic)
        ic = numel(t) + 1;
    end
    bounds = [1, ic - 1; ic, numel(t)];
    for p = 1:2
        i = bounds(p, 1);
        last = bounds(p, 2);
        while i <= last
            lo = y(i);
            hi = y(i);
            j = i;
            while j < last
                lo2 = min(lo, y(j + 1));
                hi2 = max(hi, y(j + 1));
                if hi2 - lo2 > tol
                    break
                end
                lo = lo2;
                hi = hi2;
                j = j + 1;
            end
            if t(j) - t(i) >= min_dur
                P(end + 1, :) = [t(i), t(j), mean(y(i:j))];
                i = j + 1;
            else
                i = i + 1;
            end
        end
    end
    % Merge neighbouring plateaus of the same phase whose values differ by less than tol
    k = 1;
    while k < size(P, 1)
        same_phase = (P(k, 2) < t_cut) == (P(k + 1, 1) < t_cut);
        if same_phase && P(k + 1, 1) - P(k, 2) <= 0.5 && abs(P(k + 1, 3) - P(k, 3)) <= tol
            d1 = P(k, 2) - P(k, 1);
            d2 = P(k + 1, 2) - P(k + 1, 1);
            P(k, 3) = (P(k, 3)*d1 + P(k + 1, 3)*d2)/(d1 + d2);
            P(k, 2) = P(k + 1, 2);
            P(k + 1, :) = [];
        else
            k = k + 1;
        end
    end
end

function write_gore_pdf(pdf_file, Vi, Vo, margin_cm, num_gores, add_labels, add_centre, gore_len)
% Writes a one-page vector PDF at true scale (like paracalc.py + reportlab).
% Vi: finished gore outline (sewing line, dashed) [m]; Vo: cut line (solid) [m].
% The page is sized to the cut line + margins, in PDF points (1 cm = 28.346 pt).
    cm = 72/2.54;                                    % cm -> PDF points
    xmin = min(Vo(:, 1));
    ymin = min(Vo(:, 2));
    Wc = (max(Vo(:, 1)) - xmin)*100 + 2*margin_cm;   % page width (cm)
    Hc = (max(Vo(:, 2)) - ymin)*100 + 2*margin_cm;   % page height (cm)
    W = Wc*cm;
    H = Hc*cm;
    toX = @(X) ((X - xmin)*100 + margin_cm)*cm;
    toY = @(Y) ((Y - ymin)*100 + margin_cm)*cm;

    c = sprintf('0.75 w\n');
    c = [c sprintf('[6 4] 0 d\n') pdf_path(toX(Vi(:, 1)), toY(Vi(:, 2)))];   % sewing line, dashed
    c = [c sprintf('[] 0 d\n')    pdf_path(toX(Vo(:, 1)), toY(Vo(:, 2)))];   % cut line, solid
    if add_centre
        c = [c sprintf('0.5 w\n0.6 G\n[8 3 1 3] 0 d\n')];
        c = [c sprintf('%.3f %.3f m\n%.3f %.3f l\nS\n', toX(0), toY(0), toX(0), toY(gore_len))];
        c = [c sprintf('[] 0 d\n0 G\n0.75 w\n')];
    end
    if add_labels
        c = [c pdf_text(margin_cm*cm, H - 0.8*cm, 9, ...
             sprintf('Gore pattern - cut %d gores - scale 1:1 (print at 100%%, no fit-to-page)', num_gores))];
        c = [c pdf_text(margin_cm*cm, H - 1.4*cm, 9, ...
             'Dashed line = sewing line | Solid line = cut line (1 in = 2.54 cm seam allowance)')];
        c = [c sprintf('1 w\n%.3f %.3f m\n%.3f %.3f l\nS\n', margin_cm*cm, 0.9*cm, (margin_cm + 10)*cm, 0.9*cm)];
        c = [c sprintf('%.3f %.3f m\n%.3f %.3f l\nS\n', margin_cm*cm, 0.6*cm, margin_cm*cm, 1.2*cm)];
        c = [c sprintf('%.3f %.3f m\n%.3f %.3f l\nS\n', (margin_cm + 10)*cm, 0.6*cm, (margin_cm + 10)*cm, 1.2*cm)];
        c = [c pdf_text((margin_cm + 10.5)*cm, 0.75*cm, 9, 'Scale check: this bar = 10 cm')];
    end
    c = c(1:end-1);                                  % drop the last newline (exact /Length)

    objs = cell(1, 5);
    objs{1} = '<< /Type /Catalog /Pages 2 0 R >>';
    objs{2} = '<< /Type /Pages /Kids [3 0 R] /Count 1 >>';
    objs{3} = sprintf(['<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.2f %.2f] ' ...
                       '/Contents 4 0 R /Resources << /Font << /F1 5 0 R >> >> >>'], W, H);
    objs{4} = sprintf('<< /Length %d >>\nstream\n%s\nendstream', numel(c), c);
    objs{5} = '<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>';

    pdf  = sprintf('%%PDF-1.4\n');
    offs = zeros(1, 5);
    for k = 1:5
        offs(k) = numel(pdf);
        pdf = [pdf sprintf('%d 0 obj\n', k) objs{k} sprintf('\nendobj\n')];
    end
    xref_off = numel(pdf);
    xref = sprintf('xref\n0 6\n0000000000 65535 f \n');
    for k = 1:5
        xref = [xref sprintf('%010d 00000 n \n', offs(k))];
    end
    trailer = sprintf('trailer\n<< /Size 6 /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n', xref_off);
    pdf = [pdf xref trailer];

    fid = fopen(pdf_file, 'wb');
    if fid < 0
        error('Cannot write %s (is it open in a PDF viewer?)', pdf_file);
    end
    fwrite(fid, pdf, 'char');
    fclose(fid);
end

function s = pdf_path(px, py)
% Closed, stroked PDF path through the points (px, py) in PDF points.
    px = px(:)';
    py = py(:)';
    s = sprintf('%.3f %.3f m\n', px(1), py(1));
    s = [s sprintf('%.3f %.3f l\n', [px(2:end); py(2:end)]) sprintf('h S\n')];
end

function s = pdf_text(x, y, size, str)
% One line of Helvetica text at (x, y) in PDF points (ASCII only).
    str = strrep(str, '\', '\\');
    str = strrep(str, '(', '\(');
    str = strrep(str, ')', '\)');
    s = sprintf('BT /F1 %d Tf %.3f %.3f Td (%s) Tj ET\n', size, x, y, str);
end