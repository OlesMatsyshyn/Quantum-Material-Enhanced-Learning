clear all; close all; clc
addpath(fileparts(fileparts(mfilename('fullpath'))));

%% =========================================================
% Simple strained continuum AB / BA bilayer graphene
% Based on Fal'ko-style continuum Hamiltonian from the Xiong's note.
%
% Units:
%   energy   : meV
%   momentum : 1/a
%
% Basis:
%   AB: {A1, B1, A2, B2}
%   BA: {A1, B1, A2, B2}
%
% Low-energy convention:
%   pi = xi*p_x + i*p_y
%   v  = sqrt(3)/2 * gamma0
%   v3 = sqrt(3)/2 * gamma3
%
% Since p is measured in 1/a, v and v3 have units meV.
% =========================================================
stacking = 'AB';   % choose 'AB' or 'BA'
valley_list = [+1, -1];
spin_factor = 2;

%% =========================================================
% Parameters
% =========================================================

gamma0 = 2700;     % meV
gamma1 = 390;      % meV
gamma3 = 315;      % meV
U      = 5;       % meV, layer potential difference

a_g = 0.246;       % graphene lattice constant in nm

v0 = sqrt(3)/2 * gamma0;   % meV
v3 = sqrt(3)/2 * gamma3;   % meV

%% =========================================================
% Strain parameters
% =========================================================
% strain_amp = delta, and delta' = -nu*delta.
% theta_str is the principal strain axis relative to zigzag x direction.
%
% Note convention:
%   eta0 ~ -3
%   eta3 ~ -1
%
% rho = delta_r/a is the interlayer shear deformation associated with strain.
% From the note:
%   rho = 0.336*delta*(eta0-eta3)/eta3
%   phi = pi - 2 theta
% =========================================================

strain_amp = 0.01;      % delta, e.g. 1 percent
theta_str  = 0*pi/180;  % strain direction
nu         = 0.165;

eta0 = -3.0;
eta3 = -1.0;

delta_parallel = strain_amp;
delta_perp     = -nu * strain_amp;

delta_minus_delta_p = delta_parallel - delta_perp;

phi_shear = pi - 2*theta_str;
rho = 0.336 * strain_amp * (eta0 - eta3) / eta3;

S_total = [];
Omega_z_total = [];
D_total = [];
J_total = [];
Dx_arr_dfOmega_total = [];
Dy_arr_dfOmega_total = [];
A0_by_valley = zeros(size(valley_list));
A_by_valley = zeros(size(valley_list));
valley_results = struct();

for ivalley = 1:length(valley_list)
xi = valley_list(ivalley);

A0 = (3/4) * exp(-1i*2*xi*theta_str) ...
    * delta_minus_delta_p * eta0 * gamma0;

A = (3/4) * ( ...
    exp(-1i*2*xi*theta_str) * delta_minus_delta_p * (eta3 - eta0) ...
    - 2*sqrt(3) * exp(1i*xi*phi_shear) * rho * eta3 ...
    ) * gamma3;

A0_by_valley(ivalley) = A0;
A_by_valley(ivalley) = A;

fprintf('\n===== strained continuum BLG parameters =====\n');
fprintf('stacking = %s\n', stacking);
fprintf('xi = %+d\n', xi);
fprintf('v0 = %.6f meV\n', v0);
fprintf('v3 = %.6f meV\n', v3);
fprintf('A0 = %.6f %+ .6fi meV\n', real(A0), imag(A0));
fprintf('A  = %.6f %+ .6fi meV\n', real(A),  imag(A));
fprintf('rho = %.6e\n', rho);

%% =========================================================
% Momentum grid
% =========================================================

pmax = 0.25;
Np = 1801;

px_grid = linspace(-pmax,pmax,Np);
py_grid = linspace(-pmax,pmax,Np);

dpx = px_grid(2)-px_grid(1);
dpy = py_grid(2)-py_grid(1);

nbands = 4;

E_grid     = zeros(Np,Np,nbands);
Omega_grid = zeros(Np,Np,nbands);
vxdiag     = zeros(Np,Np,nbands);
vydiag     = zeros(Np,Np,nbands);

gap_cutoff = 1e-5;

%% =========================================================
% Eigensystem, velocities, Berry curvature
% =========================================================

textprogressbar('>> calculating E, velocity, Berry curvature: ');

for ix = 1:Np
    textprogressbar(fix(100*ix/Np));

    for iy = 1:Np
        px = px_grid(ix);
        py = py_grid(iy);

        [H,Vx,Vy] = H_and_V_simple_strained_BLG(px,py,xi, ...
            v0,v3,gamma1,U,A0,A,stacking);

        [Vec,Eval] = eig(H);
        evals = real(diag(Eval));

        [evals,ord] = sort(evals,'ascend');
        Vec = Vec(:,ord);

        Vx_band = Vec' * Vx * Vec;
        Vy_band = Vec' * Vy * Vec;

        for a = 2:3
            E_grid(ix,iy,a) = evals(a);
            vxdiag(ix,iy,a) = real(Vx_band(a,a));
            vydiag(ix,iy,a) = real(Vy_band(a,a));

            Om = 0;
            Ea = evals(a);

            for b = 1:nbands
                if b == a
                    continue
                end

                dE = evals(b) - Ea;

                if abs(dE) < gap_cutoff
                    continue
                end

                Om = Om - 2 * imag(Vx_band(a,b) * Vy_band(b,a)) / dE^2;
            end

            Omega_grid(ix,iy,a) = Om;
        end
    end
end

textprogressbar(' done.');

%% =========================================================
% 1D cut
% =========================================================

py = 0;
E_cut = zeros(4,Np);

for ix = 1:Np
    px = px_grid(ix);

    H = H_simple_strained_BLG(px,py,xi, ...
        v0,v3,gamma1,U,A0,A,stacking);

    E_cut(:,ix) = sort(real(eig(H)));
end

figure('Position',[650 100 500 360])
plot(px_grid,E_cut.','LineWidth',1.5)
xlabel('$p_x$ [$1/a$]','Interpreter','latex')
ylabel('$E$ [meV]','Interpreter','latex')
title(['BLG continuum cut, ' stacking ', $p_y=0$'],'Interpreter','latex')
ylim([-20 20])
grid on
box on

%% =========================================================
% Photocurrent objects
% =========================================================
% Momentum grid is in units of 1/a.
% Physical momentum is k_phys = p/a_g in nm^{-1}.
%
% Therefore:
%   dk_phys^2 = dp_x dp_y / a_g^2
%   Omega_phys = a_g^2 Omega_code
%   d/dk_phys = a_g d/dp_code
%
% Hence:
%   S       = int dk f d_k d_k E       -> no a_g conversion
%   Omega_z = int dk f Omega_xy        -> no a_g conversion
%   D       = int dk Omega_xy d_k f    -> multiply by a_g
%   J       = int dk f d_k d_k d_k E   -> multiply by a_g
% =========================================================

mu_grid = linspace(-10,10,301);
Temp_mu = 2 * 0.08617;   % meV

fFD   = @(e,mu) 1 ./ (exp((e-mu)./Temp_mu) + 1);
dfde  = @(e,mu) -(1./(4*Temp_mu)) .* sech((e-mu)./(2*Temp_mu)).^2;

dk_area_code = dpx * dpy / (2*pi)^2;

deg_factor_response = spin_factor;

S_photo       = zeros(length(mu_grid),2,2);       % S(mu,gamma,j)
Omega_z_photo = zeros(length(mu_grid),1);         % Omega_z(mu)
D_photo       = zeros(length(mu_grid),2);         % D(mu,i)
J_photo       = zeros(length(mu_grid),2,2,2);     % J(mu,i,j,gamma)

Dx_arr_dfOmega = zeros(size(mu_grid));
Dy_arr_dfOmega = zeros(size(mu_grid));
Dx_arr_fdOmega = zeros(size(mu_grid));
Dy_arr_fdOmega = zeros(size(mu_grid));

%% =========================================================
% Precompute derivatives of energy and Berry curvature
% =========================================================

dE1 = cell(nbands,2);
dE2 = cell(nbands,2,2);


dOm_dx = cell(nbands,1);
dOm_dy = cell(nbands,1);

for a = 2:3
    En = E_grid(:,:,a);
    Om = Omega_grid(:,:,a);

    dE1{a,1} = deriv_p_nonperiodic(En,dpx,dpy,1); % d_x E
    dE1{a,2} = deriv_p_nonperiodic(En,dpx,dpy,2); % d_y E

    for gamma = 1:2
        for j = 1:2
            dE2{a,gamma,j} = deriv_p_nonperiodic(dE1{a,gamma},dpx,dpy,j);
        end
    end

    dOm_dx{a} = deriv_p_nonperiodic(Om,dpx,dpy,1);
    dOm_dy{a} = deriv_p_nonperiodic(Om,dpx,dpy,2);
end

%% =========================================================
% Integrate response objects versus chemical potential
% =========================================================

textprogressbar('>> integrating response objects over mu_grid: ');
for imu = 1:length(mu_grid)
    textprogressbar(fix(100*imu/length(mu_grid)));

    mu = mu_grid(imu);

    S_tmp  = zeros(2,2);
    D_tmp  = zeros(2,1);
    J_tmp  = zeros(2,2,2);
    Om_tmp = 0;

    Dx_dfOmega_tmp = 0;
    Dy_dfOmega_tmp = 0;

    for a = 2:3
        En = E_grid(:,:,a);
        Om = Omega_grid(:,:,a);

        vx = vxdiag(:,:,a);
        vy = vydiag(:,:,a);

        f0 = fFD(En,mu);
        df_dE = dfde(En,mu);

        dfdx = df_dE .* vx;
        dfdy = df_dE .* vy;

        % Omega_z = int f Omega_xy
        Om_tmp = Om_tmp + sum(f0 .* Om,'all') * dk_area_code;

        % BCD consistency:
        % D = - int Omega * d_k f  = int f * d_k Omega
        Dx_dfOmega_tmp = Dx_dfOmega_tmp + sum((-Om .* dfdx),'all') * dk_area_code;
        Dy_dfOmega_tmp = Dy_dfOmega_tmp + sum((-Om .* dfdy),'all') * dk_area_code;

        % D^i = - int Omega_xy d_i f
        D_tmp(1) = D_tmp(1) + sum((-Om .* dfdx),'all') * dk_area_code;
        D_tmp(2) = D_tmp(2) + sum((-Om .* dfdy),'all') * dk_area_code;

        % S^{gamma j} = - int (d_j f) (d_gamma E)
        vcell = {vx, vy};
        
        for gamma = 1:2
            for j_dir = 1:2
                S_tmp(gamma,j_dir) = S_tmp(gamma,j_dir) - ...
                    sum(df_dE .* vcell{j_dir} .* vcell{gamma}, 'all') * dk_area_code;
            end
        end

        % J^{ijgamma} = - int (d_j f) (d_i d_gamma E)
        vcell = {vx, vy};
        
        for i_dir = 1:2
            for j_dir = 1:2
                for gamma = 1:2
                    J_tmp(i_dir,j_dir,gamma) = J_tmp(i_dir,j_dir,gamma) - ...
                        sum(df_dE .* vcell{j_dir} .* dE2{a,gamma,i_dir}, 'all') * dk_area_code;
                end
            end
        end
    end

    S_photo(imu,:,:)       = deg_factor_response * S_tmp;             % meV
    Omega_z_photo(imu)     = deg_factor_response * Om_tmp;            % dimensionless
    D_photo(imu,:)         = deg_factor_response * a_g * D_tmp.';     % nm
    J_photo(imu,:,:,:)     = deg_factor_response * a_g * J_tmp;       % meV nm

    Dx_arr_dfOmega(imu) = deg_factor_response * a_g * Dx_dfOmega_tmp;
    Dy_arr_dfOmega(imu) = deg_factor_response * a_g * Dy_dfOmega_tmp;
end
textprogressbar(' done.');

if isempty(S_total)
    S_total = zeros(size(S_photo));
    Omega_z_total = zeros(size(Omega_z_photo));
    D_total = zeros(size(D_photo));
    J_total = zeros(size(J_photo));
    Dx_arr_dfOmega_total = zeros(size(Dx_arr_dfOmega));
    Dy_arr_dfOmega_total = zeros(size(Dy_arr_dfOmega));
end

S_total = S_total + S_photo;
Omega_z_total = Omega_z_total + Omega_z_photo;
D_total = D_total + D_photo;
J_total = J_total + J_photo;
Dx_arr_dfOmega_total = Dx_arr_dfOmega_total + Dx_arr_dfOmega;
Dy_arr_dfOmega_total = Dy_arr_dfOmega_total + Dy_arr_dfOmega;

if xi == +1
    valley_key = 'valley_plus';
else
    valley_key = 'valley_minus';
end

valley_results.(valley_key).S_photo = S_photo;
valley_results.(valley_key).Omega_z_photo = Omega_z_photo;
valley_results.(valley_key).D_photo = D_photo;
valley_results.(valley_key).J_photo = J_photo;
valley_results.(valley_key).Dx_arr_dfOmega = Dx_arr_dfOmega;
valley_results.(valley_key).Dy_arr_dfOmega = Dy_arr_dfOmega;

end

S_photo = S_total;
Omega_z_photo = Omega_z_total;
D_photo = D_total;
J_photo = J_total;
Dx_arr_dfOmega = Dx_arr_dfOmega_total;
Dy_arr_dfOmega = Dy_arr_dfOmega_total;

fprintf('\n===== BLG valley cancellation diagnostics =====\n');
fprintf('max abs Omega_z for xi=+1: %.6e\n', max(abs(valley_results.valley_plus.Omega_z_photo)));
fprintf('max abs Omega_z for xi=-1: %.6e\n', max(abs(valley_results.valley_minus.Omega_z_photo)));
fprintf('max abs Omega_z for valley sum: %.6e\n', max(abs(Omega_z_photo)));
fprintf('max abs D for xi=+1: %.6e\n', max(sqrt(sum(valley_results.valley_plus.D_photo.^2, 2))));
fprintf('max abs D for xi=-1: %.6e\n', max(sqrt(sum(valley_results.valley_minus.D_photo.^2, 2))));
fprintf('max abs D for valley sum: %.6e\n', max(sqrt(sum(D_photo.^2, 2))));

%% =========================================================
% Save all response objects
% =========================================================

response_data = struct();

response_data.material_name = 'bilayer_graphene';
response_data.model_name = 'simple_strained_BLG_continuum';
response_data.version = 'v1_response_objects';

response_data.mu_grid = mu_grid;
response_data.S_photo = S_photo;
response_data.Omega_z_photo = Omega_z_photo;
response_data.D_photo = D_photo;
response_data.J_photo = J_photo;

response_data.Dx_arr_dfOmega = Dx_arr_dfOmega;
response_data.Dy_arr_dfOmega = Dy_arr_dfOmega;
response_data.valley_plus = valley_results.valley_plus;
response_data.valley_minus = valley_results.valley_minus;
response_data.diagnostics.Omega_z_xi_plus = valley_results.valley_plus.Omega_z_photo;
response_data.diagnostics.Omega_z_xi_minus = valley_results.valley_minus.Omega_z_photo;
response_data.diagnostics.D_xi_plus = valley_results.valley_plus.D_photo;
response_data.diagnostics.D_xi_minus = valley_results.valley_minus.D_photo;

response_data.units.energy = 'meV';
response_data.units.length = 'nm';
response_data.units.momentum = 'nm^{-1}';
response_data.units.S = 'meV';
response_data.units.Omega_z = 'dimensionless';
response_data.units.D = 'nm';
response_data.units.J = 'meV nm';

response_data.index_convention = ...
    '1=x, 2=y. S(mu,gamma,j), D(mu,i), J(mu,i,j,gamma).';
response_data.definitions.Omega_z = 'Omega_z(mu) = sum_n int_k f_n(k,mu) Omega_n^{xy}(k)';
response_data.definitions.D = 'D^i(mu) = - sum_n int_k Omega_n^{xy}(k) partial_i f_n(k,mu)';
response_data.definitions.S = 'S^{gamma j}(mu) = - sum_n int_k f''_n(E) v_n^j v_n^gamma';
response_data.definitions.J = 'J^{i j gamma}(mu) = - sum_n int_k f''_n(E) v_n^j partial_i partial_gamma E_n';

response_data.parameters.stacking = stacking;
response_data.parameters.valley_list = valley_list;
response_data.parameters.spin_factor = spin_factor;
response_data.parameters.valley_summed_explicitly = true;
response_data.parameters.gamma0 = gamma0;
response_data.parameters.gamma1 = gamma1;
response_data.parameters.gamma3 = gamma3;
response_data.parameters.U = U;
response_data.parameters.v0 = v0;
response_data.parameters.v3 = v3;
response_data.parameters.a_g_nm = a_g;
response_data.parameters.strain_amp = strain_amp;
response_data.parameters.theta_str = theta_str;
response_data.parameters.nu = nu;
response_data.parameters.eta0 = eta0;
response_data.parameters.eta3 = eta3;
response_data.parameters.rho = rho;
response_data.parameters.phi_shear = phi_shear;
response_data.parameters.A0_by_valley = A0_by_valley;
response_data.parameters.A_by_valley = A_by_valley;
response_data.parameters.Np = Np;
response_data.parameters.pmax = pmax;
response_data.parameters.deg_factor_response = spin_factor;

response_data.grid.coordinate = 'p_code';
response_data.grid.p_units = '1/a_g';
response_data.grid.saved_momentum_unit = 'nm^{-1}';
response_data.grid.pmax = pmax;
response_data.grid.Np = Np;
response_data.grid.dpx = dpx;
response_data.grid.dpy = dpy;
response_data.grid.dk_area_code = dk_area_code;
response_data.grid.k_nm_inv_equals_p_over_a_g_nm = true;

output_name = fullfile(fileparts(mfilename('fullpath')), ...
    ['simple_strained_BLG_photocurrent_objects_' stacking '_valleysum.mat']);
validate_and_save_response_objects(output_name, response_data);

disp('Saved valley-summed simple strained BLG photocurrent objects.');

%% =========================================================
% Final plots
% =========================================================

% 1) BCD consistency check
figure('Position',[100 100 650 420]);

plot(mu_grid, Dx_arr_dfOmega, 'r-',  'LineWidth', 2); hold on
plot(mu_grid, Dy_arr_dfOmega, 'b-',  'LineWidth', 2);

xlabel('\mu [meV]')
ylabel('D [nm]')
legend( ...
    'D_x: -\Omega \partial_k f', ...
    'D_y: -\Omega \partial_k f', ...
    'Location','best')
title(['BCD consistency check, ' stacking])
box on
grid on

% 2) left axis S, right axis J
Sxx = squeeze(S_photo(:,1,1));
Sxy = squeeze(S_photo(:,1,2));
Syy = squeeze(S_photo(:,2,2));

Jxxx = squeeze(J_photo(:,1,1,1));
Jxxy = squeeze(J_photo(:,1,1,2));
Jyyy = squeeze(J_photo(:,2,2,2));

figure('Position',[150 120 760 440]);

yyaxis left
plot(mu_grid, Sxx, '-', 'LineWidth', 2); hold on
plot(mu_grid, Syy, '-', 'LineWidth', 2);
plot(mu_grid, Sxy, '--', 'LineWidth', 2);
ylabel('S [meV]')

yyaxis right
plot(mu_grid, Jxxx, '-', 'LineWidth', 2); hold on
plot(mu_grid, Jyyy, '-', 'LineWidth', 2);
plot(mu_grid, Jxxy, '--', 'LineWidth', 2);
ylabel('J [meV nm]')

xlabel('\mu [meV]')
title(['Left axis: S, Right axis: J, ' stacking])
legend( ...
    'S^{xx}','S^{yy}','S^{xy}', ...
    'J^{xxx}','J^{yyy}','J^{xxy}', ...
    'Location','best')
box on
grid on

%% =========================================================
% Helper functions
% =========================================================

function [H,Vx,Vy] = H_and_V_simple_strained_BLG(px,py,xi, ...
    v0,v3,gamma1,U,A0,A,stacking)

    H = H_simple_strained_BLG(px,py,xi,v0,v3,gamma1,U,A0,A,stacking);

    dpi_dx = xi;
    dpi_dy = 1i;

    dpidag_dx = xi;
    dpidag_dy = -1i;

    switch upper(stacking)

        case 'AB'
            Vx = [ ...
                0,              v0*dpidag_dx, 0,              v3*dpi_dx;
                v0*dpi_dx,      0,            0,              0;
                0,              0,            0,              v0*dpidag_dx;
                v3*dpidag_dx,   0,            v0*dpi_dx,      0 ];

            Vy = [ ...
                0,              v0*dpidag_dy, 0,              v3*dpi_dy;
                v0*dpi_dy,      0,            0,              0;
                0,              0,            0,              v0*dpidag_dy;
                v3*dpidag_dy,   0,            v0*dpi_dy,      0 ];

        case 'BA'
            Vx = [ ...
                0,              v0*dpidag_dx, 0,              0;
                v0*dpi_dx,      0,            v3*dpidag_dx,   0;
                0,              v3*dpi_dx,    0,              v0*dpidag_dx;
                0,              0,            v0*dpi_dx,      0 ];

            Vy = [ ...
                0,              v0*dpidag_dy, 0,              0;
                v0*dpi_dy,      0,            v3*dpidag_dy,   0;
                0,              v3*dpi_dy,    0,              v0*dpidag_dy;
                0,              0,            v0*dpi_dy,      0 ];

        otherwise
            error('stacking must be AB or BA')
    end
end

function H = H_simple_strained_BLG(px,py,xi, ...
    v0,v3,gamma1,U,A0,A,stacking)

    pi0    = xi*px + 1i*py;
    pidag0 = xi*px - 1i*py;

    pit    = pi0    + A0/v0;
    pitdag = pidag0 + conj(A0)/v0;

    switch upper(stacking)

        case 'AB'
            H = [ ...
                U/2,        v0*pitdag,       0,              v3*pit + A;
                v0*pit,     U/2,             gamma1,         0;
                0,          gamma1,         -U/2,            v0*pitdag;
                v3*pitdag + conj(A), 0,      v0*pit,        -U/2 ];

        case 'BA'
            % BA continuum analogue of the note's Eq. (22), with strain
            % inserted into pi -> pi_tilde and skew term shifted by A.
            H = [ ...
                U/2,        v0*pitdag,       0,              gamma1;
                v0*pit,     U/2,             v3*pitdag + conj(A), 0;
                0,          v3*pit + A,     -U/2,            v0*pitdag;
                gamma1,     0,               v0*pit,        -U/2 ];

        otherwise
            error('stacking must be AB or BA')
    end
end

function dFdp = deriv_p_nonperiodic(F,dpx,dpy,dir)
    % Non-periodic finite difference on a rectangular Cartesian grid.
    % dir = 1: d/dp_x
    % dir = 2: d/dp_y

    dFdp = zeros(size(F));

    switch dir
        case 1
            dFdp(2:end-1,:) = (F(3:end,:) - F(1:end-2,:)) / (2*dpx);
            dFdp(1,:)       = (F(2,:)     - F(1,:))       / dpx;
            dFdp(end,:)     = (F(end,:)   - F(end-1,:))   / dpx;

        case 2
            dFdp(:,2:end-1) = (F(:,3:end) - F(:,1:end-2)) / (2*dpy);
            dFdp(:,1)       = (F(:,2)     - F(:,1))       / dpy;
            dFdp(:,end)     = (F(:,end)   - F(:,end-1))   / dpy;

        otherwise
            error('dir must be 1 or 2')
    end
end
