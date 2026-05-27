clear all; close all; clc
addpath(fileparts(fileparts(mfilename('fullpath'))));

%% =========================================================
% MBT response-object generator
% Based on BandStructure.m, but saves only compact response tensors.
%
% Saved convention:
%   energy   : meV
%   length   : nm
%   momentum : nm^{-1}
% =========================================================

%% ------------------------- settings -----------------------
sizekx = 800;
sizeky = 800;

Lx_A_inv = 0.05;
Ly_A_inv = 0.05;
Lx = 10 * Lx_A_inv;  % nm^{-1}
Ly = 10 * Ly_A_inv;  % nm^{-1}

kx = linspace(-Lx,Lx,sizekx);
ky = linspace(-Ly,Ly,sizeky);
dkx = kx(2)-kx(1);
dky = ky(2)-ky(1);
dk_area = dkx*dky/(2*pi)^2;

mu_grid = linspace(-75,75,401); % meV
Temp_mu = 0.1;                    % meV, from old Temp = 0.0001 eV
gap_cutoff_BC = 1e-9;             % meV
deg_tol = 1e-8;                   % meV, omit internal degenerate doublet denominators
deg_factor_response = 1;

%% ------------------------- parameters ---------------------
% Native values from BandStructure.m.
alpha_eV_A = 3.2;
gamma_eV_A2 = 17;
beta0_eV_A2 = 9.4;
m0_eV = 0.04;
m1_eV = 0.05;
m2_eV = 0.09;
t1_eV = -0.0533;
t2_eV = 0.0463;
lambda_eV = 0.0557;
t1p_eV = 0.00073;
t2p_eV = 0.0011;
lambdap_eV = 0.002;
neel_order = 1;
amp = 1;
E0_eV = 0;
mass_anisotropy = 0.01;  % demonstration-only phenomenological anisotropy

% Convert Hamiltonian to meV and nm^{-1}.
alpha = alpha_eV_A * 100;     % meV nm
gamma = gamma_eV_A2 * 10;     % meV nm^2
beta0 = beta0_eV_A2 * 10;     % meV nm^2
m0 = m0_eV * 1000;            % meV
m1 = m1_eV * 1000;            % meV
m2 = m2_eV * 1000;            % meV
t1 = t1_eV * 1000;            % meV
t2 = t2_eV * 1000;            % meV
lambda = lambda_eV * 1000;    % meV
t1p = t1p_eV * 1000;          % meV, retained for provenance
t2p = t2p_eV * 1000;          % meV, retained for provenance
lambdap = lambdap_eV * 1000;  % meV, retained for provenance
E0 = E0_eV * 1000;            % meV
sx_mass = 1 + mass_anisotropy;
sy_mass = 1 - mass_anisotropy;

%% ------------------------- Hamiltonian --------------------
% This is a phenomenological anisotropic effective-mass deformation for
% demonstration and symmetry testing. It is not a microscopic fitted strain
% Hamiltonian for MBT.
kp = @(px,py) px + 1i*py;
km = @(px,py) px - 1i*py;
e0 = @(px,py) gamma*(sx_mass*px.^2 + sy_mass*py.^2);
mi = @(px,py) m0 + beta0*(sx_mass*px.^2 + sy_mass*py.^2);

hN = @(px,py) [e0(px,py) + mi(px,py),   alpha*km(px,py),        0,                       0                     ;
               alpha*kp(px,py),         e0(px,py) - mi(px,py),  0,                       0                     ;
               0,                       0,                      e0(px,py) - mi(px,py),   alpha*km(px,py)       ;
               0,                       0,                      alpha*kp(px,py),         e0(px,py) + mi(px,py)];

hAFM = [m1,  0  ,  0 ,  0   ;
        0 ,  -m2,  0 ,  0   ;
        0 ,  0  ,  m2,  0   ;
        0 ,  0  ,  0 ,  -m1] * neel_order;

h1 = @(px,py) hN(px,py) + 0.5*hAFM;
h2 = @(px,py) hN(px,py) - 0.5*hAFM;

T0 = [t1,         0,           -1i*lambda,  0           ;
      0,          t2,          0,           1i*lambda   ;
      -1i*lambda, 0,           t2,          0           ;
      0,          1i*lambda,   0,           t1         ];

T = amp * T0;
Ddisp = E0 * eye(4);

H = @(px,py) [h1(px,py)+Ddisp,  T               ;
              T',               h2(px,py)-Ddisp];

dxh = @(px,py)[2*px*sx_mass*(gamma + beta0),    alpha,                             0,                                 0                    ;
               alpha,                          2*px*sx_mass*(gamma - beta0),      0,                                 0                    ;
               0,                              0,                                 2*px*sx_mass*(gamma - beta0),      alpha                ;
               0,                              0,                                 alpha,                             2*px*sx_mass*(gamma + beta0)];

dyh = @(px,py)[2*py*sy_mass*(gamma + beta0),    -1i*alpha,                         0,                                 0                    ;
               1i*alpha,                       2*py*sy_mass*(gamma - beta0),      0,                                 0                    ;
               0,                              0,                                 2*py*sy_mass*(gamma - beta0),      -1i*alpha            ;
               0,                              0,                                 1i*alpha,                          2*py*sy_mass*(gamma + beta0)];

Vx = @(px,py) [dxh(px,py), zeros(4,4); zeros(4,4), dxh(px,py)];
Vy = @(px,py) [dyh(px,py), zeros(4,4); zeros(4,4), dyh(px,py)];

nbands = 8;
E_grid = zeros(sizekx,sizeky,nbands);
Omega_grid = zeros(sizekx,sizeky,nbands);
vx_grid = zeros(sizekx,sizeky,nbands);
vy_grid = zeros(sizekx,sizeky,nbands);

textprogressbar('>> calculating MBT E, velocity, Berry curvature: ');
for ix = 1:sizekx
    textprogressbar(fix(100*ix/sizekx));

    for iy = 1:sizeky
        [Vec,Eval] = eig(H(kx(ix),ky(iy)));
        evals = real(diag(Eval));
        [evals,ord] = sort(evals,'ascend');
        Vec = Vec(:,ord);

        Vx_band = Vec' * Vx(kx(ix),ky(iy)) * Vec;
        Vy_band = Vec' * Vy(kx(ix),ky(iy)) * Vec;

        for a = 1:nbands
            Ea = evals(a);
            Om = 0;

            for b = 1:nbands
                if b == a
                    continue
                end

                dE = Ea - evals(b);
                if abs(evals(b) - Ea) < deg_tol
                    continue
                end

                Om = Om - 2*imag(Vx_band(a,b)*Vy_band(b,a))/(dE^2);
            end

            E_grid(ix,iy,a) = Ea;
            Omega_grid(ix,iy,a) = Om;
            vx_grid(ix,iy,a) = real(Vx_band(a,a));
            vy_grid(ix,iy,a) = real(Vy_band(a,a));
        end
    end
end
textprogressbar(' done.');

%% ------------------------- doublet diagnostics ------------
doublet_pairs = [1 2; 3 4; 5 6; 7 8];
max_doublet_splitting = zeros(size(doublet_pairs,1),1);

fprintf('\n===== MBT doublet splitting diagnostics =====\n');
for ipair = 1:size(doublet_pairs,1)
    n1 = doublet_pairs(ipair,1);
    n2 = doublet_pairs(ipair,2);
    split_grid = abs(E_grid(:,:,n2) - E_grid(:,:,n1));
    max_doublet_splitting(ipair) = max(split_grid, [], 'all');
    fprintf('pair [%d %d]: max splitting = %.6e meV\n', ...
        n1, n2, max_doublet_splitting(ipair));
end

if any(max_doublet_splitting > 100 * deg_tol)
    warning('MBT:LargeDoubletSplitting', ...
        'At least one expected doublet splitting exceeds 100*deg_tol.');
end

%% ------------------------- velocity derivatives -----------
textprogressbar('>> precomputing MBT velocity derivatives: ');
M_grid = zeros(sizekx,sizeky,nbands,2,2);
for a = 1:nbands
    textprogressbar(fix(100*a/nbands));

    vx = vx_grid(:,:,a);
    vy = vy_grid(:,:,a);

    [dvx_dky, dvx_dkx] = gradient(vx, dky, dkx);
    [dvy_dky, dvy_dkx] = gradient(vy, dky, dkx);

    M_grid(:,:,a,1,1) = dvx_dkx;
    M_grid(:,:,a,1,2) = dvy_dkx;
    M_grid(:,:,a,2,1) = dvx_dky;
    M_grid(:,:,a,2,2) = dvy_dky;
end
textprogressbar(' done.');

%% ------------------------- response objects ---------------
fFD  = @(e,mu) 1 ./ (exp((e-mu)./Temp_mu) + 1);
dfde = @(e,mu) -(1./(4*Temp_mu)) .* sech((e-mu)./(2*Temp_mu)).^2;

S_photo = zeros(length(mu_grid),2,2);
Omega_z_photo = zeros(length(mu_grid),1);
D_photo = zeros(length(mu_grid),2);
J_photo = zeros(length(mu_grid),2,2,2);

textprogressbar('>> integrating response objects over mu_grid: ');
for imu = 1:length(mu_grid)
    textprogressbar(fix(100*imu/length(mu_grid)));

    mu = mu_grid(imu);

    S_tmp = zeros(2,2);
    D_tmp = zeros(2,1);
    J_tmp = zeros(2,2,2);
    Om_tmp = 0;

    for a = 1:nbands
        En = E_grid(:,:,a);
        Om = Omega_grid(:,:,a);
        vx = vx_grid(:,:,a);
        vy = vy_grid(:,:,a);
        v = {vx, vy};

        f0 = fFD(En,mu);
        fp = dfde(En,mu);

        Om_tmp = Om_tmp + sum(f0 .* Om,'all') * dk_area;

        D_tmp(1) = D_tmp(1) + sum(-Om .* fp .* vx,'all') * dk_area;
        D_tmp(2) = D_tmp(2) + sum(-Om .* fp .* vy,'all') * dk_area;

        for gamma_dir = 1:2
            for j_dir = 1:2
                S_tmp(gamma_dir,j_dir) = S_tmp(gamma_dir,j_dir) - ...
                    sum(fp .* v{j_dir} .* v{gamma_dir},'all') * dk_area;
            end
        end

        for i_dir = 1:2
            for j_dir = 1:2
                for gamma_dir = 1:2
                    M_igamma = M_grid(:,:,a,i_dir,gamma_dir);
                    J_tmp(i_dir,j_dir,gamma_dir) = J_tmp(i_dir,j_dir,gamma_dir) - ...
                        sum(fp .* v{j_dir} .* M_igamma,'all') * dk_area;
                end
            end
        end
    end

    S_photo(imu,:,:) = deg_factor_response * S_tmp;
    Omega_z_photo(imu) = deg_factor_response * Om_tmp;
    D_photo(imu,:) = deg_factor_response * D_tmp.';
    J_photo(imu,:,:,:) = deg_factor_response * J_tmp;
end
textprogressbar(' done.');

%% ------------------------- diagnostic plots ----------------
% These plots use already computed arrays only and do not affect saved data.
[~, iy0] = min(abs(ky));

figure('Name','MBT band slice');
plot(kx, squeeze(E_grid(:,iy0,:)), 'LineWidth', 1.1)
xlabel('$k_x$ [nm$^{-1}$]', 'Interpreter','latex')
ylabel('$E$ [meV]', 'Interpreter','latex')
title('MBT band slice, ky = 0')
grid on; box on

figure('Name','MBT Berry curvature dipole');
plot(mu_grid, D_photo(:,1), 'LineWidth', 2); hold on
plot(mu_grid, D_photo(:,2), 'LineWidth', 2)
xlabel('$\mu$ [meV]', 'Interpreter','latex')
ylabel('$D$ [nm]', 'Interpreter','latex')
legend('$D_x$','$D_y$', 'Interpreter','latex', 'Location','best')
title('MBT Berry curvature dipole')
grid on; box on

figure('Name','MBT S and J response objects');
Sxx = squeeze(S_photo(:,1,1));
Syy = squeeze(S_photo(:,2,2));
Sxy = squeeze(S_photo(:,1,2));
Jxxx = squeeze(J_photo(:,1,1,1));
Jyyy = squeeze(J_photo(:,2,2,2));
Jxxy = squeeze(J_photo(:,1,1,2));

yyaxis left
plot(mu_grid, Sxx, 'LineWidth', 2); hold on
plot(mu_grid, Syy, 'LineWidth', 2)
plot(mu_grid, Sxy, '--', 'LineWidth', 2)
ylabel('$S$ [meV]', 'Interpreter','latex')

yyaxis right
plot(mu_grid, Jxxx, 'LineWidth', 2); hold on
plot(mu_grid, Jyyy, 'LineWidth', 2)
plot(mu_grid, Jxxy, '--', 'LineWidth', 2)
ylabel('$J$ [meV nm]', 'Interpreter','latex')

xlabel('$\mu$ [meV]', 'Interpreter','latex')
legend('$S_{xx}$','$S_{yy}$','$S_{xy}$', ...
       '$J_{xxx}$','$J_{yyy}$','$J_{xxy}$', ...
       'Interpreter','latex', 'Location','best')
title('MBT stiffness and jerk response objects')
grid on; box on

figure('Name','MBT Omega_z response object');
plot(mu_grid, Omega_z_photo, 'LineWidth', 2)
xlabel('$\mu$ [meV]', 'Interpreter','latex')
ylabel('$\Omega_z$ [dimensionless]', 'Interpreter','latex')
title('MBT integrated Berry curvature')
grid on; box on

%% ------------------------- save ---------------------------
parameters = struct();
parameters.native_units.energy = 'eV';
parameters.native_units.length = 'Angstrom';
parameters.native_units.momentum = 'Angstrom^{-1}';
parameters.native_units.alpha = 'eV Angstrom';
parameters.native_units.gamma = 'eV Angstrom^2';
parameters.native_units.beta0 = 'eV Angstrom^2';
parameters.native_units.eV_to_meV = 1000;
parameters.native_units.Angstrom_to_nm = 0.1;
parameters.alpha_eV_A = alpha_eV_A;
parameters.gamma_eV_A2 = gamma_eV_A2;
parameters.beta0_eV_A2 = beta0_eV_A2;
parameters.m0_eV = m0_eV;
parameters.m1_eV = m1_eV;
parameters.m2_eV = m2_eV;
parameters.t1_eV = t1_eV;
parameters.t2_eV = t2_eV;
parameters.lambda_eV = lambda_eV;
parameters.t1p_eV = t1p_eV;
parameters.t2p_eV = t2p_eV;
parameters.lambdap_eV = lambdap_eV;
parameters.neel_order = neel_order;
parameters.amp = amp;
parameters.E0_eV = E0_eV;
parameters.alpha_meV_nm = alpha;
parameters.gamma_meV_nm2 = gamma;
parameters.beta0_meV_nm2 = beta0;
parameters.m0_meV = m0;
parameters.m1_meV = m1;
parameters.m2_meV = m2;
parameters.t1_meV = t1;
parameters.t2_meV = t2;
parameters.lambda_meV = lambda;
parameters.mass_anisotropy = mass_anisotropy;
parameters.mass_anisotropy_model = ...
    'demonstration-only: kx^2+ky^2 -> (1+eta)kx^2 + (1-eta)ky^2 in e0 and mi only';
parameters.Temp_mu = Temp_mu;
parameters.gap_cutoff_BC = gap_cutoff_BC;
parameters.degeneracy_handling = ...
    'Berry/Hall denominators with |E_m-E_n| < deg_tol omitted';
parameters.deg_tol_meV = deg_tol;
parameters.expected_doublet_pairs = doublet_pairs;
parameters.max_doublet_splitting_meV = max_doublet_splitting;
parameters.deg_factor_response = deg_factor_response;

grid_info = struct();
grid_info.coordinate = 'Cartesian k';
grid_info.saved_momentum_unit = 'nm^{-1}';
grid_info.sizekx = sizekx;
grid_info.sizeky = sizeky;
grid_info.Lx_nm_inv = Lx;
grid_info.Ly_nm_inv = Ly;
grid_info.Lx_A_inv_native = Lx_A_inv;
grid_info.Ly_A_inv_native = Ly_A_inv;
grid_info.dkx = dkx;
grid_info.dky = dky;
grid_info.dk_area = dk_area;

response_data = struct();
response_data.material_name = 'MBT';
response_data.model_name = 'AB_stacked_MBT_8band';
response_data.version = 'v1_response_objects';
response_data.mu_grid = mu_grid;
response_data.S_photo = S_photo;
response_data.Omega_z_photo = Omega_z_photo;
response_data.D_photo = D_photo;
response_data.J_photo = J_photo;
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
response_data.parameters = parameters;
response_data.grid = grid_info;

output_name = fullfile(fileparts(mfilename('fullpath')), 'MBT_response_objects.mat');
validate_and_save_response_objects(output_name, response_data);

disp(['Saved ' output_name])
