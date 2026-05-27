clear all; close all; clc
addpath(fileparts(fileparts(mfilename('fullpath'))));

%% =========================================================
% Monolayer 1T'-WTe2: Shi-Song six-band k.p model + delta_z
%
% Clean structure:
%   - H, Vx=dH/dkx, Vy=dH/dky are used for Berry curvature and velocities.
%   - Stiffness uses vx, vy directly.
%   - Jerk uses only derivatives of vx, vy:
%         M_xx = d_kx vx,  M_xy = d_ky vx
%         M_yx = d_kx vy,  M_yy = d_ky vy
%
% Basis:
%   (psi1 up, psi1 down, psi2 up, psi2 down, psi3 up, psi3 down)
%
% Units:
%   energy   : meV
%   momentum : Angstrom^{-1}
%   velocity : meV Angstrom
%   M_ij     : meV Angstrom^2
%   Omega    : Angstrom^2
%   D        : Angstrom
% =========================================================

%% ------------------------- settings -----------------------
Nk       = 1801;              % use 401/801 for testing; 1801 is slow
kmax     = 0.75;             % Angstrom^{-1}
Temp_mu  = 10 * 0.08617;     % meV, 10 K
mu_grid  = linspace(-120,120,401);

gap_cutoff_BC = 1e-7;        % meV

p = SS6_parameters_meV();

% Gate-induced SOC. Default: delta_z-only limit.
p.delta1z = 25.0;            % meV
p.delta3z = 0.0;             % meV

fprintf('\n===== Shi-Song six-band: simple stiffness/jerk code =====\n');
fprintf('Nk = %d, kmax = %.3f A^{-1}, Temp_mu = %.5f meV\n', Nk, kmax, Temp_mu);
fprintf('delta1z = %.4f meV, delta3z = %.4f meV\n', p.delta1z, p.delta3z);

%% ------------------------- k grid -------------------------
kx_grid = linspace(-kmax,kmax,Nk);
ky_grid = linspace(-kmax,kmax,Nk);

dkx = kx_grid(2)-kx_grid(1);
dky = ky_grid(2)-ky_grid(1);
dk_area = dkx*dky/(2*pi)^2;

nb = 6;

E_grid     = zeros(Nk,Nk,nb);
Omega_grid = zeros(Nk,Nk,nb);
vx_grid    = zeros(Nk,Nk,nb);
vy_grid    = zeros(Nk,Nk,nb);

%% ------------------------- E, velocity, Berry curvature ----
textprogressbar('>> calculating E, velocity, Berry curvature: ');
for ix = 1:Nk
    textprogressbar(fix(100*ix/Nk));

    for iy = 1:Nk
        kx = kx_grid(ix);
        ky = ky_grid(iy);

        [H,Vx,Vy] = H_SS6_and_derivatives(kx,ky,p);

        [Vec,Eval] = eig(H);
        evals = real(diag(Eval));
        [evals,ord] = sort(evals,'ascend');
        Vec = Vec(:,ord);

        Vx_band = Vec' * Vx * Vec;
        Vy_band = Vec' * Vy * Vec;

        for a = 1:nb
            Ea = evals(a);
            Om = 0;

            for b = 1:nb
                if b == a
                    continue
                end

                dE = Ea - evals(b);

                if abs(dE) < gap_cutoff_BC
                    continue
                end

                Om = Om - 2*imag(Vx_band(a,b)*Vy_band(b,a))/(dE^2);
            end

            E_grid(ix,iy,a)     = Ea;
            Omega_grid(ix,iy,a) = Om;

            % Single-band velocity: v_i = <n|dH/dk_i|n> = d_i E_n
            vx_grid(ix,iy,a) = real(Vx_band(a,a));
            vy_grid(ix,iy,a) = real(Vy_band(a,a));
        end
    end
end
textprogressbar(' done.');

%% ------------------------- velocity derivatives for jerk ----
% Array convention:
%   dimension 1 = kx index = rows
%   dimension 2 = ky index = columns
%
% MATLAB gradient(F,dky,dkx):
%   first output  = derivative along columns = d/dky
%   second output = derivative along rows    = d/dkx

textprogressbar('>> precomputing velocity derivatives: ');
M_grid = zeros(Nk,Nk,nb,2,2);

for a = 1:nb
    textprogressbar(fix(100*a/nb));

    vx = vx_grid(:,:,a);
    vy = vy_grid(:,:,a);

    [dvx_dky, dvx_dkx] = gradient(vx, dky, dkx);
    [dvy_dky, dvy_dkx] = gradient(vy, dky, dkx);

    M_grid(:,:,a,1,1) = dvx_dkx;   % M_xx = d_kx v_x
    M_grid(:,:,a,1,2) = dvx_dky;   % M_xy = d_ky v_x
    M_grid(:,:,a,2,1) = dvy_dkx;   % M_yx = d_kx v_y
    M_grid(:,:,a,2,2) = dvy_dky;   % M_yy = d_ky v_y
end
textprogressbar(' done.');

%% ------------------------- band cut -------------------------
kx_cut = 0;
E_cut_y = zeros(nb,Nk);

for iy = 1:Nk
    H = H_SS6_and_derivatives(kx_cut,ky_grid(iy),p);
    E_cut_y(:,iy) = sort(real(eig(H)));
end

figure;
plot(ky_grid,E_cut_y.','LineWidth',1.2)
xlabel('$k_y$ [$\AA^{-1}$]','Interpreter','latex')
ylabel('$E$ [meV]','Interpreter','latex')
title('Shi--Song six-band cut, $k_x=0$','Interpreter','latex')
grid on; box on
ylim([-500 1000])

%% ------------------------- Berry map ------------------------
band_plot = 4;

figure;
imagesc(ky_grid,kx_grid,Omega_grid(:,:,band_plot))
set(gca,'YDir','normal')
axis equal tight
colorbar
xlabel('$k_y$ [$\AA^{-1}$]','Interpreter','latex')
ylabel('$k_x$ [$\AA^{-1}$]','Interpreter','latex')
title(sprintf('Berry curvature, band %d',band_plot),'Interpreter','latex')
box on

%% ------------------------- response vs mu -------------------
fFD  = @(e,mu) 1 ./ (exp((e-mu)./Temp_mu) + 1);
dfde = @(e,mu) -(1./(4*Temp_mu)) .* sech((e-mu)./(2*Temp_mu)).^2;

Dx = zeros(size(mu_grid));
Dy = zeros(size(mu_grid));
Omega_int = zeros(size(mu_grid));

S_photo = zeros(length(mu_grid),2,2);       % S(mu,gamma,j), meV
J_photo_native = zeros(length(mu_grid),2,2,2); % J(mu,i,j,gamma), meV Angstrom

textprogressbar('>> integrating response objects over mu_grid: ');
for imu = 1:length(mu_grid)
    textprogressbar(fix(100*imu/length(mu_grid)));

    mu = mu_grid(imu);

    D_tmp = zeros(2,1);
    S_tmp = zeros(2,2);
    J_tmp = zeros(2,2,2);
    Om_tmp = 0;

    for a = 1:nb
        En = E_grid(:,:,a);
        Om = Omega_grid(:,:,a);

        vx = vx_grid(:,:,a);
        vy = vy_grid(:,:,a);
        v = {vx, vy};

        f0 = fFD(En,mu);
        fp = dfde(En,mu);   % df/dE, negative

        % Berry curvature integral
        Om_tmp = Om_tmp + sum(f0 .* Om,'all') * dk_area;

        % BCD: D_i = - int Omega * d_i f = - int Omega * f'(E) * v_i
        D_tmp(1) = D_tmp(1) + sum(-Om .* fp .* vx,'all') * dk_area;
        D_tmp(2) = D_tmp(2) + sum(-Om .* fp .* vy,'all') * dk_area;

        % Stiffness: S^{gamma j} = - int f'(E) v_j v_gamma
        for gamma = 1:2
            for j = 1:2
                S_tmp(gamma,j) = S_tmp(gamma,j) - ...
                    sum(fp .* v{j} .* v{gamma},'all') * dk_area;
            end
        end

        % Jerk-like tensor: J^{i j gamma} = - int f'(E) v_j M_{i gamma}
        for i = 1:2
            for j = 1:2
                for gamma = 1:2
                    M_igamma = M_grid(:,:,a,i,gamma);
                    J_tmp(i,j,gamma) = J_tmp(i,j,gamma) - ...
                        sum(fp .* v{j} .* M_igamma,'all') * dk_area;
                end
            end
        end
    end

    Omega_int(imu) = Om_tmp;
    Dx(imu) = D_tmp(1);
    Dy(imu) = D_tmp(2);
    S_photo(imu,:,:) = S_tmp;
    J_photo_native(imu,:,:,:) = J_tmp;
end
textprogressbar(' done.');

D_photo = zeros(length(mu_grid),2);
D_photo(:,1) = 0.1 * Dx(:);       % Angstrom -> nm
D_photo(:,2) = 0.1 * Dy(:);       % Angstrom -> nm
Omega_z_photo = Omega_int(:);     % dimensionless
J_photo = 0.1 * J_photo_native;   % meV Angstrom -> meV nm

%% ------------------------- response plots -------------------
figure;
plot(mu_grid,Dx,'LineWidth',2); hold on
plot(mu_grid,Dy,'LineWidth',2)
xlabel('$\mu$ [meV]','Interpreter','latex')
ylabel('$D_i$ [$\AA$]','Interpreter','latex')
legend('$D_x$','$D_y$','Interpreter','latex','Location','best')
title('Berry curvature dipole','Interpreter','latex')
grid on; box on

Sxx = squeeze(S_photo(:,1,1));
Syy = squeeze(S_photo(:,2,2));
Jxxx = squeeze(J_photo(:,1,1,1));
Jyyy = squeeze(J_photo(:,2,2,2));

figure;
plot(mu_grid,Sxx,'LineWidth',2); hold on
plot(mu_grid,Syy,'LineWidth',2)
xlabel('$\mu$ [meV]','Interpreter','latex')
ylabel('$S$','Interpreter','latex')
legend('$S_{xx}$','$S_{yy}$','Interpreter','latex','Location','best')
title('Stiffness','Interpreter','latex')
grid on; box on



%% ------------------------- save -----------------------------
response_data = struct();

response_data.material_name = 'WTe2';
response_data.model_name = 'ShiSong_6band_deltaZ';
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

response_data.parameters = p;
response_data.parameters.Temp_mu = Temp_mu;
response_data.parameters.gap_cutoff_BC = gap_cutoff_BC;
response_data.parameters.native_units.energy = 'meV';
response_data.parameters.native_units.length = 'Angstrom';
response_data.parameters.native_units.momentum = 'Angstrom^{-1}';
response_data.parameters.native_units.D = 'Angstrom';
response_data.parameters.native_units.J = 'meV Angstrom';
response_data.parameters.native_units.saved_length_unit = 'nm';
response_data.parameters.native_units.Angstrom_to_nm = 0.1;

response_data.grid.native_momentum_unit = 'Angstrom^{-1}';
response_data.grid.saved_momentum_unit = 'nm^{-1}';
response_data.grid.Nk = Nk;
response_data.grid.kmax_A_inv = kmax;
response_data.grid.kmax_nm_inv = 10 * kmax;
response_data.grid.dkx_A_inv = dkx;
response_data.grid.dky_A_inv = dky;
response_data.grid.dk_area_A_inv2 = dk_area;

output_name = fullfile(fileparts(mfilename('fullpath')), ...
    'WTe2_ShiSong_6band_deltaZ_response_objects.mat');
validate_and_save_response_objects(output_name, response_data);

debug_name = fullfile(fileparts(mfilename('fullpath')), ...
    'WTe2_ShiSong_6band_deltaZ_debug_grids.mat');
save(debug_name, ...
    'p','kx_grid','ky_grid','E_grid','Omega_grid', ...
    'vx_grid','vy_grid','M_grid', ...
    'mu_grid','Dx','Dy','Omega_int','S_photo','J_photo_native', ...
    'Temp_mu','dk_area','gap_cutoff_BC', '-v7.3')

disp('Saved WTe2_ShiSong_6band_deltaZ_response_objects.mat')
disp('Saved optional debug grids to WTe2_ShiSong_6band_deltaZ_debug_grids.mat')

%% =========================================================
% Helper functions
% =========================================================

function p = SS6_parameters_meV()
    % Shi-Song Table A-III. Energies in meV, lengths in Angstrom.
    p.c10 =  1.0*1000;
    p.c20 =  0.0*1000;
    p.c30 = -0.4*1000;

    p.c1x = -11.25*1000;  p.c1y = -6.90*1000;
    p.c2x =  -0.27*1000;  p.c2y = -1.08*1000;
    p.c3x =  -0.82*1000;  p.c3y =  0.99*1000;

    p.v1x =  1.71*1000;   p.v1y =  0.48*1000;
    p.v3x =  0.48*1000;   p.v3y = -0.48*1000;

    % Gate-induced H1 parameters.
    p.lambda1 = 0; p.lambda2 = 0; p.lambda3 = 0; % meV A

    p.alpha1x = 0; p.alpha1y = 0; % meV A
    p.alpha2x = 0; p.alpha2y = 0;
    p.alpha3x = 0; p.alpha3y = 0;

    p.delta1x = 0; p.delta1z = 25; % meV
    p.delta3x = 0; p.delta3z = 0;  % meV
end

function [H,Vx,Vy] = H_SS6_and_derivatives(kx,ky,p)
    [H0,Vx0,Vy0] = H0_SS6(kx,ky,p);
    [H1,Vx1,Vy1] = H1_SS6(kx,ky,p);

    H  = H0  + H1;
    Vx = Vx0 + Vx1;
    Vy = Vy0 + Vy1;

    H  = (H+H')/2;
    Vx = (Vx+Vx')/2;
    Vy = (Vy+Vy')/2;
end

function [H,Vx,Vy] = H0_SS6(kx,ky,p)
    e1 = p.c10 + p.c1x*kx^2 + p.c1y*ky^2;
    e2 = p.c20 + p.c2x*kx^2 + p.c2y*ky^2;
    e3 = p.c30 + p.c3x*kx^2 + p.c3y*ky^2;

    de1x = 2*p.c1x*kx; de1y = 2*p.c1y*ky;
    de2x = 2*p.c2x*kx; de2y = 2*p.c2y*ky;
    de3x = 2*p.c3x*kx; de3y = 2*p.c3y*ky;

    v1p =  p.v1x*kx + 1i*p.v1y*ky;
    v1m = -p.v1x*kx + 1i*p.v1y*ky;
    v3p =  p.v3x*kx + 1i*p.v3y*ky;
    v3m = -p.v3x*kx + 1i*p.v3y*ky;

    dv1p_x =  p.v1x; dv1p_y = 1i*p.v1y;
    dv1m_x = -p.v1x; dv1m_y = 1i*p.v1y;
    dv3p_x =  p.v3x; dv3p_y = 1i*p.v3y;
    dv3m_x = -p.v3x; dv3m_y = 1i*p.v3y;

    H = [ ...
        e1,  0,  v1p,  0,    0,   0; ...
         0, e1,    0, v1m,   0,   0; ...
      conj(v1p),0, e2,  0,  v3p,  0; ...
         0,conj(v1m),0, e2,  0,  v3m; ...
         0,  0,conj(v3p),0, e3,  0; ...
         0,  0,   0,conj(v3m),0, e3 ...
    ];

    Vx = [ ...
        de1x, 0, dv1p_x, 0, 0, 0; ...
        0, de1x, 0, dv1m_x, 0, 0; ...
        conj(dv1p_x), 0, de2x, 0, dv3p_x, 0; ...
        0, conj(dv1m_x), 0, de2x, 0, dv3m_x; ...
        0, 0, conj(dv3p_x), 0, de3x, 0; ...
        0, 0, 0, conj(dv3m_x), 0, de3x ...
    ];

    Vy = [ ...
        de1y, 0, dv1p_y, 0, 0, 0; ...
        0, de1y, 0, dv1m_y, 0, 0; ...
        conj(dv1p_y), 0, de2y, 0, dv3p_y, 0; ...
        0, conj(dv1m_y), 0, de2y, 0, dv3m_y; ...
        0, 0, conj(dv3p_y), 0, de3y, 0; ...
        0, 0, 0, conj(dv3m_y), 0, de3y ...
    ];
end

function [H,Vx,Vy] = H1_SS6(kx,ky,p)
    a1p =  1i*p.alpha1x*kx + p.alpha1y*ky;
    a1m = -1i*p.alpha1x*kx + p.alpha1y*ky;
    a2p =  1i*p.alpha2x*kx + p.alpha2y*ky;
    a2m = -1i*p.alpha2x*kx + p.alpha2y*ky;
    a3p =  1i*p.alpha3x*kx + p.alpha3y*ky;
    a3m = -1i*p.alpha3x*kx + p.alpha3y*ky;

    da1p_x =  1i*p.alpha1x; da1p_y = p.alpha1y;
    da1m_x = -1i*p.alpha1x; da1m_y = p.alpha1y;
    da2p_x =  1i*p.alpha2x; da2p_y = p.alpha2y;
    da2m_x = -1i*p.alpha2x; da2m_y = p.alpha2y;
    da3p_x =  1i*p.alpha3x; da3p_y = p.alpha3y;
    da3m_x = -1i*p.alpha3x; da3m_y = p.alpha3y;

    L1 = p.lambda1*ky;
    L2 = p.lambda2*ky;
    L3 = p.lambda3*ky;

    H = [ ...
          L1,     a1m,   1i*p.delta1z,   1i*p.delta1x,              0,              0; ...
          a1p,    -L1,   1i*p.delta1x,  -1i*p.delta1z,              0,              0; ...
     -1i*p.delta1z, -1i*p.delta1x,     L2,     a2m,   1i*p.delta3z,   1i*p.delta3x; ...
     -1i*p.delta1x,  1i*p.delta1z,     a2p,    -L2,   1i*p.delta3x,  -1i*p.delta3z; ...
                0,              0, -1i*p.delta3z, -1i*p.delta3x,     L3,     a3m; ...
                0,              0, -1i*p.delta3x,  1i*p.delta3z,     a3p,    -L3 ...
    ];

    Vx = [ ...
        0, da1m_x, 0, 0, 0, 0; ...
        da1p_x, 0, 0, 0, 0, 0; ...
        0, 0, 0, da2m_x, 0, 0; ...
        0, 0, da2p_x, 0, 0, 0; ...
        0, 0, 0, 0, 0, da3m_x; ...
        0, 0, 0, 0, da3p_x, 0 ...
    ];

    Vy = [ ...
        p.lambda1, da1m_y, 0, 0, 0, 0; ...
        da1p_y, -p.lambda1, 0, 0, 0, 0; ...
        0, 0, p.lambda2, da2m_y, 0, 0; ...
        0, 0, da2p_y, -p.lambda2, 0, 0; ...
        0, 0, 0, 0, p.lambda3, da3m_y; ...
        0, 0, 0, 0, da3p_y, -p.lambda3 ...
    ];

    H  = (H+H')/2;
    Vx = (Vx+Vx')/2;
    Vy = (Vy+Vy')/2;
end
