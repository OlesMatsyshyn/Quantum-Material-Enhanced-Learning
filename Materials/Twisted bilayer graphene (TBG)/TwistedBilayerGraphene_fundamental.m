clear all
addpath(fileparts(fileparts(mfilename('fullpath'))));
% Ref[1] is Phys. Rev. B 106, L041111 
% Giant nonlinear Hall effect in strained twisted bilayer graphene

% Parameters, constants 
theta  = 1.2 * pi / 180; % twist angle in radians
valley_list = [+1, -1];
spin_factor = 2;
alpha  = 2422.7;          % meV, (2135.4) is the hbar*v_F/a; look after Eq.(2) in 
                          % section EFFECTIVE CONTINUUM MODEL
a_g = 0.246; % graphene lattice constant in nm

eps_strain = 0.001;   % heterostrain magnitude, e.g. 0.1% = 0.001
phi_strain = 0*pi/2;  % strain direction (radians): 0 zigzag; pi/2 armchair
nu         = 0.165;   % graphene Poisson ratio
beta_gr    = 1.57;    % graphene strain gauge factor

Delta1 = 17;           % hBN staggered potential on layer 1
Delta2 = 0;            % hBN staggered potential on layer 2
Temp_mu = 2 * 0.08617; % (n) K in meV
% Select the sutable number. Controlls the accuracy of the model.
% N_moire is a positive integer. Complexity grows fast with
% N_moire. For N_moire = 10, code takes a few minutes to execute.
N_moire = 4; % indicates how many cites (forward and backward) 
             % in k space are coupled at given k_now 
             % resulting H, will have 4*(2*N_moire+1)^2 bands
N_side = 4;  % keep +- N_side bands around the two middle bands; for speed
% Density of the Brilluine zone
Nx = 321;
Ny = 323;

% Pauli Matrises
sx = [0, 1 ;
      1, 0];
sy = [ 0, -1i ;
      1i,   0];
sz = [1, 0 ;
      0,-1];
s0 = eye(2); 

% note a = 0.246 nm
% reciprocal lattice vectors in units of 1/a
b1 = 2*pi*[-1; -1/sqrt(3)]; b2 = 2*pi*[-1; 1/sqrt(3)]; 

% Anti-clock wise rotation by angle th in 2D
R2D = @(th)[ cos(th), -sin(th)  ;
             sin(th),  cos(th) ];
% Create the rotational matrises with values
R1 = R2D(-theta/2);   % bottom-layer geometry
R2 = R2D(+theta/2);   % top-layer geometry

% Uniaxial strain tensor
Estr = eps_strain * [ ...
     cos(phi_strain)^2 - nu*sin(phi_strain)^2, (1+nu)*cos(phi_strain)*sin(phi_strain); ...
     (1+nu)*cos(phi_strain)*sin(phi_strain)  , -nu*cos(phi_strain)^2 + sin(phi_strain)^2 ];

% Heterostrain: opposite strain in the two layers
% To match the paper convention, use layer 1 = -E/2, layer 2 = +E/2
E1 = Estr;        % bottom layer strained
E2 = zeros(2);    % top layer unstrained

% Strain-induced gauge fields
A1 = sqrt(3) * beta_gr * [E1(1,1)-E1(2,2); -2*E1(1,2)];
A2 = [0;0];

S_total = [];
Omega_z_total = [];
D_total = [];
J_total = [];
Dx_arr_dfOmega_nm_total = [];
Dy_arr_dfOmega_nm_total = [];
Dx_arr_fdOmega_nm_total = [];
Dy_arr_fdOmega_nm_total = [];
valley_results = struct();

for ivalley = 1:length(valley_list)
valley = valley_list(ivalley);

% =========================================================
% Giant nonlinear Hall effect in strained twisted bilayer graphene geometry
% =========================================================
% Points on the effective (dimentionless) moiré Brilluin zone
%                 top corner
%                   /  \
%                  /    \ 
%                 /      \
%     Kt_left ---          --- Kt_right
%           |                  |
%           |        Γ         |
%        M_left              M_right
%           |                  |
%           |                  |
%      Kb_left ---          --- Kb_right
%                 \      /
%                  \    / 
%                   \  /
%                bottom corner
K0 = -(b1 + b2)/3;   % +K, along +x
% Layer Dirac points before gauge field
% Apparently we strain then rotate
K_b0 = (eye(2)-E1') * K0;   % bottom layer
K_t0 = (eye(2)-E2') * K0;   % top layer

Kb_geom = valley * R1 * K_b0;   % bottom geometric Dirac point, no A
Kt_geom = valley * R2 * K_t0;   % top geometric Dirac point, no A

% Gauge-shifted Dirac points before the rotation
Kb_right = R1*(valley*K_b0 - valley*A1);   % bottom
Kt_right = R2*(valley*K_t0 - valley*A2);   % top

% SI momentum transfers
qb  = R1*(eye(2)-E1')*K0        - R2*(eye(2)-E2')*K0;
qtl = R1*(eye(2)-E1')*(K0 + b1) - R2*(eye(2)-E2')*(K0 + b1);
qtr = R1*(eye(2)-E1')*(K0 + b2) - R2*(eye(2)-E2')*(K0 + b2);

% Moire reciprocal vectors defined by SI transfer differences
Gtl = qtl - qb;   % associated with T3 = T_qtl
Gtr = qtr - qb;   % associated with T2 = T_qtr

% Use code names:
G1M = Gtl;
G2M = Gtr;

% The old horizontal moire vector is therefore:
Gpath = G1M - G2M; % negative to put Gamma to the left wrt Kt/b

% Gamma is the center of the moire BZ cell
M_right = 0.5 * (Kt_right + Kb_right);
Gamma  = M_right + 0.5 * Gpath;

Kb_left = Kb_right + Gpath;   % bottom
Kt_left = Kt_right + Gpath;   % top

% Path: K_t -> K_b -> Gamma -> K_t
N1 = 80;
N2 = 80;
N3 = 80;

k_grid = [];
k_grid = path_add(k_grid, Kt_left,   Kb_left,    N1);
k_grid = path_add(k_grid, Kb_left,   Gamma, N2);
k_grid = path_add(k_grid, Gamma, Kt_left,   N3);

k_length = size(k_grid,2);

fprintf('\n===== geometry check =====\n');
fprintf('G1M = [% .6f, % .6f]\n', G1M(1), G1M(2));
fprintf('G2M = [% .6f, % .6f]\n', G2M(1), G2M(2));
fprintf('Gpath = [% .6f, % .6f]\n', Gpath(1), Gpath(2));
fprintf('Kb right = [% .6f, % .6f]\n', Kb_right(1), Kb_right(2));
fprintf('Kt right = [% .6f, % .6f]\n', Kt_right(1), Kt_right(2));
fprintf('Kt_right-Kb_right = [% .6f, % .6f]\n', ...
    Kt_right(1)-Kb_right(1), ...
    Kt_right(2)-Kb_right(2));
fprintf('Kb left = [% .6f, % .6f]\n', Kb_left(1), Kb_left(2));
fprintf('Kt left = [% .6f, % .6f]\n', Kt_left(1), Kt_left(2));
fprintf('Kt_left-Kb_left = [% .6f, % .6f]\n', ...
    Kt_left(1)-Kb_left(1), ...
    Kt_left(2)-Kb_left(2));
fprintf('det([G1M,G2M]) = %.6e\n', det([G1M,G2M]));

% EFFECTIVE CONTINUUM HAMILTONIAN
% We employ the convention and values from Ref[1]
% note u1 = u, u2 = u^\prime
% with interlayer tunnellings:
u1  = 110; % meV
u2  = 110; % meV               
fai = valley*2*pi/3;
% We inititalise the tunnelling matrix as in Eq(3) of the paper.
% The position dependent factors are irrelevet here, as they dictate 
% how we fill up the full Hamiltonian, with these constant matrises:
T1 = [u1, u2 ;
      u2, u1];
% T_tr 
T2 = [u1,             u2*exp(-1i*fai);
      u2*exp(1i*fai), u1            ];
% T_tl
T3 = [u1,              u2*exp(1i*fai);
      u2*exp(-1i*fai), u1           ];

% Constructing the Hamiltonian
% Number of moiré lattice jumps back and forward to consider
m = generate_moire_lattice_indexing(N_moire);
N_insert = (2*N_moire+1)^2; % Number of matrix insertions                          
N_bands = 4*N_insert;
Energy = zeros(N_bands, k_length);

textprogressbar('>> calculating the band structure: ');
% For each point in K space we want, we need to solve the effective H
for iK = 1 : k_length
   textprogressbar(fix(100*iK/k_length));
   k_now = k_grid(:,iK); 
   U  = zeros(2*N_insert,2*N_insert);
   H1 = zeros(2*N_insert,2*N_insert);
   H2 = zeros(2*N_insert,2*N_insert);
   for i = 1 : N_insert

       % define the momentum 
       k_shift = k_now + m(i,1)*G1M + m(i,2)*G2M;

       q1 = R2 * ( (eye(2)+E1') * (k_shift - Kb_geom) + valley*A1 );
       q2 = R1 * ( (eye(2)+E2') * (k_shift - Kt_geom) + valley*A2 );
       
       % Diagonal components of the blocks
       H1(2*i-1:2*i,2*i-1:2*i) = ...
            alpha * (q1(1) * valley * sx + q1(2) * sy) + Delta1 * sz;
        
       H2(2*i-1:2*i,2*i-1:2*i) = ...
            alpha * (q2(1) * valley * sx + q2(2) * sy) + Delta2 * sz;
       U(2*i-1:2*i,2*i-1:2*i)  = T1;
       for j = 1 : N_insert
            if m(j,1)-m(i,1)==0 && m(j,2)-m(i,2)==-valley
                U(2*i-1:2*i,2*j-1:2*j) = T2; 
            end
            if m(j,1)-m(i,1)==-valley && m(j,2)-m(i,2)==0
                U(2*i-1:2*i,2*j-1:2*j) = T3; 
            end
       end
   end    
   Energy(:,iK) = eig([H1 U' ;
                       U H2]);
end
textprogressbar(' done. plotting...');


set(gcf,'Position',[100 300 400 500])
for i = 1:height(Energy)
    hold on;
    if valley == 1
        plot(Energy(i,:),'-','Linewidth',1.1,'Color','b');
    else
        plot(Energy(i,:),'--','Linewidth',1.1,'Color','r');
    end
    ylabel('E(meV)');
    set(gca,'xticklabel',[])
    box on
end
% axis([0 301 -80 100]);
xticks = [ ...
    1, ...
    N1+1, ...
    N1+N2+1, ...
    N1+N2+N3+1 ];

set(gca,'XTick', xticks);
xticklabels_latex = {'$K_t$', '$K_b$', '$\Gamma$', '$K_t$'};
set(gca, 'XTickLabel', xticklabels_latex);
set(gca,'TickLabelInterpreter','latex')
ylim([-20 20])
xlim([0 k_length])
set(gca,'YTick',-20:10:20);
set(gca,'TickLabelInterpreter','latex')

hold on 
for xline_pos = xticks(2:end-1)
    hold on
    plot([xline_pos xline_pos], [-300 300], '--', 'LineWidth', 0.9, 'Color', [0.8 0.8 0.8]);
end

%% Berry curvature and velocity matrix elements on the BZ grid
u_grid = linspace(-0.5, 0.5 - 1/Nx, Nx);
v_grid = linspace(-0.5, 0.5 - 1/Ny, Ny);
du = u_grid(2)-u_grid(1);
dv = v_grid(2)-v_grid(1);

GammaBZ = Gamma;

BZ = zeros(2,Nx,Ny);
for iu = 1:Nx
    for iv = 1:Ny
        BZ(:,iu,iv) = GammaBZ + u_grid(iu)*G1M + v_grid(iv)*G2M;
    end
end

dk_area = abs(det([G1M, G2M])) * du * dv / (2*pi)^2;

% analytic velocity blocks
dqdx1 = R2 * (eye(2) + E1') * [1;0];
dqdy1 = R2 * (eye(2) + E1') * [0;1];

dqdx2 = R1 * (eye(2) + E2') * [1;0];
dqdy2 = R1 * (eye(2) + E2') * [0;1];

dHx1 = alpha * (valley*dqdx1(1)*sx + dqdx1(2)*sy);
dHy1 = alpha * (valley*dqdy1(1)*sx + dqdy1(2)*sy);

dHx2 = alpha * (valley*dqdx2(1)*sx + dqdx2(2)*sy);
dHy2 = alpha * (valley*dqdy2(1)*sx + dqdy2(2)*sy);

mid1 = 2*N_insert;
mid2 = 2*N_insert + 1;

n_keep_below = 4;
n_keep_above = 4;

band_ids = (mid1 - n_keep_below) : (mid2 + n_keep_above);
N_keep = numel(band_ids);

E_grid     = zeros(Nx,Ny,N_keep,'single');
Omega_grid = zeros(Nx,Ny,N_keep,'single');
vxdiag     = zeros(Nx,Ny,N_keep,'single');
vydiag     = zeros(Nx,Ny,N_keep,'single');

textprogressbar('>> calculating eigensystem, velocities, Berry curvature: ');
for iu = 1:Nx
    textprogressbar(fix(100*iu/Nx));
    for iv = 1:Ny

        k_now = BZ(:,iu,iv);

        U  = zeros(2*N_insert,2*N_insert);
        H1 = zeros(2*N_insert,2*N_insert);
        H2 = zeros(2*N_insert,2*N_insert);

        for i = 1:N_insert
            k_shift = k_now + m(i,1)*G1M + m(i,2)*G2M;

            q1 = R2 * ( (eye(2)+E1') * (k_shift - Kb_geom) + valley*A1 );
            q2 = R1 * ( (eye(2)+E2') * (k_shift - Kt_geom) + valley*A2 );
       
            H1(2*i-1:2*i,2*i-1:2*i) = ...
                alpha * (q1(1)*valley*sx + q1(2)*sy) + Delta1*sz;

            H2(2*i-1:2*i,2*i-1:2*i) = ...
                alpha * (q2(1)*valley*sx + q2(2)*sy) + Delta2*sz;

            U(2*i-1:2*i,2*i-1:2*i) = T1;

            for j = 1 : N_insert

                if m(j,1)-m(i,1)==0 && m(j,2)-m(i,2)==-valley
                    U(2*i-1:2*i,2*j-1:2*j) = T2; 
                end
                
                if m(j,1)-m(i,1)==-valley && m(j,2)-m(i,2)==0
                    U(2*i-1:2*i,2*j-1:2*j) = T3; 
                end
           end
        end

        Hfull = [H1 U'; U H2];

        Vx = zeros(4*N_insert,4*N_insert);
        Vy = zeros(4*N_insert,4*N_insert);
        for i = 1:N_insert
            Vx(2*i-1:2*i,2*i-1:2*i) = dHx1;
            Vy(2*i-1:2*i,2*i-1:2*i) = dHy1;

            j0 = 2*N_insert + 2*i - 1;
            Vx(j0:j0+1,j0:j0+1) = dHx2;
            Vy(j0:j0+1,j0:j0+1) = dHy2;
        end

        [Vec, Eval] = eig(Hfull);
        evals = real(diag(Eval));
        [evals, ord] = sort(evals);
        Vec = Vec(:,ord);

        % rotate velocity matrices to eigenbasis
        Vx_band = Vec' * Vx * Vec;
        Vy_band = Vec' * Vy * Vec;

        % store only what you actually need
        for ia = 1:N_keep
            a = band_ids(ia);
        
            E_grid(iu,iv,ia) = single(evals(a));
            vxdiag(iu,iv,ia) = single(real(Vx_band(a,a)));
            vydiag(iu,iv,ia) = single(real(Vy_band(a,a)));
        
            Om = 0;
            Ea = evals(a);
            for b = 1:N_bands
                if a == b
                    continue
                end
                dE = evals(b) - Ea;
                Om = Om - 2 * imag(Vx_band(a,b) * Vy_band(b,a)) / (dE^2 + 1e-12);
            end
            Omega_grid(iu,iv,ia) = single(Om);
        end

    end
end
textprogressbar(' done.');

mu_grid = linspace(-20,20,141);
%% =========================================================
% Direct BCD from f * d_k Omega  for the middle two bands
% =========================================================

% Fermi function
fFD = @(e,mu) 1 ./ (exp((e-mu)./Temp_mu) + 1);

% locate the two narrow bands inside the kept-band arrays
ia_mid1 = find(band_ids == mid1, 1);
ia_mid2 = find(band_ids == mid2, 1);

mid_keep_ids = (ia_mid1 - N_side) : (ia_mid2 + N_side);
mid_keep_ids = mid_keep_ids(mid_keep_ids >= 1 & mid_keep_ids <= N_keep);
mid_keep_ids = unique(mid_keep_ids);

% arrays vs chemical potential
Dx_arr_fdOmega = zeros(size(mu_grid));
Dy_arr_fdOmega = zeros(size(mu_grid));

% ---------- derivatives: convert (u,v) derivatives -> (kx,ky) derivatives
% k = Gamma + u*G1M + v*G2M
% so [d/du; d/dv] = [G1x G1y; G2x G2y] [d/dkx; d/dky]
Muv_to_k = [G1M(1) G1M(2); ...
            G2M(1) G2M(2)];
Minv = inv(Muv_to_k);

%% =========================================================
% Photocurrent objects:
%   S^{gamma j}       = - int f'(E) v_j v_gamma
%   Omega_z           = int f Omega^{xy}
%   D^i               = - int Omega^{xy} d_i f
%   J^{ijgamma}       = - int f'(E) v_j d_i d_gamma E
%
% Units:
%   k in code is dimensionless, measured in 1/a_g.
%   physical k is in nm^{-1}: k_phys = k_code/a_g.
%
% Therefore:
%   S has units meV and needs no a_g conversion.
%   Omega_z is dimensionless and needs no a_g conversion.
%   D has units nm and gets one factor a_g.
%   J has units meV nm and gets one factor a_g.
% =========================================================

deg_factor_response = spin_factor;

% Use the same band subset as BCD
bands_for_response = mid_keep_ids;

% Arrays versus chemical potential
S_photo = zeros(length(mu_grid), 2, 2);       % S(mu, gamma, j)
D_photo = zeros(length(mu_grid), 2);          % D(mu, i)
J_photo = zeros(length(mu_grid), 2, 2, 2);    % J(mu, i, j, gamma)
Omega_z_photo = zeros(length(mu_grid), 1);    % Omega_z(mu)

% Precompute Cartesian derivatives of energy for each retained band
% dE1{ia,dir}
% dE2{ia,gamma,j}
dE1 = cell(N_keep, 2);
dE2 = cell(N_keep, 2, 2);

textprogressbar('>> precomputing energy derivatives for S and J: ');

for ia = bands_for_response
    textprogressbar(fix(100*(find(bands_for_response==ia,1))/numel(bands_for_response)));

    En = double(E_grid(:,:,ia));

    % First derivatives
    dE1{ia,1} = deriv_k_periodic(En, du, dv, Minv, 1); % d_x E
    dE1{ia,2} = deriv_k_periodic(En, du, dv, Minv, 2); % d_y E

    % Second derivatives
    for gamma = 1:2
        for j = 1:2
            dE2{ia,gamma,j} = deriv_k_periodic(dE1{ia,gamma}, du, dv, Minv, j);
        end
    end

end

textprogressbar(' done.');

%% =========================================================
% Integrate S, Omega_z, D, J versus chemical potential
% =========================================================
dfde  = @(e,mu) -(1./(4*Temp_mu)) .* sech((e-mu)./(2*Temp_mu)).^2;
textprogressbar('>> integrating response objects over mu_grid: ');
for imu = 1:length(mu_grid)
    textprogressbar(fix(100*imu/length(mu_grid)));

    mu = mu_grid(imu);

    S_tmp = zeros(2,2);
    D_tmp = zeros(2,1);
    J_tmp = zeros(2,2,2);
    Om_tmp = 0;

    for ia = bands_for_response
        En = double(E_grid(:,:,ia));
        Om = double(Omega_grid(:,:,ia));

        vx = double(vxdiag(:,:,ia));
        vy = double(vydiag(:,:,ia));

        f0 = fFD(En, mu);
        df_dE = dfde(En, mu);

        % d_i f = df/dE * d_i E
        dfdx = df_dE .* vx;
        dfdy = df_dE .* vy;

        % ---- Omega_z = int f Omega_xy
        Om_tmp = Om_tmp + sum(f0 .* Om, 'all') * dk_area;

        % ---- D^i = - int Omega_xy d_i f
        D_tmp(1) = D_tmp(1) + sum((-Om .* dfdx), 'all') * dk_area;
        D_tmp(2) = D_tmp(2) + sum((-Om .* dfdy), 'all') * dk_area;

        % ---- S^{gamma j} = - int f'(E) v_j v_gamma
        vcell = {vx, vy};
        for gamma = 1:2
            for j_dir = 1:2
                S_tmp(gamma,j_dir) = S_tmp(gamma,j_dir) - ...
                    sum(df_dE .* vcell{j_dir} .* vcell{gamma}, 'all') * dk_area;
            end
        end

        % ---- J^{ijgamma} = - int f'(E) v_j d_i d_gamma E
        for i_dir = 1:2
            for j_dir = 1:2
                for gamma = 1:2
                    J_tmp(i_dir,j_dir,gamma) = J_tmp(i_dir,j_dir,gamma) - ...
                        sum(df_dE .* vcell{j_dir} .* dE2{ia,gamma,i_dir}, 'all') * dk_area;
                end
            end
        end
    end

    % Unit conversions and degeneracy factor
    S_photo(imu,:,:) = deg_factor_response * S_tmp;             % meV
    Omega_z_photo(imu) = deg_factor_response * Om_tmp;          % dimensionless
    D_photo(imu,:) = deg_factor_response * a_g * D_tmp.';       % nm
    J_photo(imu,:,:,:) = deg_factor_response * a_g * J_tmp;     % meV nm
end
textprogressbar(' done.');

%% =========================================================
% BCD consistency check arrays:
%   D = - int Omega * d_k f
%   D =   int f * d_k Omega
% =========================================================

Dx_arr_dfOmega = zeros(size(mu_grid));
Dy_arr_dfOmega = zeros(size(mu_grid));
Dx_arr_fdOmega = zeros(size(mu_grid));
Dy_arr_fdOmega = zeros(size(mu_grid));

% Precompute Berry-curvature derivatives for kept response bands
dOm_dx = cell(N_keep,1);
dOm_dy = cell(N_keep,1);

textprogressbar('>> precomputing Berry-curvature derivatives: ');
for ia = bands_for_response
    textprogressbar(fix(100*(find(bands_for_response==ia,1))/numel(bands_for_response)));

    Om = double(Omega_grid(:,:,ia));

    dOm_dx{ia} = deriv_k_periodic(Om, du, dv, Minv, 1);
    dOm_dy{ia} = deriv_k_periodic(Om, du, dv, Minv, 2);
end
textprogressbar(' done.');

textprogressbar('>> integrating BCD consistency over mu_grid: ');
for imu = 1:length(mu_grid)
    textprogressbar(fix(100*imu/length(mu_grid)));

    mu = mu_grid(imu);

    Dx_df_tmp = 0;
    Dy_df_tmp = 0;
    Dx_fd_tmp = 0;
    Dy_fd_tmp = 0;

    for ia = bands_for_response
        En = double(E_grid(:,:,ia));
        Om = double(Omega_grid(:,:,ia));

        vx = double(vxdiag(:,:,ia));
        vy = double(vydiag(:,:,ia));

        f0 = fFD(En, mu);
        df_dE = dfde(En, mu);

        dfdx = df_dE .* vx;
        dfdy = df_dE .* vy;

        % D = - int Omega * d_k f
        Dx_df_tmp = Dx_df_tmp + sum((-Om .* dfdx), 'all') * dk_area;
        Dy_df_tmp = Dy_df_tmp + sum((-Om .* dfdy), 'all') * dk_area;

        % D = int f * d_k Omega
        Dx_fd_tmp = Dx_fd_tmp + sum(f0 .* dOm_dx{ia}, 'all') * dk_area;
        Dy_fd_tmp = Dy_fd_tmp + sum(f0 .* dOm_dy{ia}, 'all') * dk_area;
    end

    Dx_arr_dfOmega(imu) = Dx_df_tmp;
    Dy_arr_dfOmega(imu) = Dy_df_tmp;
    Dx_arr_fdOmega(imu) = Dx_fd_tmp;
    Dy_arr_fdOmega(imu) = Dy_fd_tmp;
end
textprogressbar(' done.');

% Convert to nm and include degeneracy factor
Dx_arr_dfOmega_nm = deg_factor_response * a_g * Dx_arr_dfOmega;
Dy_arr_dfOmega_nm = deg_factor_response * a_g * Dy_arr_dfOmega;
Dx_arr_fdOmega_nm = deg_factor_response * a_g * Dx_arr_fdOmega;
Dy_arr_fdOmega_nm = deg_factor_response * a_g * Dy_arr_fdOmega;

if isempty(S_total)
    S_total = zeros(size(S_photo));
    Omega_z_total = zeros(size(Omega_z_photo));
    D_total = zeros(size(D_photo));
    J_total = zeros(size(J_photo));
    Dx_arr_dfOmega_nm_total = zeros(size(Dx_arr_dfOmega_nm));
    Dy_arr_dfOmega_nm_total = zeros(size(Dy_arr_dfOmega_nm));
    Dx_arr_fdOmega_nm_total = zeros(size(Dx_arr_fdOmega_nm));
    Dy_arr_fdOmega_nm_total = zeros(size(Dy_arr_fdOmega_nm));
end

S_total = S_total + S_photo;
Omega_z_total = Omega_z_total + Omega_z_photo;
D_total = D_total + D_photo;
J_total = J_total + J_photo;
Dx_arr_dfOmega_nm_total = Dx_arr_dfOmega_nm_total + Dx_arr_dfOmega_nm;
Dy_arr_dfOmega_nm_total = Dy_arr_dfOmega_nm_total + Dy_arr_dfOmega_nm;
Dx_arr_fdOmega_nm_total = Dx_arr_fdOmega_nm_total + Dx_arr_fdOmega_nm;
Dy_arr_fdOmega_nm_total = Dy_arr_fdOmega_nm_total + Dy_arr_fdOmega_nm;

if valley == +1
    valley_key = 'valley_plus';
else
    valley_key = 'valley_minus';
end

valley_results.(valley_key).S_photo = S_photo;
valley_results.(valley_key).Omega_z_photo = Omega_z_photo;
valley_results.(valley_key).D_photo = D_photo;
valley_results.(valley_key).J_photo = J_photo;
valley_results.(valley_key).Dx_arr_dfOmega_nm = Dx_arr_dfOmega_nm;
valley_results.(valley_key).Dy_arr_dfOmega_nm = Dy_arr_dfOmega_nm;
valley_results.(valley_key).Dx_arr_fdOmega_nm = Dx_arr_fdOmega_nm;
valley_results.(valley_key).Dy_arr_fdOmega_nm = Dy_arr_fdOmega_nm;

end

S_photo = S_total;
Omega_z_photo = Omega_z_total;
D_photo = D_total;
J_photo = J_total;
Dx_arr_dfOmega_nm = Dx_arr_dfOmega_nm_total;
Dy_arr_dfOmega_nm = Dy_arr_dfOmega_nm_total;
Dx_arr_fdOmega_nm = Dx_arr_fdOmega_nm_total;
Dy_arr_fdOmega_nm = Dy_arr_fdOmega_nm_total;

fprintf('\n===== TBG valley cancellation diagnostics =====\n');
fprintf('max abs Omega_z for xi=+1: %.6e\n', max(abs(valley_results.valley_plus.Omega_z_photo)));
fprintf('max abs Omega_z for xi=-1: %.6e\n', max(abs(valley_results.valley_minus.Omega_z_photo)));
fprintf('max abs Omega_z for valley sum: %.6e\n', max(abs(Omega_z_photo)));
fprintf('max abs D for xi=+1: %.6e\n', max(sqrt(sum(valley_results.valley_plus.D_photo.^2, 2))));
fprintf('max abs D for xi=-1: %.6e\n', max(sqrt(sum(valley_results.valley_minus.D_photo.^2, 2))));
fprintf('max abs D for valley sum: %.6e\n', max(sqrt(sum(D_photo.^2, 2))));

%% =========================================================
% Save all four key quantities
% =========================================================

response_data = struct();

response_data.material_name = 'twisted_bilayer_graphene';
response_data.model_name = 'fundamental_continuum';
response_data.version = 'v1_response_objects';

response_data.mu_grid = mu_grid;
response_data.units.energy = 'meV';
response_data.units.length = 'nm';
response_data.units.momentum = 'nm^{-1}';
response_data.units.S = 'meV';
response_data.units.Omega_z = 'dimensionless';
response_data.units.D = 'nm';
response_data.units.J = 'meV nm';

response_data.S_photo = S_photo;             % S(mu, gamma, j)
response_data.Omega_z_photo = Omega_z_photo; % Omega_z(mu)
response_data.D_photo = D_photo;             % D(mu, i)
response_data.J_photo = J_photo;             % J(mu, i, j, gamma)

response_data.index_convention = ...
    '1=x, 2=y. S(mu,gamma,j), D(mu,i), J(mu,i,j,gamma).';
response_data.definitions.Omega_z = 'Omega_z(mu) = sum_n int_k f_n(k,mu) Omega_n^{xy}(k)';
response_data.definitions.D = 'D^i(mu) = - sum_n int_k Omega_n^{xy}(k) partial_i f_n(k,mu)';
response_data.definitions.S = 'S^{gamma j}(mu) = - sum_n int_k f''_n(E) v_n^j v_n^gamma';
response_data.definitions.J = 'J^{i j gamma}(mu) = - sum_n int_k f''_n(E) v_n^j partial_i partial_gamma E_n';

response_data.parameters.theta = theta;
response_data.parameters.valley_list = valley_list;
response_data.parameters.spin_factor = spin_factor;
response_data.parameters.valley_summed_explicitly = true;
response_data.parameters.valley_note = ...
    'TBG response explicitly sums xi=+1 and xi=-1; degeneracy factor contains spin only.';
response_data.parameters.alpha = alpha;
response_data.parameters.a_g_nm = a_g;
response_data.parameters.eps_strain = eps_strain;
response_data.parameters.phi_strain = phi_strain;
response_data.parameters.nu = nu;
response_data.parameters.beta_gr = beta_gr;
response_data.parameters.Delta1 = Delta1;
response_data.parameters.Delta2 = Delta2;
response_data.parameters.Temp_mu = Temp_mu;
response_data.parameters.N_moire = N_moire;
response_data.parameters.Nx = Nx;
response_data.parameters.Ny = Ny;
response_data.parameters.dk_area = dk_area;
response_data.parameters.deg_factor_response = spin_factor;
response_data.parameters.band_ids = band_ids;
response_data.parameters.bands_for_response = bands_for_response;

response_data.grid.coordinate = 'moire_fractional_uv';
response_data.grid.saved_momentum_unit = 'nm^{-1}';
response_data.grid.Nx = Nx;
response_data.grid.Ny = Ny;
response_data.grid.du = du;
response_data.grid.dv = dv;
response_data.grid.G1M = G1M;
response_data.grid.G2M = G2M;
response_data.grid.dk_area = dk_area;

response_data.diagnostics.Omega_z_xi_plus = valley_results.valley_plus.Omega_z_photo;
response_data.diagnostics.Omega_z_xi_minus = valley_results.valley_minus.Omega_z_photo;
response_data.diagnostics.D_xi_plus = valley_results.valley_plus.D_photo;
response_data.diagnostics.D_xi_minus = valley_results.valley_minus.D_photo;
response_data.Dx_arr_dfOmega_nm = Dx_arr_dfOmega_nm;
response_data.Dy_arr_dfOmega_nm = Dy_arr_dfOmega_nm;
response_data.Dx_arr_fdOmega_nm = Dx_arr_fdOmega_nm;
response_data.Dy_arr_fdOmega_nm = Dy_arr_fdOmega_nm;

output_name = fullfile(fileparts(mfilename('fullpath')), 'tbg_photocurrent_objects_valleysum.mat');
validate_and_save_response_objects(output_name, response_data);

disp('Saved valley-summed S_photo, Omega_z_photo, D_photo, J_photo to tbg_photocurrent_objects_valleysum.mat');

%% =========================================================
% Final diagnostic plots
% 1) D_x, D_y consistency check
% 2) S on left axis, J on right axis
% =========================================================

% ---------------------------------------------------------
% 1) D_x, D_y consistency check
% ---------------------------------------------------------
figure('Position',[100 100 650 420]);

plot(mu_grid, Dx_arr_dfOmega_nm, 'r-',  'LineWidth', 2); hold on
plot(mu_grid, Dx_arr_fdOmega_nm, 'r--', 'LineWidth', 2);
plot(mu_grid, Dy_arr_dfOmega_nm, 'b-',  'LineWidth', 2);
plot(mu_grid, Dy_arr_fdOmega_nm, 'b--', 'LineWidth', 2);

xlabel('\mu [meV]')
ylabel('D [nm]')
legend( ...
    'D_x: -\Omega \partial_k f', ...
    'D_x: f \partial_k \Omega', ...
    'D_y: -\Omega \partial_k f', ...
    'D_y: f \partial_k \Omega', ...
    'Location','best')
title('BCD consistency check')
box on
grid on

% ---------------------------------------------------------
% 2) Compact diagnostic: |D|, rho_xx / Sxx, and J_xxx
% ---------------------------------------------------------

% Magnitude of BCD from the two equivalent formulas
Dmag_dfOmega = sqrt(Dx_arr_dfOmega_nm.^2 + Dy_arr_dfOmega_nm.^2);
Dmag_fdOmega = sqrt(Dx_arr_fdOmega_nm.^2 + Dy_arr_fdOmega_nm.^2);

% Drude/stiffness-like object
% In your current notation, S_photo(:,1,1) is the FS/intraband S_xx object
rho_xx = squeeze(S_photo(:,1,1));

% Jerk component
Jxxx = squeeze(J_photo(:,1,1,1));

figure('Position',[150 120 780 460]);

yyaxis left
plot(mu_grid, Dmag_dfOmega, 'k-',  'LineWidth', 2); hold on
plot(mu_grid, Dmag_fdOmega, 'k--', 'LineWidth', 2);
ylabel('$|\mathbf D|$ [nm]', 'Interpreter','latex')

yyaxis right
plot(mu_grid, rho_xx, 'b-', 'LineWidth', 2); hold on

ylabel('$\rho_{xx}$ [meV]', 'Interpreter','latex')

xlabel('$\mu$ [meV]', 'Interpreter','latex')
title('BCD, stiffness, and jerk diagnostics', 'Interpreter','latex')

legend( ...
    '$|\mathbf D|: -\Omega\,\partial_k f$', ...
    '$|\mathbf D|: f\,\partial_k\Omega$', ...
    '$\rho_{xx}$', ...
    'Interpreter','latex', ...
    'Location','best')

box on
grid on

function dFdk = deriv_k_periodic(F, du, dv, Minv, dir)
    % Periodic central derivative in Cartesian k.
    %
    % k = Gamma + u G1M + v G2M
    %
    % [d/du; d/dv] = [G1x G1y; G2x G2y] [d/dkx; d/dky]
    % so
    % [d/dkx; d/dky] = Minv [d/du; d/dv].
    %
    % dir = 1 gives d/dkx
    % dir = 2 gives d/dky

    dFdu = (circshift(F,[-1,0]) - circshift(F,[1,0])) / (2*du);
    dFdv = (circshift(F,[0,-1]) - circshift(F,[0,1])) / (2*dv);

    dFdk = Minv(dir,1)*dFdu + Minv(dir,2)*dFdv;
end

%% Generate moiré reciprocal lattice points indexing

function lattice_points = generate_moire_lattice_indexing(N)
    lattice_points = zeros((2*N + 1)^2, 2);
    index = 1;
    for i = -N:N
        for j = -N:N
            lattice_points(index, :) = [i, j];
            index = index + 1;
        end
    end
end

%% Linspace between two vectors that doesn't double include points

function path = path_add(A,p1,p2,N)
    if isempty(A)
        path = [linspace(p1(1),p2(1),N+1); linspace(p1(2),p2(2),N+1)];
        return 
    end
    x = linspace(p1(1),p2(1),N+1);
    y = linspace(p1(2),p2(2),N+1);
    path = [A(1,:) x(2:end); A(2,:) y(2:end)];
end

% Just a fancy progressbar. Delete/cpmment it if complicated.

% function textprogressbar(c)
%     % Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
%     % Version: 1.0
%     % Changes tracker:  29.06.2010  - First version
%     %% Initialization
%     persistent strCR;          % Carriage return pesistent variable
%     
%     % Vizualization parameters
%     strPercentageLength = 10;  % Length of percentage string (must be >5)
%     strDotsMaximum      = 10;  % The total number of dots in a progress bar
%     
%     %% Main 
%     
%     if isempty(strCR) && ~ischar(c)
%         % Progress bar must be initialized with a string
%         error('The text progress must be initialized with a string');
%     elseif isempty(strCR) && ischar(c)
%         % Progress bar - initialization
%         fprintf('%s',c);
%         strCR = -1;
%     elseif ~isempty(strCR) && ischar(c)
%         % Progress bar  - termination
%         strCR = [];  
%         fprintf([c '\n']);
%     elseif isnumeric(c)
%         % Progress bar - normal progress
%         c = floor(c);
%         percentageOut = [num2str(c) '%%'];
%         percentageOut = [percentageOut repmat(' ',1,...
%                              strPercentageLength-length(percentageOut)-1)];
%         nDots = floor(c/100*strDotsMaximum);
%         dotOut = ['[' repmat('.',1,nDots) ...
%                                    repmat(' ',1,strDotsMaximum-nDots) ']'];
%         strOut = [percentageOut dotOut];
%         
%         % Print it on the screen
%         if strCR == -1
%             % Don't do carriage return during first run
%             fprintf(strOut);
%         else
%             % Do it during all the other runs
%             fprintf([strCR strOut]);
%         end
%         
%         % Update carriage return
%         strCR = repmat('\b',1,length(strOut)-1);
%         
%     else
%         % Any other unexpected input
%         error('Unsupported argument type');
%     end
% end
