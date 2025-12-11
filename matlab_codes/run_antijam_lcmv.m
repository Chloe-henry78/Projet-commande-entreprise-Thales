function run_antijam_lcmv()
% --- Anti-brouillage LCMV avec convention encadrant (spatialSteeringVector) ---
% u = [sin(az)*cos(el); sin(el); cos(az)*cos(el)], az/el en radians
% a(az,el) = exp( j*(2*pi/lambda)*POS*u ), POS en mètres

% ===== Réglages =====
matFile = 'ant_IMAGING_RX.mat';
fc  = 77e9;                       % Hz
c   = 3e8;                        % m/s
lambda = c/fc;
positions_in_wavelength = false;  % true si POS est déjà en λ (sinon POS en mètres)

% angles (degrés) – "az" = circulaire, "el" = élévation
az_des = 0;   el_des = 0;         % signal d'intérêt
az_jam = -25; el_jam = 0;         % brouilleur

% ===== Chargement géométrie + snapshots (si dispo) =====
[POS_m, Xsnap] = load_geometry_and_snapshots_meters(matFile, lambda, positions_in_wavelength);

% ===== Vecteurs directionnels (utile + jam) — convention encadrant =====
aU = steer_enc(POS_m, lambda, az_des, el_des);   % Mx1
aJ = steer_enc(POS_m, lambda, az_jam, el_jam);   % Mx1
C  = [aU, aJ];                                   % contraintes
g  = [1; 0];                                     % gain unité utile, nul jam

% ===== Estimation de R (ou simulation) =====
if isempty(Xsnap)
    fprintf('Aucun snapshot trouvé : simulation pour démo…\n');
    K = 4000; SNRdB = 0; JNRdB = 20; noisevar = 1;
    Xsnap = simulate_snapshots_enc(POS_m, lambda, az_des, el_des, az_jam, el_jam, SNRdB, JNRdB, noisevar, K);
end

Rx = (Xsnap*Xsnap')/size(Xsnap,2);           % SMI
alpha = 1e-3 * trace(Rx)/size(Rx,1);         % diagonal loading
Rx = Rx + alpha*eye(size(Rx));


% ====== DIAGNOSTIC aliasage / faisabilité des contraintes ======
% Corrélation géométrique entre l'utile et le jam (1 = colinéaires)
rho = abs(aU' * aJ) / (norm(aU)*norm(aJ));
fprintf('corr(aU,aJ) = %.4f\n', rho);

% Conditionnement de (C^H R^{-1} C)
G = C' * (Rx \ C);
kappa = cond(G);
fprintf('cond(C^H R^{-1} C) = %.2e\n', kappa);

if rho > 0.98 || kappa > 1e6
    warning('Contraintes quasi colinéaires / mal conditionnées : risque de scission du lobe principal.');
end
% ================================================================

% ===== Poids LCMV =====
w = (Rx\C) * ((C'*(Rx\C))\g);                % pas d'inv() explicite

w = w / (aU' * w);   % impose numériquement w^H aU = 1

% ===== Évaluation du diagramme sur grille =====
az_span = -90:0.5:90;
el_span = -60:0.5:60;
[AZ, EL] = meshgrid(az_span, el_span);

PAT = pattern_grid_enc(POS_m, w, lambda, AZ, EL);   % |w^H a|^2
PATdB = 10*log10(PAT / max(PAT(:)));
dB_floor = -60;  PATdB = max(PATdB, dB_floor);

% ===== Tracés =====
figure('Name','LCMV 2D (conv. encadrant)');
imagesc(az_span, el_span, PATdB); axis xy; colormap parula; colorbar;
xlabel('Azimut (°)'); ylabel('Élévation (°)');
title('Diagramme LCMV (dB) — convention spatialSteeringVector');
hold on;
plot(az_des, el_des,'wo','MarkerFaceColor','w','MarkerSize',6);
plot(az_jam, el_jam,'wx','LineWidth',1.2,'MarkerSize',10);
hold off;

% Coupe 1D à l'élévation de l'utile
[~, idx_el] = min(abs(el_span - el_des));
figure('Name','Coupe 1D (EL=utile)');
plot(az_span, PATdB(idx_el,:),'LineWidth',1.6); grid on;
xlabel('Azimut (°)'); ylabel('Gain (dB)');
title(sprintf('Coupe @ élévation = %g°', el_des));
ylim([dB_floor 0]);
xline(az_des,'--','Utile'); xline(az_jam,':','Jam');

% Affiche le niveau au jam
[~, idx_az_jam] = min(abs(az_span - az_jam));
gain_jam_dB = PATdB(idx_el, idx_az_jam);
fprintf('Gain au jam ~ %.1f dB (doit être << 0 dB)\n', gain_jam_dB);
end

% ===================== FONCTIONS AIDE (même fichier) =====================

function [POS_m, Xsnap] = load_geometry_and_snapshots_meters(matFile, lambda, in_lambda)
% Renvoie POS en mètres (convertit si positions en λ)
POS_m = []; Xsnap = [];
S = load(matFile);

% Géométrie
if isfield(S,'ant_IMAGING_RX')
    A = S.ant_IMAGING_RX;
    if isnumeric(A) && size(A,2) >= 3
        P = double(A(:,1:3));
    elseif isstruct(A) && all(isfield(A,{'X','Y','Z'}))
        P = [A.X(:), A.Y(:), A.Z(:)];
    else
        error('Format ant_IMAGING_RX non reconnu.');
    end
else
    error('ant_IMAGING_RX absent dans %s.', matFile);
end

if in_lambda
    POS_m = P * lambda;   % λ -> mètres
else
    POS_m = P;            % déjà en mètres
end
fprintf('Geometry loaded: %d elements. Positions en mètres.\n', size(POS_m,1));

% Snapshots (optionnels)
cands = {'Xsnap','X','SNAPSHOTS','snap','data','rx'};
for k=1:numel(cands)
    if isfield(S, cands{k}), Xsnap = double(S.(cands{k})); break; end
end
if ~isempty(Xsnap) && size(Xsnap,1) > size(Xsnap,2)
    Xsnap = Xsnap.'; % MxK
end
if isempty(Xsnap), fprintf('No snapshots found — simulation will be used.\n'); end
end

function a = steer_enc(POS_m, lambda, az_deg, el_deg)
% Steering (convention encadrant) pour 1 angle ou un vecteur d'angles.
% Az/El en degrés (scalaires ou vecteurs). POS en mètres. Sortie: MxN.
az = deg2rad(az_deg(:));
el = deg2rad(el_deg(:));
u  = [ sin(az).*cos(el), ...   % x
       sin(el),           ...  % y
       cos(az).*cos(el)   ...  % z
     ];                         % N×3
phase = (2*pi/lambda) * (POS_m * u.');   % M×N   (+j, même que spatialSteeringVector)
a = exp(1j*phase);                        % M×N
end

function PAT = pattern_grid_enc(POS_m, w, lambda, AZ, EL)
% |w^H a(az,el)|^2 sur grille (AZ, EL en degrés)
[H,W] = size(AZ);
PAT = zeros(H,W);
for i = 1:H
    a_row = steer_enc(POS_m, lambda, AZ(i,:), EL(i,:));  % M×W
    PAT(i,:) = abs(w' * a_row).^2;                       % 1×W
end
end

function X = simulate_snapshots_enc(POS_m, lambda, azU, elU, azJ, elJ, SNRdB, JNRdB, noisevar, K)
% Simulation cohérente avec la même convention de steering
M = size(POS_m,1);
aU = steer_enc(POS_m, lambda, azU, elU);  % M×1
aJ = steer_enc(POS_m, lambda, azJ, elJ);  % M×1
s  = (randn(1,K)+1j*randn(1,K))/sqrt(2);  % utile
j  = (randn(1,K)+1j*randn(1,K))/sqrt(2);  % jam
n  = sqrt(noisevar/2)*(randn(M,K)+1j*randn(M,K));
Ps = 10^(SNRdB/10) * noisevar;
Pj = 10^(JNRdB/10) * noisevar;
s = s * sqrt(Ps); j = j * sqrt(Pj);
X = aU*s + aJ*j + n;
end

