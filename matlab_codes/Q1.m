%% TP Traitement d'antenne 2D – Diagramme de réception (RX) — VERSION CORRIGÉE
clear; clc; close all;

%% Paramètres généraux
matFile = 'ant_IMAGING_RX.mat';
if ~exist(matFile,'file'); error('Fichier introuvable : %s', matFile); end

c  = 3e8;                   % m/s
fc = 77e9;                  % Hz
lambda = c/fc;

% Si X,Y,Z sont déjà en longueurs d'onde, passe à true
positions_in_wavelength = false;

% Grille d'angles (en degrés pour l'IHM)
az_span = -90:0.5:90;       % "circulaire" selon la convention de ton encadrant
el_span = -60:0.5:60;       % "élévation"
[AZ, EL] = meshgrid(az_span, el_span);

% Pointage (boresight)
az0_deg = 0;
el0_deg = 0;
dB_floor = -60;

%% Chargement positions (on n'utilise QUE X,Y,Z)
S = load(matFile);
if ~isfield(S,'ant_IMAGING_RX'), error('Variable ant_IMAGING_RX absente.'); end
A = S.ant_IMAGING_RX;
if size(A,2) < 3, error('ant_IMAGING_RX doit contenir au moins X,Y,Z.'); end

pos_m = A(:,1:3);           % positions en mètres (cf. capture)
N = size(pos_m,1);
fprintf('Réseau chargé : %d capteurs RX.\n', N);

% Si les positions étaient fournies en λ, on les reconvertit en mètres pour spatialSteeringVector
if positions_in_wavelength
    pos_m = pos_m * lambda; % -> mètres
    fprintf('Positions initialement en λ, reconverties en mètres (λ = %.6f m).\n', lambda);
else
    fprintf('Positions déjà en mètres (λ = %.6f m).\n', lambda);
end

%% Formation de voies (steering + pondération uniforme) avec la CONVENTION DEMANDÉE
% Attention : spatialSteeringVector attend des angles en RADIANS et des positions EN MÈTRES.
a0 = spatialSteeringVector(pos_m, lambda, deg2rad(az0_deg), deg2rad(el0_deg)); % N×1
w  = a0 / N;                               % poids uniformes + steering

% Carte 2D du diagramme |w^H a(az,el)|^2 (normalisée)
AF_2D = array_factor_2D_convention(pos_m, lambda, w, AZ, EL); % puissance
AF_2D = AF_2D / (max(AF_2D(:)) + eps);

%% Tracés
figure('Name','Diagramme 2D (AZ×EL)');
imagesc(az_span, el_span, to_dB(AF_2D));
axis xy; cb = colorbar; caxis([dB_floor 0]);
xlabel('Azimut (deg)'); ylabel('Élévation (deg)'); ylabel(cb,'Gain (dB)');
title(sprintf('Diagramme RX (az=%g°, el=%g°)', az0_deg, el0_deg));

[~, idx_el0] = min(abs(el_span - 0));
AF_az_cut = AF_2D(idx_el0, :);
figure('Name','Coupe azimut (el=0°)');
plot(az_span, to_dB(AF_az_cut), 'LineWidth',1.5);
ylim([dB_floor 0]); grid on;
xlabel('Azimut (deg)'); ylabel('Gain (dB)');

[~, idx_az0] = min(abs(az_span - 0));
AF_el_cut = AF_2D(:, idx_az0);
figure('Name','Coupe élévation (az=0°)');
plot(el_span, to_dB(AF_el_cut), 'LineWidth',1.5);
ylim([dB_floor 0]); grid on;
xlabel('Élévation (deg)'); ylabel('Gain (dB)');

%% ======================== SOUS-FONCTIONS ========================

% -- Convention exigée par l'encadrant --
% antenna : Nx3 positions EN METRES, lambda (m), circulaire (rad), elevation (rad)
% Repère : antenne dans XOY, OX à gauche, OY en haut, OZ vers l'avant (boresight).
% Vecteur direction u = [sin(circ)*cos(el);  sin(el);  cos(circ)*cos(el)]
function spVect = spatialSteeringVector(antenna, lambda, circulaire, elevation)
    vect = [ sin(circulaire)*cos(elevation); ...
             sin(elevation); ...
             cos(circulaire)*cos(elevation) ];      % 3×1
    dpl = (2*pi/lambda);
    x   = (antenna(:,1:3) * vect);                  % N×1, en mètres·(sans unité)
    spVect = exp(1i * dpl * x);                     % N×1
end

% Array factor 2D avec la même convention d'angles (AZ,EL en DEGRÉS)
% Renvoie la puissance |w^H a|^2
function AF = array_factor_2D_convention(pos_m, lambda, w, AZ_deg, EL_deg)
    [nr, nc] = size(AZ_deg);
    AF = zeros(nr, nc);
    for i = 1:nr
        for j = 1:nc
            a_ij = spatialSteeringVector(pos_m, lambda, deg2rad(AZ_deg(i,j)), deg2rad(EL_deg(i,j)));
            AF(i,j) = abs(w' * a_ij).^2;
        end
    end
end

function y = to_dB(x), y = 10*log10(abs(x) + eps); end
