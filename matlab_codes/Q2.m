%% TP Traitement d'antenne 2D – Diagramme de réception (RX)
clear; clc; close all;

%% Paramètres généraux
matFile = 'ant_IMAGING_RX.mat';
load(matFile)
if ~exist(matFile,'file')
    error('Fichier introuvable : %s', matFile);
end

c  = 3e8;        % m/s
fc = 77e9;       % Hz (adapter si besoin)
lambda = c/fc;

positions_in_wavelength = false;   % true si X,Y,Z sont déjà en λ

az_span = -90:0.5:90;
el_span = -60:0.5:60;
[AZ, EL] = meshgrid(az_span, el_span);

az0_deg = 0;
el0_deg = 0;
dB_floor = -60;

%% essai
%dimension de l'antenne élémentaire  Rx en vertical
HRX=0.01;
%dimension de l'antenne élémentaire  Rx en horizontal
LRX=0.0039;


angle_brouilleur = 7.5;
angle_pointage_azimut = 0;
s0 = spatialSteeringVector(ant_IMAGING_RX, lambda, deg2rad(angle_pointage_azimut), 0);
sj = spatialSteeringVector(ant_IMAGING_RX, lambda, deg2rad(angle_brouilleur), 0);
R = 200*(sj*sj') +  eye(16);
w_mvdr = (R \ s0) / (s0' * (R \ s0));

for iboucle=1:length(az_span)
spVect = spatialSteeringVector(ant_IMAGING_RX, lambda, deg2rad(az_span(iboucle)), 0);
diag(iboucle) = spVect'*s0/(spVect'*spVect);
% diag elementaire
urxaz=pi*LRX*sin(deg2rad(az_span(iboucle)))/lambda;
diagrxaz=20*log10(abs(sinc(urxaz/pi))+0.1)+14.9-0.0025*abs((deg2rad(az_span(iboucle))*180/pi).^2);
diagrxaz=diagrxaz-15.727;%nomalisation du gain à 0 dB
diagrxaz=(10).^(diagrxaz/20);%champ électrique
diagc(iboucle) = diagrxaz * diag(iboucle);

diag_abf(iboucle)= abs(w_mvdr' * spVect);
diag_abf_c(iboucle)= diagrxaz * diag_abf(iboucle);
end
figure, plot(az_span,10*log10(abs(diag))), hold on
plot(az_span,10*log10(abs(diag_abf)))
xline(angle_brouilleur, '--m', 'LineWidth', 2); % vertical dashed red line at x=3
xlabel('Azimut (degres)');
ylabel('Gain normalise');
legend('Diagramme DBF','Diagramme ABF','position brouilleur');

figure, plot(az_span,10*log10(abs(diagc))), hold on
plot(az_span,10*log10(abs(diag_abf_c)))
xline(angle_brouilleur, '--m', 'LineWidth', 2); % vertical dashed red line at x=3
title('Gain relatif')
xlabel('Azimut (degres)');
ylabel('Gain normalise');
legend('Diagramme DBF avec gain elementaire','Diagramme ABF avec gain elementaire','position brouilleur');

figure, plot(az_span,10*log10(abs(diagc)*16*16)), hold on
plot(az_span,10*log10(abs(diag_abf_c)*16*16))
xline(angle_brouilleur, '--m', 'LineWidth', 2); % vertical dashed red line at x=3
title('Gain absolu')
xlabel('Azimut (degres)');
ylabel('Gain normalise');
legend('Diagramme DBF avec gain elementaire','Diagramme ABF avec gain elementaire','position brouilleur');


%% Chargement et extraction (on n'utilise QUE X,Y,Z)
S = load(matFile);
if ~isfield(S,'ant_IMAGING_RX'), error('Variable ant_IMAGING_RX absente.'); end
A = S.ant_IMAGING_RX;
if size(A,2) < 6, error('ant_IMAGING_RX doit faire Nx6.'); end

pos = A(:,1:3);              % [X Y Z]
N   = size(pos,1);
disp(sprintf('Réseau chargé : %d capteurs RX.', N));

% Conversion en longueurs d'onde (MATLAB n'a pas l'opérateur ?:)
if positions_in_wavelength
    pos_lambda = pos;
    disp('Positions déjà en λ (pas de conversion).');
else
    pos_lambda = pos / lambda;
    disp(sprintf('Positions en mètres converties en λ (lambda = %.4g m).', lambda));
end

%% Formation de voies (SANS pondérations RX)
vec_guidage0 = steering_vec(pos_lambda, az0_deg, el0_deg);   % steering vers (0°,0°)
w  = vec_guidage0 / N;                                       % poids uniformes + steering

% Carte 2D du diagramme |w^H a(az,el)|^2 (normalisée)
AF_2D = array_factor_2D(pos_lambda, w, AZ, EL);
AF_2D = AF_2D / max(AF_2D(:) + eps);

%% Tracés
figure('Name','Diagramme 2D (AZ×EL)');
imagesc(az_span, el_span, to_dB(AF_2D));
axis xy; colorbar; caxis([dB_floor 0]);
xlabel('Azimut (deg)'); ylabel('Élévation (deg)');
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

%% Sous-fonctions
function vec_guidage = steering_vec(pos_lambda, az_deg, el_deg)
    az = deg2rad(az_deg); el = deg2rad(el_deg);
    u  = [cos(el)*cos(az), cos(el)*sin(az), sin(el)];
    phase = 2*pi * (pos_lambda * u.');
    vec_guidage = exp(1j * phase);
end

function AF = array_factor_2D(pos_lambda, w, AZ, EL)
    [nr, nc] = size(AZ); AF = zeros(nr,nc);
    for i = 1:nr
        for j = 1:nc
            vec_guidage = steering_vec(pos_lambda, AZ(i,j), EL(i,j));
            AF(i,j) = abs(w' * vec_guidage).^2;
        end
    end
end

function y = to_dB(x)
    y = 10*log10(abs(x) + eps);
end
