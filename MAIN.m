clear all;
close all;

addpath(genpath('.'));

% -------------------------------------------------------------------------
% CONSTANTES
% -------------------------------------------------------------------------
CST.RADAR_IMAGING = 1;

% -------------------------------------------------------------------------
% Set data file
% -------------------------------------------------------------------------

radar_type = CST.RADAR_IMAGING;  % Radar type : IMAGING
repData='C:\Users\chloe\OneDrive\Bureau\IMTA 2A\Commande entreprise\Correction_steer_abf (1)\Correction_steer_abf';
filtreFic = [{'2025_11_10_21_5_47_FO6_30M_45deg.data'}]; % Select only somes files in the folder


% Look for the files verifying the selection criteria (filtreFic)
% filtreFic = [];
[listeFic] = rechercheListeFichiers(repData, '*.data', filtreFic);

if ~isempty(listeFic), fprintf(1, "Liste des fichiers a traiter :\n");fprintf(1, "\t%s\n", listeFic{:});
else,                  fprintf(1, "Erreur, aucun fichier trouve\n");
end

% -------------------------------------------------------------------------
% FILE PARAMETERS
% -------------------------------------------------------------------------
Debut                   = 1;  % Start index
Nb_Data_Dec             = 24; % Number of bundles to read (from index "debut");
% Nb_Data_Dec=-1 => will read all bundles

% Display configuration
% ------------------------------------------------------------------
visu.antenne            = -1;   % Visu de la configuration antennaire
visu.paramFO            = -2;   % Display waveform paramters

% Video camera (alone)
% ------------------------------------------------------------------
visu.camera             = 50;    % Display camera

% Range Doppler display
% ------------------------------------------------------------------
visu.rangeDoppler             = 3;    % Display camera
visu.zeroDoppler              = 4;

% Global display
% ------------------------------------------------------------------
visu.global             = 5;    % Display camera
visu.sauveVideoCarteXYZ = false; % to save image in movie

% Visu SNR
visu.num_cd_cible = 65;

%%

% -------------------------------------------------------------------------
% PROCESS OF ALL SELECTED FILES
% -------------------------------------------------------------------------
for iFic=1:length(listeFic)
    [DBF, paramFO, Analyse] = Processing(CST, repData, listeFic{iFic}, radar_type, Debut, Nb_Data_Dec, visu);
end

%% Faisceau forme
% suivi du bruit

Pbruit_faisceaux = zeros(Nb_Data_Dec,paramFO.nb_dir_tx);

for iboucle=1:Nb_Data_Dec
    for jboucle=1:paramFO.nb_dir_tx %nb faisceaus Tx
        Pbruit_faisceaux(iboucle,jboucle) = Analyse(iboucle).bruit(jboucle);
    end
end
figure,
hold on,
axis xy;

if paramFO.nb_dir_tx > 1
    for jboucle=1:paramFO.nb_dir_tx %nb faisceaus Tx
        subplot(4,4,jboucle)
        hold on
        plot(Pbruit_faisceaux(:,jboucle),'bx-');
        %     xlabel('Numero rafale')
        %     ylabel('Puissance bruit')
        axis([0 Nb_Data_Dec+1 60 90])
        title(['azimut : ' num2str(rad2deg(paramFO.tab_dir_ffc_tx(iboucle))) ' °'])

    end
    sgtitle('Niveau de bruit')
    

else
    plot(Pbruit_faisceaux,'bx-')
    axis([0 Nb_Data_Dec+1 60 90])
    xlabel('Numero rafale')
    ylabel('Puissance bruit')
    title('Puissance du bruit')
end

% suivi SNR point fixe

SNR_faisceaux = zeros(Nb_Data_Dec,paramFO.nb_dir_tx);

for iboucle=1:Nb_Data_Dec
    for jboucle=1:paramFO.nb_dir_tx %nb faisceaus Tx
        SNR_faisceaux(iboucle,jboucle) = Analyse(iboucle).SNR_cible(jboucle);
    end
end

figure,
hold on,
title('SNR point fixe')
if paramFO.nb_dir_tx > 1
    for jboucle=1:paramFO.nb_dir_tx %nb faisceaus Tx
        subplot(4,4,jboucle)
        hold on
        grid on
        plot(SNR_faisceaux(:,jboucle),'kx-')
        axis([0 Nb_Data_Dec+1 0 100])
        title(['azimut : ' num2str(rad2deg(paramFO.tab_dir_ffc_tx(jboucle,1))) ' °'])

    end
    sgtitle('SNR point fixe')
else
    plot(SNR_faisceaux,'kx-')
    axis([0 Nb_Data_Dec+1 60 90])
    xlabel('Numero rafale')
    ylabel('SNR (dB)')
    title('SNR cible')
end


% SNR par pointage
if paramFO.nb_dir_tx > 1
    figure, hold on
    for iboucle=1:Nb_Data_Dec
        plot(rad2deg(squeeze(paramFO.tab_dir_ffc_tx(:,1))),SNR_faisceaux(iboucle,:))
        xlabel('Azimut (deg)')
        ylabel('SNR (dB)')
        title('SNR point fixe');
    end
end




%% par voies RX
if 1
    % niveau de bruit par RX
    Pbruit_RX = zeros(Nb_Data_Dec,paramFO.nb_dir_tx,paramFO.nbRX);
    for iboucle=1:Nb_Data_Dec
        for jboucle=1:paramFO.nb_dir_tx %nb faisceaus Tx
            for kboucle=1:paramFO.nbRX
                Pbruit_RX(iboucle,jboucle,kboucle) = Analyse(iboucle).bruitRX(jboucle,kboucle);
            end
        end
    end


    for zboucle=1:paramFO.nb_dir_tx
        figure(100+zboucle),
        hold on,
        title('Bruit ')
        for iboucle=1:paramFO.nbRX

            subplot(4,4,iboucle)
            hold on;
            plot(squeeze(Pbruit_RX(:,zboucle,iboucle)),'mx-')
            axis xy;
            axis([0 Nb_Data_Dec+1 40 80])
            title(['RX no ' num2str(iboucle)])



        end
        sgtitle(['Puissance de bruit par RX - Azimut  ' num2str(rad2deg(paramFO.tab_dir_ffc_tx(zboucle,1))) ' °' ])
    end


    % SNR par RX
    SNR_RX = zeros(Nb_Data_Dec,paramFO.nb_dir_tx,paramFO.nbRX);

    for iboucle=1:Nb_Data_Dec
        for jboucle=1:paramFO.nb_dir_tx %nb faisceaus Tx
            for kboucle=1:paramFO.nbRX
                SNR_RX(iboucle,jboucle,kboucle) = Analyse(iboucle).SNR_cible_RX(jboucle,kboucle);
            end
        end
    end


    for zboucle=1:paramFO.nb_dir_tx
        figure(200+zboucle),
        hold on,
        title('SNR ')
        for iboucle=1:paramFO.nbRX

            subplot(4,4,iboucle)
            hold on;
            plot(squeeze(SNR_RX(:,zboucle,iboucle)),'cx-')
            axis xy;
            axis([0 Nb_Data_Dec+1 10 100])
            title(['RX no ' num2str(iboucle)])

        end
        sgtitle(['SNR par RX - Azimut  ' num2str(rad2deg(paramFO.tab_dir_ffc_tx(zboucle,1))) ' °' ])
    end
end

%% affichage matrice de covariance
for iboucle=1:Nb_Data_Dec
    figure(300+iboucle)
%     hold on
    if paramFO.nb_dir_tx == 1
        imagesc(10*log10(abs(Analyse(iboucle).R_covar))-max(max(10*log10(abs(Analyse(iboucle).R_covar)))))
        title(['Matrice de covariance - azimut : ' num2str(rad2deg(paramFO.tab_dir_ffc_tx(1))) ' °'])
        clim([-20 0]);
        colorbar;
    else
        for jboucle=1:paramFO.nb_dir_tx
            subplot(4,4,jboucle);


            imagesc(10*log10(abs(squeeze(Analyse(iboucle).R_covar(:,:,jboucle))))-max(max(10*log10(abs(squeeze(Analyse(iboucle).R_covar(:,:,jboucle)))))))
            title(['azimut : ' num2str(rad2deg(paramFO.tab_dir_ffc_tx(jboucle))) ' °'])
            clim([-20 0]);
        end
        sgtitle('Matrice de covariance bth')

    end
end

%% ================== SPECTRES CAPON & MUSIC POUR TOUTES LES RAFALES ==================
% Faisceau TX sur lequel on travaille
ib_tx = 3;  % à adapter : entre 1 et paramFO.nb_dir_tx

% Sécurité sur les bornes
ib_tx = min(max(1, ib_tx), paramFO.nb_dir_tx);

% Nombre de rafales effectivement traitées
nbRaf = Nb_Data_Dec;

% On recharge la géométrie antennaire (même fonction que dans Processing)
[antenne_RX_loc] = fLoadAntennaRX(CST, radar_type, visu);

% Grille d'angles (en degrés / radians)
theta_deg = linspace(-90, 90, 721);   % par ex. pas de 0.25°
theta_rad = deg2rad(theta_deg);
nTheta    = numel(theta_deg);

% Exemple de matrice de covariance pour connaître le nombre d'antennes
R_example = Analyse(1).R_covar(:,:,ib_tx);
M = size(R_example, 1);               % nb d'antennes RX

% Pré-calcul des vecteurs de pointage pour tous les angles (M x nTheta)
A = zeros(M, nTheta);
for it = 1:nTheta
    % On suppose élévation = 0 rad (à adapter si besoin)
    A(:, it) = spatialSteeringVector(antenne_RX_loc.ANT, paramFO.lambda, theta_rad(it), 0);
    % Si tu veux inclure la pondération sigma comme dans le BF :
    % A(:, it) = A(:, it) .* antenne_RX_loc.pondS;
end

% Allocation des matrices de résultats : [rafale x angle]
P_capon_all = zeros(nbRaf, nTheta);
P_music_all = zeros(nbRaf, nTheta);

% Paramètre MUSIC : nombre de sources (cible + brouilleur)
n_sources = 2;

for k = 1:nbRaf
    R = Analyse(k).R_covar(:,:,ib_tx);
    
    % Hermitisation pour être sûr (symétrie Hermitienne)
    R = (R + R')/2;
    
    % Diagonal loading pour stabiliser l'inversion
    dl  = 1e-3 * trace(R) / M;
    Rdl = R + dl * eye(M);
    
    % Inverse (pour CAPON)
    Ri = inv(Rdl);
    
    % ================= CAPON =================
    % denom_capon(theta) = a^H R^{-1} a
    % On vectorise : chaque colonne de A est un a(theta)
    denom_capon = sum(conj(A) .* (Ri * A), 1);   % 1 x nTheta
    denom_capon = real(denom_capon);             % on enlève les résidus imaginaires
    denom_capon = max(denom_capon, eps);         % éviter division par 0
    P_capon = 1 ./ denom_capon;                  % 1 x nTheta
    
    % ================= MUSIC =================
    % Décomposition en valeurs propres
    [EVEC, EVAL] = eig(Rdl);
    lambda = diag(EVAL);
    [lambda_sorted, idx] = sort(lambda, 'descend'); %#ok<ASGLU>
    EVEC = EVEC(:, idx);
    
    % Sélection du nombre de sources (borné entre 1 et M-1)
    n_src_eff = min(max(1, n_sources), M-1);
    En = EVEC(:, n_src_eff+1:end);   % sous-espace bruit
    
    % proj = E_n^H a(theta) pour tous les theta
    Proj = En' * A;                  % (M-n_src_eff) x nTheta
    denom_music = sum(abs(Proj).^2, 1);   % 1 x nTheta
    denom_music = max(denom_music, eps);
    P_music = 1 ./ denom_music;           % 1 x nTheta
    
    % ================= Normalisation & passage en dB =================
    P_capon_all(k, :) = 10*log10(P_capon / max(P_capon));
    P_music_all(k, :) = 10*log10(P_music / max(P_music));
end

%% ====== AFFICHAGES : RAFALE vs AZIMUT ======

% CAPON
hFigCapon = figure;
imagesc(theta_deg, 1:nbRaf, P_capon_all);
set(gca, 'YDir', 'normal');
xlabel('Azimut (deg)');
ylabel('Indice de rafale');
title(sprintf('Spectre CAPON (TX %d)', ib_tx));
colorbar;
caxis([-30 0]);   % échelle en dB, à ajuster si besoin

% MUSIC
hFigMusic = figure;
imagesc(theta_deg, 1:nbRaf, P_music_all);
set(gca, 'YDir', 'normal');
xlabel('Azimut (deg)');
ylabel('Indice de rafale');
title(sprintf('Spectre MUSIC (TX %d)', ib_tx));
colorbar;
caxis([-40 0]);   % échelle en dB, à ajuster si besoin

% Option : marquer l'azimut théorique de ce faisceau TX (si dispo)
if size(paramFO.tab_dir_ffc_tx,1) >= ib_tx
    th_tx_deg = rad2deg(paramFO.tab_dir_ffc_tx(ib_tx,1));

    % Sur la figure CAPON
    figure(hFigCapon); hold on;
    yL = ylim;
    plot([th_tx_deg th_tx_deg], yL, 'w--', 'LineWidth', 1.2);
    ylim(yL);

    % Sur la figure MUSIC
    figure(hFigMusic); hold on;
    yL = ylim;
    plot([th_tx_deg th_tx_deg], yL, 'w--', 'LineWidth', 1.2);
    ylim(yL);
end

%% ================================
% VALEURS PROPRES DE LA MATRICE DE COVARIANCE AU COURS DES RAFALES
% ================================

disp('Calcul des valeurs propres de la matrice de covariance...')

NbVP = size(Analyse(1).R_covar,1);   % nombre de valeurs propres
VP = zeros(Nb_Data_Dec, NbVP);

if paramFO.nb_dir_tx == 1
    % === CAS 1 FAISCEAU ===
    for iboucle = 1:Nb_Data_Dec
        R = Analyse(iboucle).R_covar;
        vals = sort(eig(R), 'descend');   % valeurs propres triées
        VP(iboucle, :) = real(vals);
    end

    figure;
    hold on;
    grid on;
    plot(10*log10(VP), 'LineWidth', 1.5);
    xlabel('Numéro de rafale');
    ylabel('Valeurs propres (dB)');
    title('Évolution des valeurs propres de la matrice de covariance');
    legend(arrayfun(@(x) ['VP ' num2str(x)], 1:NbVP, 'UniformOutput', false));

else
    % === CAS MULTI-FAISCEAUX ===
    for jboucle = 1:paramFO.nb_dir_tx
        VP = zeros(Nb_Data_Dec, NbVP);

        for iboucle = 1:Nb_Data_Dec
            R = squeeze(Analyse(iboucle).R_covar(:,:,jboucle));
            vals = sort(eig(R), 'descend');
            VP(iboucle, :) = real(vals);
        end

        figure;
        hold on;
        grid on;
        plot(10*log10(VP), 'LineWidth', 1.5);
        xlabel('Numéro de rafale');
        ylabel('Valeurs propres (dB)');
        title(['Valeurs propres - Azimut ' num2str(rad2deg(paramFO.tab_dir_ffc_tx(jboucle,1))) '°']);
        legend(arrayfun(@(x) ['VP ' num2str(x)], 1:NbVP, 'UniformOutput', false));
    end
end

%% ================================
%   AFFICHAGE RANGE–DOPPLER LCMV ADAPTATIF
% ================================

% Rafale à afficher (ex : 15 si brouilleur activé à partir de 15)
k_rafale_aff = 15;
ib_tx_aff    = 1;          % faisceau TX à afficher

% Sécurisation
nbRaf = numel(Analyse);
k_rafale_aff = min(max(1, k_rafale_aff), nbRaf);
ib_tx_aff    = min(max(1, ib_tx_aff), paramFO.nb_dir_tx);

% Récupération du RD adaptatif
if isfield(Analyse(k_rafale_aff), 'sig_DBF_LCMV_adapt_dB')
    RDtmp = Analyse(k_rafale_aff).sig_DBF_LCMV_adapt_dB;

    if ndims(RDtmp) == 2
        RD_LCMV = RDtmp;  % cas 1 faisceau
    else
        RD_LCMV = RDtmp(:,:,ib_tx_aff); % cas multi-faisceaux
    end
else
    error('Le champ sig_DBF_LCMV_adapt_dB n''existe pas dans Analyse.');
end

% Axes Range & Doppler
range_axis = (0:paramFO.ncd-1)*paramFO.DistResol;
speed_axis = (0:paramFO.nrec_dir-1)*paramFO.Vresol;

% ================== AFFICHAGE =====================
figure;
imagesc(speed_axis - paramFO.Vitamb/2, range_axis, fftshift(RD_LCMV, 2));
axis xy;
colormap(jet2);
colorbar;

% Choisis ton échelle selon ton dataset (à ajuster si besoin)
clim([65 155]);  

xlabel('Speed (m/s)');
ylabel('Range (m)');

title(sprintf('Range–Doppler LCMV adaptatif\nRafale %d   |   TX %d', ...
      k_rafale_aff, ib_tx_aff));




%% ============================================================
%   MUSIC / CAPON APRES SUPPRESSION ADAPTATIVE DES BROUILLEURS
%   (visualisation d'un spectre nettoyé avec idéalement une seule raie)
% ============================================================

% Choix de la rafale & du faisceau TX à analyser
k_rafale = 22;      % à adapter (par ex. après activation du brouilleur)
ib_tx    = 1;       % FO6 -> 1 faisceau, FO2 -> adapter

% Sécurisation
nbRaf = numel(Analyse);
k_rafale = min(max(1,k_rafale), nbRaf);
ib_tx    = min(max(1,ib_tx),   paramFO.nb_dir_tx);

% Matrice de covariance d'origine pour cette rafale & ce TX
R = Analyse(k_rafale).R_covar(:,:,ib_tx);
M = size(R,1);

% Hermitisation + diagonal loading
R = (R + R')/2;
dl = 1e-3 * trace(R)/M;
Rdl = R + dl*eye(M);

% ========= Grille angulaire & vecteurs de pointage ==========
theta_deg = -90:0.5:90;
theta_rad = deg2rad(theta_deg);
nTheta    = numel(theta_deg);

[antenne_RX_loc] = fLoadAntennaRX(CST, radar_type, visu);

A_steer = zeros(M,nTheta);
for it = 1:nTheta
    a_tmp = spatialSteeringVector(antenne_RX_loc.ANT, paramFO.lambda, theta_rad(it), 0);
    A_steer(:,it) = a_tmp .* antenne_RX_loc.pondS(:);
end

% ========= MUSIC AVANT NETTOYAGE (pour comparaison) ==========
[EVEC, EVAL] = eig(Rdl);
[lambda, idx] = sort(diag(EVAL), 'descend');
EVEC = EVEC(:,idx);

n_sources = 2;                          % cible + 1 brouilleur (à adapter)
n_sources = min(max(1,n_sources),M-1);
En = EVEC(:,n_sources+1:end);

Proj = En' * A_steer;                   % (M-n_sources) x nTheta
denom = sum(abs(Proj).^2,1);
denom = max(denom,eps);
P_music_orig = 1 ./ denom;
P_music_orig_dB = 10*log10(P_music_orig / max(P_music_orig));

% ========= Détection automatique des raies (DOA) ==========
[pk, loc_idx] = findpeaks(P_music_orig_dB, theta_deg, ...
                          'SortStr','descend', ...
                          'MinPeakProminence', 3);   % seuil à ajuster

theta_peaks = loc_idx(:);      % DOA estimées en degrés
nPeaks      = numel(theta_peaks);

if nPeaks == 0
    warning('Pas de raie nette détectée par MUSIC pour cette rafale.');
    return;
end

% Cible = pic le plus proche de la direction TX théorique
th_tx_deg = rad2deg(paramFO.tab_dir_ffc_tx(ib_tx,1));
[~,i_cible] = min(abs(theta_peaks - th_tx_deg));
theta_cible_deg = theta_peaks(i_cible);

% Brouilleurs = les autres pics, en dehors d'une zone de garde autour de la cible
zone_cible = 1/4;             % ±5° autour de la cible
theta_jam_deg = theta_peaks;
theta_jam_deg(i_cible) = [];

theta_jam_deg = theta_jam_deg( abs(theta_jam_deg - theta_cible_deg) > zone_cible );

% ========= Construction du projecteur qui enlève le(s) brouilleur(s) ==========
if isempty(theta_jam_deg)
    warning('Aucun brouilleur nettement détecté. On ne modifie pas R.');
    R_clean = Rdl;
else
    % Matrice des steering vectors des brouilleurs
    nJam = numel(theta_jam_deg);
    A_jam = zeros(M,nJam);
    for jj = 1:nJam
        a_j = spatialSteeringVector(antenne_RX_loc.ANT, paramFO.lambda, ...
                                    deg2rad(theta_jam_deg(jj)), 0);
        A_jam(:,jj) = a_j .* antenne_RX_loc.pondS(:);
    end

    % Projecteur sur l'orthogonal du sous-espace brouilleur
    % P_perp = I - A_jam (A_jam^H A_jam)^(-1) A_jam^H
    G = A_jam' * A_jam;
    P_perp = eye(M) - A_jam * (G \ A_jam');

    % Covariance nettoyée
    R_clean = P_perp * Rdl * P_perp';
    R_clean = (R_clean + R_clean')/2;
end

% ========= MUSIC APRES NETTOYAGE ==========
[EVECc, EVALc] = eig(R_clean);
[lambda_c, idx_c] = sort(diag(EVALc), 'descend');
EVECc = EVECc(:,idx_c);

En_c = EVECc(:, n_sources+1:end);

Proj_c = En_c' * A_steer;
denom_c = sum(abs(Proj_c).^2,1);
denom_c = max(denom_c,eps);
P_music_clean = 1 ./ denom_c;
P_music_clean_dB = 10*log10(P_music_clean / max(P_music_clean));

% (Optionnel) CAPON après nettoyage
Ri_clean = inv(R_clean);
den_capon = sum(conj(A_steer) .* (Ri_clean * A_steer), 1);
den_capon = real(max(den_capon,eps));
P_capon_clean = 1 ./ den_capon;
P_capon_clean_dB = 10*log10(P_capon_clean / max(P_capon_clean));

%% ========= AFFICHAGES ==========
% MUSIC avant / après sur la même figure
figure;
plot(theta_deg, P_music_orig_dB, 'b-', 'LineWidth',1.2); hold on;
plot(theta_deg, P_music_clean_dB, 'r-', 'LineWidth',1.2);
yline(-3,'k:');
grid on;
xlabel('Azimut (deg)');
ylabel('Spectre (dB, normalisé)');
title(sprintf('MUSIC avant / après antibrouillage (Rafale %d, TX %d)', ...
       k_rafale, ib_tx));
legend('MUSIC original','MUSIC nettoyé','Location','best');

% Marque la DOA de la cible et des brouilleurs estimés
xline(theta_cible_deg, 'k--', 'LineWidth',1.2);
for jj = 1:numel(theta_jam_deg)
    xline(theta_jam_deg(jj), 'm--', 'LineWidth',1.0);
end

% (Optionnel) CAPON nettoyé seul
figure;
plot(theta_deg, P_capon_clean_dB, 'g-', 'LineWidth',1.2); grid on;
xlabel('Azimut (deg)');
ylabel('Spectre CAPON (dB, normalisé)');
title(sprintf('CAPON après antibrouillage (Rafale %d, TX %d)', ...
       k_rafale, ib_tx));

%% =================== PARTIE LCMV ADAPTATIF (Cas 1 Faisceau ou Multi-Faisceaux) ===================================

% --- Paramètres communs ---
azimut_annulation_deg = 4.5; % Azimut connu du brouilleur
azimut_annulation_rad = deg2rad(azimut_annulation_deg);
NbRX = paramFO.nbRX;
element_indices = (0:NbRX-1)'; 
% Fonction du vecteur de pointage (Steering Vector)
steering_vector = @(theta) exp( -1i * pi * element_indices * sin(theta) );
% Grille d'angles pour le tracé
azimuts_deg_plot = -90:0.1:90;
azimuts_rad_plot = deg2rad(azimuts_deg_plot);


if paramFO.nb_dir_tx == 1
    
    disp('Calcul et affichage du diagramme LCMV (Cas 1 Faisceau), une figure par rafale...');
    
    % --- Cas 1 Faisceau : Initialisation ---
    jboucle_cible = 1; % Le seul faisceau existant
    azimut_look = paramFO.tab_dir_ffc_tx(jboucle_cible, 1);
    
    % --- Définition des contraintes ---
    a_look = steering_vector(azimut_look);
    a_null = steering_vector(azimut_annulation_rad); 
    C = [a_look, a_null];
    g = [1; 0];

    % --- Boucle LCMV Rafale par Rafale (Cas 1 Faisceau) ---
    for iboucle = 1:Nb_Data_Dec
        % 1. Récupération de la Matrice de Covariance (sans 3ème dimension)
        R = Analyse(iboucle).R_covar; 
        
        % 2. Calcul du Vecteur de Poids Optimal LCMV
        epsilon = 1e-3 * trace(R) / size(R, 1);
        R_inv = inv(R + epsilon * eye(size(R)));
        
        R_inv_C = R_inv * C;
        inv_term = inv(C' * R_inv_C); 
        w_lcmv = R_inv_C * inv_term * g;
        
        % 3. Calcul du Diagramme 
        P_lcmv_dB_current = zeros(size(azimuts_deg_plot));
        for idx = 1:length(azimuts_rad_plot)
            a_theta = steering_vector(azimuts_rad_plot(idx));
            gain_linear = abs(w_lcmv' * a_theta)^2;
            P_lcmv_dB_current(idx) = 10 * log10(gain_linear);
        end
        
        % 4. Normalisation
        P_lcmv_dB_current = P_lcmv_dB_current - max(P_lcmv_dB_current);
        
        % 5. Tracé de la Figure
        figure(400 + iboucle);
        plot(azimuts_deg_plot, P_lcmv_dB_current, 'LineWidth', 2);
        
        grid on;
        xlabel('Azimut (\circ)');
        ylabel('Gain Normalisé (dB)');
        title(['Diagramme Anti-Brouillage LCMV - Rafale N° ' num2str(iboucle) ' (Cas 1 Faisceau)' ]);
        ylim([-60 0]);
        xlim([-90 90]);
        
        % Ajout de la ligne pour le zéro forcé
        line([azimut_annulation_deg azimut_annulation_deg], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
        text(azimut_annulation_deg + 2, -10, ['Null @ ' num2str(azimut_annulation_deg) '\circ'], 'Color', 'r', 'FontWeight', 'bold');
        
    end
    disp('Affichage figure par figure terminé pour le cas 1 Faisceau.');


else % === CAS MULTI-FAISCEAUX (paramFO.nb_dir_tx > 1) ===
    
    disp('Calcul et affichage du diagramme LCMV (Multi-faisceaux), en ciblant le faisceau le plus près de 0°...');
    
    % --- 1. IDENTIFICATION DU FAISCEAU CIBLE (Le plus proche de 0°) ---
    azimuts_rad = paramFO.tab_dir_ffc_tx(:, 1);
    azimuts_deg = rad2deg(azimuts_rad);
    [~, jboucle_cible] = min(abs(azimuts_deg)); % Indice du faisceau le plus proche de 0°
    azimut_look = azimuts_rad(jboucle_cible);
    
    disp(['Faisceau cible sélectionné : N° ' num2str(jboucle_cible) ' (Azimut: ' num2str(azimuts_deg(jboucle_cible)) '°)']);

    % --- 2. Définition des contraintes ---
    a_look = steering_vector(azimut_look);           % Contrainte 1 : Maintenir le gain (Gain=1)
    a_null = steering_vector(azimut_annulation_rad); % Contrainte 2 : Forcer le zéro (Gain=0)
    C = [a_look, a_null];
    g = [1; 0];
    
    % --- 3. Boucle LCMV Rafale par Rafale (Cas Multi-Faisceaux) ---
    for iboucle = 1:Nb_Data_Dec
        % Récupération de la Matrice de Covariance pour le faisceau cible (jboucle_cible)
        R = squeeze(Analyse(iboucle).R_covar(:,:,jboucle_cible)); 
        
        % Calcul du Vecteur de Poids Optimal LCMV
        epsilon = 1e-3 * trace(R) / size(R, 1);
        R_inv = inv(R + epsilon * eye(size(R)));
        
        R_inv_C = R_inv * C;
        inv_term = inv(C' * R_inv_C); 
        w_lcmv = R_inv_C * inv_term * g;
        
        % Calcul du Diagramme 
        P_lcmv_dB_current = zeros(size(azimuts_deg_plot));
        for idx = 1:length(azimuts_rad_plot)
            a_theta = steering_vector(azimuts_rad_plot(idx));
            gain_linear = abs(w_lcmv' * a_theta)^2;
            P_lcmv_dB_current(idx) = 10 * log10(gain_linear);
        end
        
        % Normalisation
        P_lcmv_dB_current = P_lcmv_dB_current - max(P_lcmv_dB_current);
        
        % Tracé de la Figure
        figure(400 + iboucle); 
        plot(azimuts_deg_plot, P_lcmv_dB_current, 'LineWidth', 2);
        
        grid on;
        xlabel('Azimut (\circ)');
        ylabel('Gain Normalisé (dB)');
        title(['Diagramme Anti-Brouillage LCMV - Rafale N° ' num2str(iboucle) ' (Faisceau N° ' num2str(jboucle_cible) ')' ]);
        ylim([-60 0]);
        xlim([-90 90]);
        
        % Ajout des lignes pour la visée et l'annulation
        line([azimuts_deg(jboucle_cible) azimuts_deg(jboucle_cible)], ylim, 'Color', 'b', 'LineStyle', ':', 'LineWidth', 1);
        line([azimut_annulation_deg azimut_annulation_deg], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
        text(azimut_annulation_deg + 2, -10, ['Null @ ' num2str(azimut_annulation_deg) '\circ'], 'Color', 'r', 'FontWeight', 'bold');
        
    end
    disp('Affichage figure par figure terminé pour le cas Multi-Faisceaux.');
end

