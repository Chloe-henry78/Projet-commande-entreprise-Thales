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
repData='C:\Users\chloe\Downloads\Correction_steer_abf (1)\Correction_steer_abf';
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
Nb_Data_Dec             = 25; % Number of bundles to read (from index "debut");
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
visu.sauveVideoCarteXYZ = true; % to save image in movie

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

%% ======================= MUSIC 2D : AZIMUT x RAFALE =======================

disp('Calcul MUSIC 2D (Azimut x Rafale)...')

% Paramètres
azimuts_deg_music = -90:0.5:90;
azimuts_rad_music = deg2rad(azimuts_deg_music);
Nb_az = length(azimuts_rad_music);

NbRX = paramFO.nbRX;
element_indices = (0:NbRX-1)'; 

% Steering vector ULA (lambda/2)
steering_vector = @(theta) exp( -1i * pi * element_indices * sin(theta) );

% Initialisation carte MUSIC
MUSIC_MAP = zeros(Nb_Data_Dec, Nb_az);

for iboucle_music = 1:Nb_Data_Dec

    % On prend la covariance du premier faisceau Tx
    R_music = squeeze(Analyse(iboucle_music).R_covar(:,:,1));
    R_music = (R_music + R_music')/2;

    % Décomposition spectrale
    [Vecs, D] = eig(R_music);
    eigvals = real(diag(D));
    [eigvals_sorted, idxs] = sort(eigvals, 'descend');
    Vecs_sorted = Vecs(:, idxs);

    % Estimation du nombre de sources (séparation signal/bruit)
    ratios = eigvals_sorted(1:end-1) ./ (eigvals_sorted(2:end) + eps);
    [~, kmax] = max(ratios);
    num_sources = kmax;

    if num_sources < 1
        num_sources = 1;
    elseif num_sources > NbRX-1
        num_sources = NbRX-1;
    end

    % Sous-espace de bruit
    En = Vecs_sorted(:, num_sources+1:end);

    % Calcul MUSIC pour tous les azimuts
    for iaz = 1:Nb_az
        a = steering_vector(azimuts_rad_music(iaz));
        MUSIC_MAP(iboucle_music, iaz) = ...
            1 / ( real( a' * (En*En') * a ) + eps );
    end
end

% Normalisation globale
MUSIC_MAP = MUSIC_MAP ./ max(MUSIC_MAP(:));
MUSIC_MAP_dB = 10*log10(MUSIC_MAP + eps);

% ======================= AFFICHAGE 2D =======================
figure(360)
imagesc(azimuts_deg_music, 1:Nb_Data_Dec, MUSIC_MAP_dB)
axis xy
xlabel('Azimut (°)')
ylabel('Numéro de rafale')
title('MUSIC Continu : Carte Azimut × Rafale')
colorbar
clim([-40 0])
grid on



%% VALEURS PROPRES DE LA MATRICE DE COVARIANCE AU COURS DES RAFALES
disp('Calcul des valeurs propres de la matrice de covariance...')
NbVP = size(Analyse(1).R_covar,1);
VP = zeros(Nb_Data_Dec, NbVP);
if paramFO.nb_dir_tx == 1
    for iboucle = 1:Nb_Data_Dec
        R_vp = Analyse(iboucle).R_covar;
        vals = sort(eig(R_vp), 'descend');
        VP(iboucle, :) = real(vals);
    end
    figure; hold on; grid on;
    plot(10*log10(VP), 'LineWidth', 1.5);
    xlabel('Numéro de rafale'); ylabel('Valeurs propres (dB)');
    title('Évolution des valeurs propres de la matrice de covariance');
    legend(arrayfun(@(x) ['VP ' num2str(x)], 1:NbVP, 'UniformOutput', false));
else
    for jboucle = 1:paramFO.nb_dir_tx
        VP = zeros(Nb_Data_Dec, NbVP);
        for iboucle = 1:Nb_Data_Dec
            R_vp = squeeze(Analyse(iboucle).R_covar(:,:,jboucle));
            vals = sort(eig(R_vp), 'descend');
            VP(iboucle, :) = real(vals);
        end
        figure; hold on; grid on;
        plot(10*log10(VP), 'LineWidth', 1.5);
        xlabel('Numéro de rafale'); ylabel('Valeurs propres (dB)');
        title(['Valeurs propres - Azimut ' num2str(rad2deg(paramFO.tab_dir_ffc_tx(jboucle,1))) '°']);
        legend(arrayfun(@(x) ['VP ' num2str(x)], 1:NbVP, 'UniformOutput', false));
    end
end


%% PARTIE LCMV 

azimut_annulation_deg = -5; 
azimut_annulation_rad = deg2rad(azimut_annulation_deg);

% 1. Initialisation pour LCMV
iboucle_sel = 1; 
jboucle_sel = 1; 
R = squeeze(Analyse(iboucle_sel).R_covar(:,:,jboucle_sel));
azimut_look = paramFO.tab_dir_ffc_tx(jboucle_sel, 1);
NbRX = paramFO.nbRX;
element_indices = (0:NbRX-1)'; 

% Steering Vector Function (Assumes ULA lambda/2)
steering_vector = @(theta) exp( -1i * pi * element_indices * sin(theta) );

% Matrice de covariance inverse
epsilon = 1e-6 * trace(R) / size(R, 1);
R_inv = inv(R + epsilon * eye(size(R)));

% 2. Définition des Contraintes (Look Direction et Null Direction)
a_look = steering_vector(azimut_look);           % Contrainte 1 : Maintenir le gain (Gain=1)
a_null = steering_vector(azimut_annulation_rad); % Contrainte 2 : Forcer le zéro (Gain=0)

A = [a_look, a_null];
c = [1; 0]; 

% 3. Calcul du Vecteur de Poids Optimal LCMV
R_inv_A = R_inv * A;
inv_term = inv(A' * R_inv_A);
w_lcmv = R_inv_A * inv_term * c;

% 4. Calcul du Diagramme 
azimuts_deg = -90:0.1:90;
azimuts_rad = deg2rad(azimuts_deg);
P_lcmv_dB = zeros(size(azimuts_deg));

for idx = 1:length(azimuts_rad)
    a_theta = steering_vector(azimuts_rad(idx));
    gain_linear = abs(w_lcmv' * a_theta)^2;
    P_lcmv_dB(idx) = 10 * log10(gain_linear);
end

P_lcmv_dB = P_lcmv_dB - max(P_lcmv_dB); % Normalisation

% 5. Tracé du Diagramme
figure(400); % Nouvelle figure pour le LCMV
plot(azimuts_deg, P_lcmv_dB, 'LineWidth', 2);
grid on;
xlabel('Azimut (\circ)');
ylabel('Gain (dB)');
title(['Diagramme Anti-Brouillage LCMV  \circ']);
ylim([-60 0]);
xlim([-100 100]);

% Ajout de la ligne pour visualiser le zéro forcé
line([azimut_annulation_deg azimut_annulation_deg], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
text(azimut_annulation_deg + 2, -10, ['Null @ ' num2str(azimut_annulation_deg) '\circ'], 'Color', 'r', 'FontWeight', 'bold');% -------------------------------------------------------------------------