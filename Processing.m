function [DBF, paramFO, Analyse] = Processing(CST, rep_fic, nom_fic, radar_type, Debut, Nb_Data_Dec, visu)

DBF     = [];
Analyse = [];

% Remove file extension if exists
% -------------------------------------------------
nom_fic=strrep(nom_fic, ".data", "");
nom_fic=strrep(nom_fic, ".index", "");

LectureDataSg = false;

% Get screen size
% --------------------------------------------------------------------------------
visu.scrSize            = get(0, 'screensize');

% Panneau antennaire
% (Position des antennes et coefficients de ponderation permettant de
%  former les voies S et delta 2D)
% -------------------------------------------------------------------
[antenne_RX] = fLoadAntennaRX(CST, radar_type, visu);

% Load waveform and processing parameters
% -------------------------------------------------------------------------
antenne_RX.nbRX = 16; % number of RX antenna
[paramFO] = fInitParamFO(CST, rep_fic, nom_fic, radar_type, antenne_RX, visu); % Waveform parameters
paramT     = fInitParamT(CST, paramFO, antenne_RX, visu);
% -------------------------------------------------------------------------
% Load steering vectors + calibration
% -------------------------------------------------------------------------

if paramFO.num_fo == 1
    load steeringCalTables2D_FO1.mat
elseif paramFO.num_fo == 2
    load steeringCalTables2D_FO2.mat
end


    tabCalRX = (paramFO.tab_cal_RX).';
        
        % Composer la table des vecteurs de pointage pour chacun 
        % des faisceaux formes a la reception
        % ------------------------------------------------------
%         steer_table_cal_S = zeros(paramFO.nbRX, paramT.nb_dir_ffc_rx);
        for iFRX=1:paramT.nb_dir_ffc_rx
            CircDir = paramT.tab_dir_ffc_rx(iFRX, 1);
            ElDir   = paramT.tab_dir_ffc_rx(iFRX, 2);
            % Le vecteur de pointage
            [steer_table]     = spatialSteeringVector(antenne_RX.ANT, paramFO.lambda, CircDir, ElDir);
            % Calibration
            steer_table_cal = steer_table;
            % Ponderation Sigma/Delta
            steeringCalTables2D.S(:,iFRX) = steer_table_cal.*antenne_RX.pondS;
             
        end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                     Read data file
%                  Signal processing
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Gestion des indices de debut et de fin des bundles a traiter.
% ---------------------------------------------------------------------
if( isempty(Debut) || Debut<0 )
    debut = 1;
end
if( Nb_Data_Dec<0 || (Debut+Nb_Data_Dec)>(paramFO.index.nBundles-1) )
    Nb_Data_Dec=paramFO.index.nBundles-Debut;
end

% Go to 1st bundle to process
% --------------------------------------
if( GoToDataFileBundle(paramFO.fileDataID, Debut) && ((LectureDataSg==false) || ...
        GoToDataSgFileBundle(fileDataSgID, structIndexSg, Debut)) )

    % for movie
    if(visu.sauveVideoCarteXYZ)
        v = VideoWriter("movie.avi");
        v.FrameRate = 5;
        open(v)
    end

    % Process "Nb_Data_Dec" bundles.
    % --------------------------------------
    txtRetourDebutChamp = repmat('\b',1, 10); % Pour l'affichage console uniquement
    fprintf('[%s] : ', nom_fic);
    for nbDataDec=1:Nb_Data_Dec

        iBundle = Debut + (nbDataDec-1); % Numero du "bundle" par rapport au debut du fichier
        txt = sprintf('%4d/%-4d', iBundle, paramFO.index.nBundles);
        if( nbDataDec>1), fprintf(1, txtRetourDebutChamp);  end
        fprintf('%10s', txt);

        % Read bundle (signal(I,Q), camera, ...)
        % ---------------------------------------------------------------
        [HrmSyn]    = ReadData(paramFO.fileDataID);

        % --------------------------------------------------
        % SIGNAL PROCESSING
        % --------------------------------------------------
        % Formating radar data
        % ------------
        sig = permute(complex(HrmSyn.Hammer.Matrix(:,:,:,1), HrmSyn.Hammer.Matrix(:,:,:,2)),[3,2,1]);
        % First dimension = number of samples = fast-time dimension
        % Second dimension = number of chrips = slow-time dimension
        % (caution : it contains differents emission beams)
        % Third dimension = number of RX antenna

        % Reformating radar data and calibration
        sig = reshape(sig,paramFO.ncd*paramFO.nrec_tot,antenne_RX.nbRX);
        sig_desentr = zeros(paramFO.ncd*paramFO.nrec_dir,antenne_RX.nbRX,paramFO.nb_dir_tx);

        if paramFO.ptgs_entrelaces
            for iboucle=1:paramFO.nb_dir_tx
                for jboucle=1:antenne_RX.nbRX
                    for kboucle=1:paramFO.nrec_dir
                        sig_desentr((kboucle-1)*paramFO.ncd+1:kboucle*paramFO.ncd,jboucle,iboucle) = conj(tabCalRX(jboucle)) * sig((kboucle-1)*(paramFO.ncd*paramFO.nb_dir_tx)+1+paramFO.ncd*(iboucle-1):(kboucle-1)*(paramFO.ncd*paramFO.nb_dir_tx)+paramFO.ncd+paramFO.ncd*(iboucle-1),jboucle);
                    end
                end
            end
        else
            for iboucle=1:paramFO.nb_dir_tx
                for jboucle=1:antenne_RX.nbRX
                    sig_desentr(:,jboucle,iboucle) = conj(tabCalRX(jboucle)) * sig((iboucle-1)*paramFO.ncd*paramFO.nrec_dir+1:iboucle*paramFO.ncd*paramFO.nrec_dir,jboucle);
                end
            end
        end

        sig_desentr = reshape(sig_desentr,paramFO.ncd,paramFO.nrec_dir,antenne_RX.nbRX,paramFO.nb_dir_tx);
        % Dimension 1 = number of samples = fast-time dimension
        % Dimension 2 = number of chrips in one beam = slow-time dimension
        % Dimension 3 = number of RX antenna
        % Dimension 4 = number of TX beams

   

        % Range FFT
        sig_FFT1 = zeros(paramFO.ncd,paramFO.nrec_dir,antenne_RX.nbRX,paramFO.nb_dir_tx);
        w = hamming(paramFO.ncd);
%         w = ones(paramFO.ncd);
        for iboucle=1:paramFO.nb_dir_tx
            for jboucle=1:antenne_RX.nbRX
              for kboucle=1:paramFO.nrec_dir
                sig_FFT1(:,kboucle,jboucle,iboucle) = fft(w.*sig_desentr(:,kboucle,jboucle,iboucle));
              end
            end
        end

%         Doppler FFT
        sig_FFT1_FFT2 = zeros(paramFO.ncd,paramFO.nrec_dir,antenne_RX.nbRX,paramFO.nb_dir_tx);
%         w = hamming(paramFO.nrec_dir)';
%         w = ones(paramFO.nrec_dir)';
        nbre_point=paramFO.nrec_dir;
        w=(0.35875-0.48829.*cos(2.*pi./nbre_point.*(0:nbre_point-1))...
        +0.14128.*cos(4.*pi./nbre_point.*(0:nbre_point-1))...
        -0.01168.*cos(6.*pi./nbre_point.*(0:nbre_point-1)));
        for iboucle=1:paramFO.nb_dir_tx
            for jboucle=1:antenne_RX.nbRX
              for kboucle=1:paramFO.ncd
                sig_FFT1_FFT2(kboucle,:,jboucle,iboucle) = fft(w.*sig_FFT1(kboucle,:,jboucle,iboucle));
              end
            end
        end
                % Complex to dB power
        sig_FFT1_FFT2_2 = abs(sig_FFT1_FFT2).^2;
        sig_FFT1_FFT2_dB = 20*log10(abs(sig_FFT1_FFT2));
%         
%         for iboucle=1:16
%         Analyse(nbDataDec).norm1 = (1./squeeze(mean(mean(abs(sig_FFT1_FFT2(250:450,10:26,:,8)),2),1)));
%         end

        % Adaptive beamforming (a faire)
        % il faut utiliser les donnees *sig_FFT1_FFT2*
        
        % Beamforming
        sig_DBF_FFT1_FFT2 = zeros(paramFO.ncd,paramFO.nrec_dir,paramFO.nb_dir_tx);

        for iboucle=1:paramFO.nb_dir_tx
            for jboucle=1:paramFO.ncd
                for kboucle=1:paramFO.nrec_dir
                    sig_DBF_FFT1_FFT2(jboucle,kboucle,iboucle) = steeringCalTables2D.S(:,iboucle)'*squeeze(sig_FFT1_FFT2(jboucle,kboucle,:,iboucle));
                end
            end
        end
    
       

        % Complex to dB power
        sig_DBF_FFT1_FFT2_2 = abs(sig_DBF_FFT1_FFT2).^2;
        sig_DBF_FFT1_FFT2_dB = 20*log10(abs(sig_DBF_FFT1_FFT2));
        
        % calcul du niveau de bruit pour tous les pointages TX
        
        Analyse(nbDataDec).bruit = squeeze(squeeze(mean(mean(sig_DBF_FFT1_FFT2_dB(paramT.bruit.iCDmin:paramT.bruit.iCDmax,paramT.bruit.iDopplerMin:paramT.bruit.iDopplerMax,:),1),2)));

        % estimation matrice de covariance
        R_covar = zeros(paramFO.nbRX,paramFO.nbRX,paramFO.nb_dir_tx);

        for iboucle=1:paramFO.nb_dir_tx
            
            data_entrainement = sig_FFT1_FFT2(paramT.bruit.iCDmin:paramT.bruit.iCDmax,paramT.bruit.iDopplerMin:paramT.bruit.iDopplerMax,:,iboucle);
            vect_entrainement = reshape(data_entrainement,size(data_entrainement,1)*size(data_entrainement,2),size(data_entrainement,3)).';
            R_covar(:,:,iboucle) = (1/size(vect_entrainement,2)) * (vect_entrainement*vect_entrainement');
        end
        Analyse(nbDataDec).R_covar = R_covar;

        % calcul SNR reflecteur
        num_cd_cible = visu.num_cd_cible;
        SNR_cible=zeros(paramFO.nb_dir_tx,1);
        for iboucle=1:paramFO.nb_dir_tx
        SNR_cible(iboucle) = sig_DBF_FFT1_FFT2_dB(num_cd_cible,1,iboucle) - Analyse(nbDataDec).bruit(iboucle);
        end
        Analyse(nbDataDec).SNR_cible = SNR_cible;

        % calcul du niveau de bruit pour les voies RX
        for iboucle=1:paramFO.nb_dir_tx
             for jboucle=1:paramFO.nbRX
                 Analyse(nbDataDec).bruitRX(iboucle,jboucle) = squeeze(squeeze(mean(mean(sig_FFT1_FFT2_dB(paramT.bruit.iCDmin:paramT.bruit.iCDmax,paramT.bruit.iDopplerMin:paramT.bruit.iDopplerMax,jboucle,iboucle),1),2)));
                 Analyse(nbDataDec).SNR_cible_RX(iboucle,jboucle) = sig_FFT1_FFT2_dB(num_cd_cible,1,jboucle,iboucle) - Analyse(nbDataDec).bruitRX(iboucle,jboucle);
             end
        end




        % -------------------------------------------------------------------------
        % -------------------------------------------------------------------------
        % DISPLAY
        % -------------------------------------------------------------------------
        % -------------------------------------------------------------------------
        % Display zero Doppler range profile
        if visu.zeroDoppler > 0
            figure(visu.zeroDoppler);
            plot(sig_DBF_FFT1_FFT2_dB(:,1))
            xlabel('Case distance')
            ylabel('niveau (dB)')
            title('Profil distance Doppler 0')
        end

        % Display range Doppler map for 8 beams out of 16
        if visu.rangeDoppler > 0
            figure(visu.rangeDoppler);
            ind = 1;
            if paramFO.num_fo == 6
                    range_axis = (0:paramFO.ncd-1)*paramFO.DistResol;
                    speed_axis = (0:paramFO.nrec_dir-1)*paramFO.Vresol;
                    imagesc(speed_axis-paramFO.Vitamb/2,range_axis, fftshift(sig_DBF_FFT1_FFT2_dB,2));
                    axis xy;
                    colormap(jet2);
                    clim([65 170]);
                    title(['azimut : ' num2str(rad2deg(paramFO.tab_dir_ffc_tx(iboucle))) ' °'])
                    xlabel('Speed (m/s)')
                    ylabel('Range (m)')
            else
                for iboucle=1:2:16
                    subplot(2,4,ind)
                    range_axis = (0:paramFO.ncd-1)*paramFO.DistResol;
                    speed_axis = (0:paramFO.nrec_dir-1)*paramFO.Vresol;
                    imagesc(speed_axis-paramFO.Vitamb/2,range_axis, fftshift(sig_DBF_FFT1_FFT2_dB(:,:,iboucle),2));
                    axis xy;
                    colormap(jet2);
                    clim([65 155]);
                    title(['azimut : ' num2str(rad2deg(paramFO.tab_dir_ffc_tx(iboucle))) ' °'])
                    xlabel('Speed (m/s)')
                    ylabel('Range (m)')
                    ind = ind + 1;
                end
            end
        end
        % Display angle-range map of moving objects

        if visu.global > 0
            if paramFO.num_fo ~= 6
            max_sig_S=zeros(paramFO.ncd,paramFO.nb_dirCirc_tx);

            for iFRXc=1:paramFO.nb_dirCirc_tx
                [max_sig_S(:, iFRXc)]=max(squeeze(sig_DBF_FFT1_FFT2_dB(:,3:31,iFRXc, :)), [], [2,3]);
            end

            R=repmat(paramFO.tab_dist',[1 16]);
            Theta=repmat(paramFO.az_ffc_tx.',[paramFO.ncd,1]);
            X=R.*sin(Theta);
            Y=R.*cos(Theta);

            fig = figure(visu.global);
            fig.Position = [100 100 1600 900];

            subplot(1,2,1)
            img1D=uint8(HrmSyn.Camera(1,1).data);
            img=reshape(img1D,3,1920,1080);
            img2=permute(img,[3,2,1]);
            img3=img2;
            img3(:,:,1)=img2(:,:,3);
            img3(:,:,2)=img2(:,:,2);
            img3(:,:,3)=img2(:,:,1);
            imagesc(img3);
            subplot(1,2,2)
            surf(X,Y,fliplr(max_sig_S));
            ylim([0 20]);
            xlim([-2 2]);
            xlabel('Azimuth (m)');
            ylabel('Range (m)');
            title("Range angle map");
            clim([75 150]);
            colormap(jet);
            colorbar;
            shading interp;
            %axis equal;
            view(0,90);
            if(visu.sauveVideoCarteXYZ)
                frame = getframe(gcf);
                writeVideo(v,frame)
            end
            end
        end


        % Display camera video
        % ----------------------------------------------------
        if( visu.camera > 0 )
            fPlotCamera(HrmSyn, visu.camera);
        end

    end % End of loop on bundles
    fprintf(1, '\n');

    % close movie
    if(visu.sauveVideoCarteXYZ)
        close(v);
    end

else % "Else" of access test of 1st bundle to process
    fprintf(1, "Error : Impossible to read the bundle %d in the file %s\n", Debut, nom_fic);
end
