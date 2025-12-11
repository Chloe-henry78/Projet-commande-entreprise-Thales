function [spVect] = spatialSteeringVector(antenna, lambda, circulaire, elevation)

    % Definition des vecteurs directeurs dans un repere OXYZ
    % tel que 
    %     - l'antenne est definie dans le plan XOY
    %     - OX est positif a gauche
    %     - OY est positif a haut
    %     - Z est normal au plan de l'antenne, vers l'avant
    % --------------------------------------------------------
    vect(1,1) = sin(circulaire)*cos(elevation);
    vect(2,1) = sin(elevation);
    vect(3,1) = cos(circulaire)*cos(elevation);

    % Contruire les vecteurs directeurs
    % ----------------------------------
    dpl    = (2*pi/lambda);
    x      = antenna(:, 1:3) * vect;
    x      = x * dpl;
    spVect = exp(1i * x);

end
