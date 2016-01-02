function [B,C,Montage] = demix_image(Images,Mix_matrix,method)
% à la fin, C est un dictionnaire (ou une sorte de liste) qui contient dans
% chaque valeur une image "démixée". B lui, contient les différents
% mélanges de sources (autant que de sources).

Nb_sources = length(Images); % Images est un dictionnaire (cell) d'images : Images{1} = lena par exemple
N_pixels = 512; % Vu qu'on en a N² à la fin, 256 (65536) c'est ambitieux
Sources = zeros(Nb_sources, N_pixels^2);
for i = 1 : Nb_sources
    Sources(i,:) = im2double(reshape( imresize(Images{i}, [N_pixels N_pixels], 'bicubic'), 1 , N_pixels^2));
    Sources(i,:) = Sources(i,:) - mean(Sources(i,:)); 
end

% On mélange
Mixed_images = Mix_matrix * Sources;
[Mixed_images,~,M_white] = whitening(Mixed_images);
B = cell(1,Nb_sources);
for i = 1 : Nb_sources
    B{i} = mat2gray(reshape(Mixed_images(i,:), N_pixels, N_pixels));
    B{i} = floor(B{i} * 255./max(max(B{i})));
end



% On démélange
Real_M_demix = inv( M_white * Mix_matrix ); % La vraie
switch method
    case 'JADE'
        M_demix = JADE(Mixed_images,Nb_sources);
    case {'FastICA','fastica','fastICA'}
        [~ , M_demix] = fastica(Mixed_images);
    case {'HJ','hj'}
        M_demix = hj(Mixed_images,Nb_sources);
        size(M_demix)
    case {'KICA','kICA','kica'}
        M_demix =  kernel_ica(Mixed_images);
    otherwise
        warning('Mauvaise méthode sélectionnée, implémente-la ou prends ce qu on te donne !')
end
fprintf('Amari error:%s',dBach(Real_M_demix,M_demix))
Estimates = M_demix*Mixed_images;
C = cell(1,Nb_sources);
for i = 1 : Nb_sources
    C{i} = mat2gray(reshape(Estimates(i,:), N_pixels, N_pixels));
    C{i} = floor(C{i} * 255./max(max(C{i})));
end
Montage = [];
for i = 1 : Nb_sources
    Montage = [Montage,[imresize(Images{i},[N_pixels N_pixels],'bicubic');uint8(B{i});uint8(C{i})]];
end
end