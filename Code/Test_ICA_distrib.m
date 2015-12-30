%% TEST ICA

addpath( 'Test_ICA' );

%% %%%%%%% %% %% %  %                                      %  % %% %% %%%%%%%%%%
%%%%% %% %  %                      PARAMETRES                      %  % %% %%%%%
%%%%%%%%%% %% %% %  %                                      %  % %% %% %%%%%%% %%

source_size = 250;
% Taille de la source
% Bach utilise : 250, 1000, 2000 et 4000.
% En sachant qu'il augmente cette taille en même temps que le nombre de
% source (voir les couples utilisés dans nb_iter).

nb_source = 2;
% Nombre de source.
% Bach utilise 2 sources lorsque l'on prend une même distribution
% (type='same') et il utilise 2, 4, 8 ou 16 sources.

type = 'random';
% On a 2 types de données :
% 'same' : on prend la même distribution
%          ça correspond au tableau 1 de Bach
%          (on prendra de préférence nb_source = 2 du coup)
% 'random' : on prend des distributions différentes, prises aléatoirement
%            parmi toutes les distributions disponibles.
%            (colonne 'rand' du tableau 1 de Bach, pour nb_source = 2)

param = 5;
% Cela selectionne le numéro de la distribution pour un type 'same'.
% Dans ce cas, il faut un entier entre 1 et 12.
% Dans le cas du type 'random', ce paramètre est inutilisé.
% Ci après la correspondance entre les lettres représentant les
% distributions dans l'article de Bach et les entiers utilisés ici :
%  1 -> m
%  2 -> n
%  3 -> o
%  4 -> p
%  5 -> q
%  6 -> r
%  7 -> g
%  8 -> h
%  9 -> i
% 10 -> j
% 11 -> k
% 12 -> l

algo_name = 'JADE'; % Utiliser pour créer le fichier texte du graphe

ICAfun = @( signal, m ) JADE( signal, m );
% On doit absolument avoir une fonction qui prend le signal mixé et le
% nombre de source à extraire et qui renvoie la matrice de demixage.
% Donc forcement une fonction de la forme :
% ICAfun = @( signal, m ) fct(arg1,...,argn);
%
% Au niveau des algos, on a :
% ICAfun = @( signal, m ) JADE( signal, m );

nb_iter = 10;
% Le nombre de fois qu'on itère le calcul de l'erreur. Bach prend les
% valeurs suivantes en fonction de nb_source(m) et source_size(N) :
% (m,N) = (2,250)     =>   nb_iter = 1000;
% (m,N) = (2,1000)    =>   nb_iter = 1000;
% (m,N) = (4,1000)    =>   nb_iter = 100;
% (m,N) = (4,4000)    =>   nb_iter = 100;
% (m,N) = (8,2000)    =>   nb_iter = 50;
% (m,N) = (8,4000)    =>   nb_iter = 50;
% (m,N) = (16,4000)   =>   nb_iter = 25;
% Il faut noter que Bach utilise 16 distributions différentes alors qu'ici,
% on en a que 12.


% On considère ici le nombre de composante que l'on va modifier pour en
% faire des outliers. L'article de Bach en considère entre 0 et 25 et trace
% le graphe. (voir la fonction que j'ai certainement dû appeler un truc
% comme plot_error_outliers)
% Si on veut une seule valeur :
nb_outliers = 120;
% Si on veut tracer le graphe entier, à partir de 0 jusqu'à 'max_outliers'
% par pas de 'out_step' :
out_step = 5;
max_outliers = 25;


%% %%%%%%% %% %% %  %                                      %  % %% %% %%%%%%%%%%
%%%%% %% %  %                      EXECUTIONS                      %  % %% %%%%%
%%%%%%%%%% %% %% %  %                                      %  % %% %% %%%%%%% %%

% % Calculer l'erreur d'une seule configuration
% mean_err = iter_test_ICA( source_size, nb_source, ICAfun, nb_iter,...
%                                                      type, param, nb_outliers );

% Calculer et tracer l'erreur en fonction du nombre d'outliers
v_err = plot_ICA_error_outliers( source_size, nb_source, ICAfun, ...
                                 nb_iter, type, param, 0:out_step:max_outliers);

file_name = sprintf('errOut_%s(m%i,N%i)%s%i_iter%i(outstep%i).data', ...
             algo_name, nb_source, source_size, type, param, nb_iter, out_step);
export_plot_in_text( 0:out_step:max_outliers, v_err, file_name );

