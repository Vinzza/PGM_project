

%% Parameter

source_size = 2000;
nb_source = 4;


fprintf('construction of the initial signals...');
initial_signal = create_signal( nb_source, source_size, 'random');

%%

% Create mixing matrix
fprintf('\nconstruction of the mixing signal...\n');
M_mix = 1+rand(nb_source);

N = size(M_mix, 2);
P = size(M_mix, 1);

% Mix signal
mix_signal = M_mix * initial_signal;

% Whitening of the mixed signals
fprintf('whitening of the data...\n');
[mix_white_signal, mu, M_white] = whitening( mix_signal );

% applying ICA
fprintf('ICA algorithm...\n');
M_demix = JADE( mix_white_signal, P );


fprintf('Εvaluation...\n');

% Demix signals
Estimated_signal = M_demix * mix_signal;


% Evaluation of the demixing
Real_M_demix = inv( M_white * M_mix );
err = dBach( Real_M_demix, M_demix );
disp(err);
