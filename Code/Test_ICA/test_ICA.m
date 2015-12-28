

%% Parameter

source_size = 500;
nb_source = 4;

initial_signal = zeros( nb_source, source_size );
perm = randperm(12);

% sources
fprintf('construction of the initial signals... %2i',0);
for i = 1:nb_source
    switch perm(i)
        case 1
        initial_signal(i,:) = mixture_gaussian( source_size, [1 2 2 1], [-1 -.33 .33 1], [.16  .16 .16 .16] ); % m
        case 2
        initial_signal(i,:) = mixture_gaussian( source_size, [1 2 2 1], [-1 -.2 .2 1], [.2 .3 .3 .2] ); % n
        case 3
        initial_signal(i,:) = mixture_gaussian( source_size, [1 2 2 1], [-.7 -.2 .2 .7], [.2 .3 .3 .2] ); % o
        case 4
        initial_signal(i,:) = mixture_gaussian( source_size, [1 1 2 1], [-1 .3 -.3 1.1], [.2 .2 .2 .2] ); % p
        case 5
        initial_signal(i,:) = mixture_gaussian( source_size, [1 3 2 .5], [1 3 2 .5], [.2 .3 .2 .2] ); % q
        case 6
        initial_signal(i,:) = mixture_gaussian( source_size, [1 2 2 1], [-.8 -.2 .2 .5], [.22 .3 .3 .2] ); % r
        case 7
        initial_signal(i,:) = mixture_gaussian( source_size, [1 1], [-.5 .5 ], [.15 .15] ); % g
        case 8
        initial_signal(i,:) = mixture_gaussian( source_size, [1 1], [-.5 .5 ], [.4 .4] ); % h
        case 9
        initial_signal(i,:) = mixture_gaussian( source_size, [1 1], [-.5 .5 ], [.5 .5] ); % i
        case 10
        initial_signal(i,:) = mixture_gaussian( source_size, [1 3], [-.5 .5 ], [.15 .15] ); % j
        case 11
        initial_signal(i,:) = mixture_gaussian( source_size, [1 2], [-.7 .5 ], [.4 .4] ); % k
        case 12
        initial_signal(i,:) = mixture_gaussian( source_size, [1 2], [-.7 .5 ], [.5 .5] ); % l
    end
    fprintf('\b\b%2i',i);
end

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
Error = dBach( Real_M_demix, M_demix );
