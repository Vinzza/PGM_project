function [mean_err, v_err] = iter_test_ICA( source_size, nb_source, ICAfun, nb_iter, type, param )
% ICAfun must take two arguments: the mixing signal and the number of sources,
% and must return the demixing matrix. => Demix = ICAfun( Signal, n );

  if nargin < 6
    param = 0;
  end
  
v_err = zeros( nb_iter, 1 );

for i = 1:nb_iter
    fprintf('\n------------------\nIteration %i:', i);

    initial_signal = create_signal( nb_source, source_size, type, param );

    % Create mixing matrix
    M_mix = 1+rand(nb_source);

    P = size(M_mix, 1);

    % Mix signal
    mix_signal = M_mix * initial_signal;

    % Whitening of the mixed signals
    [mix_white_signal, ~, M_white] = whitening( mix_signal );

    % applying ICA
    M_demix = ICAfun( mix_white_signal, P );

    % Evaluation of the demixing
    Real_M_demix = inv( M_white * M_mix );
    v_err(i) = dBach( Real_M_demix, M_demix );
    fprintf('\nerr = %2.4f', v_err(i));
end

mean_err = mean(v_err);

end