function [mean_err, v_err] = iter_test_ICA( source_size, nb_source, ICAfun, ...
                                             nb_iter, type, param, nb_outliers )
% ICAfun must take two arguments: the mixing signal and the number of sources,
% and must return the demixing matrix. => Demix = ICAfun( Signal, n );

  outliers_error = 5;
  accepted_error = nb_iter;
  
v_err = zeros( nb_iter, 1 );
i = 0;

while( i < nb_iter && accepted_error > 0 )
    fprintf('\n------------------\nIteration number %i:\n', i);

    initial_signal = create_signal( nb_source, source_size, type, param );

    % Adding outliers
    if( nb_outliers > 0 )
        perm = randperm(source_size);
        for l = 1:nb_outliers
            dum = randi([1 nb_source]);
            initial_signal(dum,perm(l)) = initial_signal(dum,perm(l)) + ...
                        outliers_error * ( 2*randi([0 1]) - 1 );
        end
    end
    
    
    % Create mixing matrix
    M_mix = 1+rand(nb_source);

    % Mix signal
    mix_signal = M_mix * initial_signal;

    % Whitening of the mixed signals
    [mix_white_signal, ~, M_white] = whitening( mix_signal );

    try
        % applying ICA
        M_demix = ICAfun( mix_white_signal, nb_source );
        
        % Evaluation of the demixing
        Real_M_demix = inv( M_white * M_mix );
        i = i + 1;
        v_err(i) = dBach( Real_M_demix, M_demix );
        fprintf('\nerr = %2.4f', v_err(i));
    catch
        fprintf(' ... error (%i left)... ', accepted_error);
        accepted_error = accepted_error - 1;
    end
end

v_err = v_err(1:i);

mean_err = mean(v_err);
fprintf( '\n\n(raised error accepted %i)\n---\nMean Error = %2.4f\n\n', ...
                                             nb_iter-accepted_error, mean_err );

end