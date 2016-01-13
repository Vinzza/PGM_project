function [v_err] = plot_ICA_error_outliers( source_size, nb_source, ICAfun, ...
                                             nb_iter, type, param, outliers_nbs )

    v_err = zeros(size(outliers_nbs));

    for i = 1:length(outliers_nbs)
        v_err(i) = iter_test_ICA( source_size, nb_source, ICAfun, ...
                                        nb_iter, type, param, outliers_nbs(i) );
    end

    figure(1); clf; hold on;
    plot( outliers_nbs, v_err );

end
