function export_plot_in_text( X, Y, name )

    fid = fopen( ['Export_Data/' name], 'w' );
    for i = 1:length(X)
        fprintf(fid,'%2.4f\t%2.4f\n', X(i), Y(i) );
    end
    fclose( fid );


end