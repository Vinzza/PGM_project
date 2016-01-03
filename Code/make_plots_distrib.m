
addpath( 'Test_ICA' );

X = -2:.1:2;

for i = 1:12

  switch i
    case 1 % m
      s = mixture_gaussian_plot( X, [1 2 2 1], [-1 -.33 .33 1], [.16  .16 .16 .16] );
    case 2 % n
      s = mixture_gaussian_plot( X, [1 2 2 1], [-1 -.2 .2 1], [.2 .3 .3 .2] );
    case 3 % o
      s = mixture_gaussian_plot( X, [1 2 2 1], [-.7 -.2 .2 .7], [.2 .3 .3 .2] );
    case 4 % p
      s = mixture_gaussian_plot( X, [1 1 2 1], [-1 .3 -.3 1.1], [.2 .2 .2 .2] );
    case 5 % q
      s = mixture_gaussian_plot( X, [1 3 2 .5], [1 3 2 .5], [.2 .3 .2 .2] );
    case 6 % r
      s = mixture_gaussian_plot( X, [1 2 2 1], [-.8 -.2 .2 .5], [.22 .3 .3 .2] );
    case 7 % g
      s = mixture_gaussian_plot( X, [1 1], [-.5 .5 ], [.15 .15] );
    case 8 % h
      s = mixture_gaussian_plot( X, [1 1], [-.5 .5 ], [.4 .4] );
    case 9 % i
      s = mixture_gaussian_plot( X, [1 1], [-.5 .5 ], [.5 .5] );
    case 10 % j
      s = mixture_gaussian_plot( X, [1 3], [-.5 .5 ], [.15 .15] );
    case 11 % k
      s = mixture_gaussian_plot( X, [1 2], [-.7 .5 ], [.4 .4] );
    case 12 % l
      s = mixture_gaussian_plot( X, [1 2], [-.7 .5 ], [.5 .5] );
  end
  name = sprintf('distrib_%i.data', i);
  export_plot_in_text( X, s, name );

end
