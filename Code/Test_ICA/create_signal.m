function signal = create_signal( nb_source, source_size, type, param )

  if nargin < 4
    param = 0;
  end

  signal = zeros( nb_source, source_size );

  switch type
    
    % TYPE RANDOM
    case 'random'

        if( nb_source > 12 )
        error('not enough distributions available');
        end

        perm = randperm(12); % We only have 12 distributions types

        % sources
        for i = 1:nb_source
            signal(i,:) = my_gaussian_mixture( source_size, perm(i) );
        end

    % TYPE SAME
    case 'same'
        if( param > 12 || param < 1 )
            error(sprintf('distribution %i not available',param));
        end
        
        for i = 1:nb_source
            signal(i,:) = my_gaussian_mixture( source_size, param );
        end
  end

end

% LOCAL FUNCTION

function s = my_gaussian_mixture( m, nb )
% the letter correspond to the indice in the Jordan and Bach paper.
  switch nb
    case 1 % m
      s = mixture_gaussian( m, [1 2 2 1], [-1 -.33 .33 1], [.16  .16 .16 .16] );
    case 2 % n
      s = mixture_gaussian( m, [1 2 2 1], [-1 -.2 .2 1], [.2 .3 .3 .2] );
    case 3 % o
      s = mixture_gaussian( m, [1 2 2 1], [-.7 -.2 .2 .7], [.2 .3 .3 .2] );
    case 4 % p
      s = mixture_gaussian( m, [1 1 2 1], [-1 .3 -.3 1.1], [.2 .2 .2 .2] );
    case 5 % q
      s = mixture_gaussian( m, [1 3 2 .5], [1 3 2 .5], [.2 .3 .2 .2] );
    case 6 % r
      s = mixture_gaussian( m, [1 2 2 1], [-.8 -.2 .2 .5], [.22 .3 .3 .2] );
    case 7 % g
      s = mixture_gaussian( m, [1 1], [-.5 .5 ], [.15 .15] );
    case 8 % h
      s = mixture_gaussian( m, [1 1], [-.5 .5 ], [.4 .4] );
    case 9 % i
      s = mixture_gaussian( m, [1 1], [-.5 .5 ], [.5 .5] );
    case 10 % j
      s = mixture_gaussian( m, [1 3], [-.5 .5 ], [.15 .15] );
    case 11 % k
      s = mixture_gaussian( m, [1 2], [-.7 .5 ], [.4 .4] );
    case 12 % l
      s = mixture_gaussian( m, [1 2], [-.7 .5 ], [.5 .5] );
  end
end