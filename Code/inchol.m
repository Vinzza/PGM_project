function [P,G,i] = inchol(K,eta)
% Compute the incomplete cholesky factorization given
% - K : the initial matrix
% - eta : the error acceptable on the decomposition
%
% Gives P a permutation matrix, G the lower triangular matrix and i the
% number of column of G such that
% norm( P*K*P' - G*G' ) < eta

	[M,N] = size(K);

    if N ~= M
        error('K is not a square matrix');
    end
    
    i = 0;
	P = eye(N);
	G = diag(diag(K));
	K1 = K;

    while sum(diag(G)) > eta && i < N
		i = i+1;
        h = diag(G);
        [~,j] = max(h(i:N)); j = j+i-1;
        if i ~= j
            % Update the permutation matrix
            P(i,i) = 0; P(j,j) = 0; P(i,j) = 1; P(j,i) = 1;

            % Apply transposition on K1
            dum = K1(:,i); K1(:,i) = K1(:,j); K1(:,j)= dum;     
            dum = K1(i,:); K1(i,:) = K1(j,:); K1(j,:)=dum;
            % Apply transposition on G
            dum = G(i,1:i); G(i,1:i) = G(j,1:i); G(i,1:i) = dum;
        end

        G(i,i) = sqrt(K1(i,i));
		G(i+1:N,i) = ( K1(i+1:N,i) - G(i+1:N,1:i-1)*G(i,1:i-1)' ) ./ G(i,i);

        for k = i+1:N
            G(k,k) = K(k,k) - sum( G(k,1:i).^2 );
            % g = K(k,k) - sum(G(k,1:i).^2);
            % G(k,k) = g;
        end
    end
    
    G = G( :, 1:i);
    
end