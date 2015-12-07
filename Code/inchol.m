function [P,G,i] = inchol(K,eta)
	[M,N] = size(K);
    if N ~= M
        'matrice pas carrée'
		RETURN 
    end
	i = 1;
	P = eye(N);
	G = diag(diag(K));
	K1 = K;
    while sum(diag(G)) > eta && i < N
        h= diag(G);
		[~,j] = max(h(i:N));j = j+i-1;
		P(i,i) = 0;
        P(j,j) = 0;
        P(i,j) = 1;
        P(j,i) = 1;
        k1 = K1(:,i);k2 = K1(:,j);
        K1(:,j)=k1;K1(:,i)=k2;
        k1 = K1(i,:);k2 = K1(j,:);
        K1(j,:)=k1;K1(i,:)=k2;
        g1 = G(i,1:i);g2 = G(j,1:i);
        G(i,1:i) = g2;G(j,1:i) = g1;
		G(i,i) = sqrt(K1(i,i));
		G(i+1:N,i) = (K1(i+1:N,i)- G(i+1:N,1:i-1)*G(i,1:i-1)')/G(i,i);
        for k =i+1:N
            g = K(k,k) - sum(G(k,1:i).^2);
            G(k,k) = g;
        end
        G
		i = i+1;
    end
end