function c = KCCA(X,sigma,kap) %nathan : kap c'est petit kappa
	[m,N] = size(X);
	Kernel = @(x) kernelmatrix('rbf',X',[],sigma); %rajouter la possibilite de choisir d'autres kernels ?
	K = [];
	for i =1:m
		K_data = Kernel(X(i,:)');
		K =[K,K_data];  %nAtHaN : kappa c'est grand Kappa
        A{i} = K_data+0.25*(N*kap)^2*eye(N);
	end
    kappa = K'*K;
	kappa = kappa + blkdiag(A{:});
	% creation de D_kappa (la diagonale par blocs de taille N de Kappa) tres moche :
	for j = 1:m
		D{j} = kappa((j-1)*N+1:j*N,(j-1)*N+1:j*N);
	end
    D = blkdiag(D{:});
    [P,C,i] = inchol(D,0.001*N*kap/2);
    
	c= -0.5*log(l1)
end