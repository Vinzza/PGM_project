function c = KCCA(X,sigma,kap) %nathan : kap c'est petit kappa
	[m,N] = size(X);
	Kernel = @(x) kernelmatrix('rbf',X',[],sigma); %rajouter la possibilite de choisir d'autres kernels ?
	M = [];
    Z = [];
	for i =1:m
		K_data = Kernel(X(i,:)');
        [~,G,~] = inchol(K_data,0.01*N*kap/2.);
        [U,S,~] = svd(G);
        L = diag(S).^2;
        R = L./(L+N*kap/2);
        [Mi,~] = size(L); 
        M = [M, Mi]; %Taille du futur bloc
        Z = [Z,U*R]; %Futur demi bloc : un bloc = Z'Z
    end
    R = Z'*Z;
    R = mat2cell(R,M,M);
    for k = 1:m
        R{k,k} = eye(M(k));
    end
	R = cell2mat(R);
    c = -0.5*log(det(R));
end