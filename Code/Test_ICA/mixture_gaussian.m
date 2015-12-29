function res = mixture_gaussian( nb, prob, mu, cov )

prob = prob ./ sum(prob);
cumprob = cumsum(prob(:));
res = zeros(1,nb);

dum = rand(1,nb);
r = ones(1,nb);

for i = 1:length(prob)
    r = r + ( dum >= cumprob(i) );
end

for i = 1:length(prob)
    dum = random('normal', mu(i), cov(i), 1, nb);
    res(r==i) = dum(r==i);
end

% % Old version
% for k = 1:nb
%     i = sum( rand >= cumprob ) + 1;
%     res(k) = random('normal', mu(i), cov(i));
% end

end