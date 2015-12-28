function res = mixture_gaussian( nb, prob, mu, cov )

prob = prob ./ sum(prob);
cumprob = cumsum(prob(:));
res = zeros(1,nb);

for k = 1:nb
    i = sum( rand >= cumprob ) + 1;
    res(k) = random('normal', mu(i), cov(i));
end

end