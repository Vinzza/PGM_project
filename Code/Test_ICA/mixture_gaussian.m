function res = mixture_gaussian( nb, prob, mu, cov )
% Sample the given mixture of gaussian

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


end
