function res = mixture_gaussian_plot(x, prob, mu, cov )
% Sample the pdf of the given mixture of gaussian

  prob = prob ./ sum(prob);

  res = zeros(size(x));

  for i = 1:length(prob)
    res = res + prob(i) * normpdf( x, mu(i), cov(i) );
  end

end
