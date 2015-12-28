function d = dBach(V,W)
% Computes the Bach distance between arrays U and V)
A = V/W;
[m,~] = size(V);
d1 = sum(sum(abs(A),2)./max(abs(A'))'-1);
d2 = sum(sum(abs(A),1)./max(abs(A))-1);
d = (d1+d2)/(2*m);
end
