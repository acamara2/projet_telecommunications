function z_bruit = bruitage(x,hr,Ns,M,RSB)
sigma_carre = (mean(x.^2)*Ns)/(2 * log2(M) * 10^(RSB/10));
bruit = randn(1,length(x))*sqrt(sigma_carre);
x_bruite = x + bruit;
z_bruit = conv(hr,x_bruite);
end

