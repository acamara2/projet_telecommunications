function x_bruite = bruitage(x, Ns, M, RSB, equivalent)
    sigma_carre = (mean(abs(x).^2)*Ns)/(2 * log2(M) * 10^(RSB/10));
    bruit1 = randn(1,length(x))*sqrt(sigma_carre);
    if equivalent
        bruit2 = randn(1,length(x))*sqrt(sigma_carre);
        x_bruite = x + bruit1 + 1j * bruit2;
    else
        x_bruite = x + bruit1;
    end
end