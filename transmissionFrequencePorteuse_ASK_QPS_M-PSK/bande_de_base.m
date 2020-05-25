    function zb = bande_de_base(x, t, fp, fe, N)
    %   Demodulation
    Im = x.*cos(2*pi*fp*t);
    Jm = x.*sin(2*pi*fp*t);

    %   Filtrage passe-bas
    k = (-N:N)/fe;

    hb =  2*fp/fe*sinc(2*fp*k);

    Ib = conv(Im,hb,'same');
    Jb = -conv(Jm,hb,'same');
    
    zb = Ib + 1j*Jb;
end