function taux = TEB(bits_decides, bits, nb_bits)
    nb_erreurs = length(find(bits ~= bits_decides));
    taux = nb_erreurs/nb_bits;
end