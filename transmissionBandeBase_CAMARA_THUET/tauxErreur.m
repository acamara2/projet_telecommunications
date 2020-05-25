function taux = tauxErreur(attendu,trouve)
nb_erreurs = length(find(attendu ~= trouve));
taux = nb_erreurs/length(attendu);
end

