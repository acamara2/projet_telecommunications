function decision = decision4aire(z_ech_bruit, Ns)
% Symboles decides
decisionMoins3 = (z_ech_bruit <= (-2*Ns))*(-3);
decisionMoins1 = ((-2*Ns<z_ech_bruit)&(z_ech_bruit<=0))*(-1);
decisionPlus1  = (0<z_ech_bruit)&(z_ech_bruit<=2*Ns);
decisionPlus3  =   (z_ech_bruit >= (2*Ns))*(3);
decision = decisionMoins3+decisionMoins1+decisionPlus1+decisionPlus3;
end

