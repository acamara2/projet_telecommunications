function bits_decides_bruit = decision2aire(z_ech_bruit)
% Symboles decides
symboles_decides_bruit=sign(z_ech_bruit);
% Demapping
bits_decides_bruit = (1/2) *( symboles_decides_bruit + 1);
end

