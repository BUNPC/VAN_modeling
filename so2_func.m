function so2 = so2_func( po2, species )
if ~exist('species','var')
    species=1;
end
if species==1 %human
    if 1
        % taken from Lobdell Am Physiol Soc 1981
        % table 2 and eq 1
        a = 0.34332;
        b = 0.64073;
        c = 0.34128;
        n = 1.58678;

        P50 = 26.6; % Roughton 1973

        %x = (po2/P50)*10^(0.024*(37-T) + 0.40*(pH-7.4) + 0.06 * log(40/pco2));
        x = po2/P50;

        so2 = (a*x.^n + b*x.^(2*n))./(1 + c*x.^n + b*x.^(2*n));

    elseif 0
        % Empirical model relating SO2 and pO2
        % Ref: http://www.ventworld.com/resources/oxydisso/dissoc.html#model
        a1 = -8.532e3;
        a2 =  2.121e3;
        a3 = -6.707e1;
        a4 =  9.360e5;
        a5 = -3.135e4;
        a6 =  2.396e3;
        a7 = -6.710e1;

        so2 = (a1*po2+a2*po2.^2+a3*po2.^3+po2.^4) ./ (a4+a5*po2+a6*po2.^2+a7*po2.^3+po2.^4);

    elseif 0
        % another source
        so2=(((po2.^3+150*po2).^-1*23400)+1).^-1;
    end
    
elseif species==2 %rat
    hillC = 2.7; p50 = 37; % Rat, Ellis C.G. et al., Am. J. Phys Heart Circ Phys 258, H1216-, 2002
    so2=po2.^hillC ./ ( po2.^hillC + p50.^hillC  );
    
elseif species==3 %mouse
    hillC = 2.59; p50 = 40.2; % C57BL/6 mice, Uchida K. et al., Zoological Science 15, 703-706, 1997
    so2=po2.^hillC ./ ( po2.^hillC + p50.^hillC  );
    
end