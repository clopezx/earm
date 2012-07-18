function conc = earm2_embedded_init_conds()

conc = zeros(74, 1);
conc(1) = 3000; % L_0; L(bf=None)
conc(2) = 200; % R_0; R(bf=None)
conc(3) = 100; % flip_0; flip(bf=None)
conc(4) = 20000; % C8_0; C8(bf=None, state=pro)
conc(5) = 1000; % BAR_0; BAR(bf=None)
conc(6) = 100000; % Apaf_0; Apaf(bf=None, state=I)
conc(7) = 10000; % C3_0; C3(bf=None, state=pro)
conc(8) = 10000; % C6_0; C6(bf=None, state=pro)
conc(9) = 100000; % C9_0; C9(bf=None)
conc(10) = 1e+06; % PARP_0; PARP(bf=None, state=U)
conc(11) = 100000; % XIAP_0; XIAP(bf=None)
conc(12) = 40000; % Bid_0; Bid(bf=None, state=U)
conc(13) = 80000; % Bax_0; Bax(bf=None, s1=None, s2=None, state=C)
conc(14) = 20000; % Bak_0; Bak(bf=None, s1=None, s2=None, state=M)
conc(15) = 20000; % Bcl2_0; Bcl2(bf=None)
conc(16) = 20000; % BclxL_0; BclxL(bf=None, state=C)
conc(17) = 20000; % Mcl1_0; Mcl1(bf=None, state=M)
conc(18) = 500000; % CytoC_0; CytoC(bf=None, state=M)
conc(19) = 100000; % Smac_0; Smac(bf=None, state=M)

end

