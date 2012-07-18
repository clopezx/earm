function out = earm2_embedded_odes(t, input, param)

% input(1) = L(bf=None);
% input(2) = R(bf=None);
% input(3) = flip(bf=None);
% input(4) = C8(bf=None, state=pro);
% input(5) = BAR(bf=None);
% input(6) = Apaf(bf=None, state=I);
% input(7) = C3(bf=None, state=pro);
% input(8) = C6(bf=None, state=pro);
% input(9) = C9(bf=None);
% input(10) = PARP(bf=None, state=U);
% input(11) = XIAP(bf=None);
% input(12) = Bid(bf=None, state=U);
% input(13) = Bax(bf=None, input(2)=None, input(3)=None, state=C);
% input(14) = Bak(bf=None, input(2)=None, input(3)=None, state=M);
% input(15) = Bcl2(bf=None);
% input(16) = BclxL(bf=None, state=C);
% input(17) = Mcl1(bf=None, state=M);
% input(18) = CytoC(bf=None, state=M);
% input(19) = Smac(bf=None, state=M);
% input(20) = L(bf=1) % R(bf=1);
% input(21) = Bax(bf=None, input(2)=None, input(3)=None, state=M);
% input(22) = BclxL(bf=None, state=M);
% input(23) = DISC(bf=None);
% input(24) = C8(bf=1, state=pro) % DISC(bf=1);
% input(25) = DISC(bf=1) % flip(bf=1);
% input(26) = C8(bf=None, state=A);
% input(27) = Bid(bf=1, state=U) % C8(bf=1, state=A);
% input(28) = BAR(bf=1) % C8(bf=1, state=A);
% input(29) = C3(bf=1, state=pro) % C8(bf=1, state=A);
% input(30) = Bid(bf=None, state=T);
% input(31) = C3(bf=None, state=A);
% input(32) = C3(bf=1, state=A) % XIAP(bf=1);
% input(33) = C3(bf=1, state=A) % PARP(bf=1, state=U);
% input(34) = C3(bf=1, state=A) % C6(bf=1, state=pro);
% input(35) = Bid(bf=None, state=M);
% input(36) = C3(bf=None, state=ub);
% input(37) = PARP(bf=None, state=C);
% input(38) = C6(bf=None, state=A);
% input(39) = Bax(bf=1, input(2)=None, input(3)=None, state=M) % Bid(bf=1, state=M);
% input(40) = Bak(bf=1, input(2)=None, input(3)=None, state=M) % Bid(bf=1, state=M);
% input(41) = Bcl2(bf=1) % Bid(bf=1, state=M);
% input(42) = BclxL(bf=1, state=M) % Bid(bf=1, state=M);
% input(43) = Bid(bf=1, state=M) % Mcl1(bf=1, state=M);
% input(44) = BclxL(bf=1, state=C) % Bid(bf=1, state=M);
% input(45) = C6(bf=1, state=A) % C8(bf=1, state=pro);
% input(46) = Bax(bf=None, input(2)=None, input(3)=None, state=A);
% input(47) = Bak(bf=None, input(2)=None, input(3)=None, state=A);
% input(48) = Bax(bf=1, input(2)=None, input(3)=None, state=A) % Bax(bf=1, input(2)=None, input(3)=None, state=C);
% input(49) = Bak(bf=1, input(2)=None, input(3)=None, state=A) % Bak(bf=1, input(2)=None, input(3)=None, state=M);
% input(50) = Bax(bf=1, input(2)=None, input(3)=None, state=A) % BclxL(bf=1, state=C);
% input(51) = Bax(bf=1, input(2)=None, input(3)=None, state=A) % Bcl2(bf=1);
% input(52) = Bax(bf=1, input(2)=None, input(3)=None, state=A) % BclxL(bf=1, state=M);
% input(53) = Bak(bf=1, input(2)=None, input(3)=None, state=A) % BclxL(bf=1, state=M);
% input(54) = Bak(bf=1, input(2)=None, input(3)=None, state=A) % Mcl1(bf=1, state=M);
% input(55) = Bax(bf=None, input(2)=1, input(3)=None, state=A) % Bax(bf=None, input(2)=None, input(3)=1, state=A);
% input(56) = Bak(bf=None, input(2)=1, input(3)=None, state=A) % Bak(bf=None, input(2)=None, input(3)=1, state=A);
% input(57) = Bax(bf=None, input(2)=1, input(3)=2, state=A) % Bax(bf=None, input(2)=2, input(3)=3, state=A) % Bax(bf=None, input(2)=3, input(3)=1, state=A);
% input(58) = Bak(bf=None, input(2)=1, input(3)=2, state=A) % Bak(bf=None, input(2)=2, input(3)=3, state=A) % Bak(bf=None, input(2)=3, input(3)=1, state=A);
% input(59) = Bax(bf=None, input(2)=1, input(3)=2, state=A) % Bax(bf=None, input(2)=2, input(3)=3, state=A) % Bax(bf=None, input(2)=3, input(3)=4, state=A) % Bax(bf=None, input(2)=4, input(3)=1, state=A);
% input(60) = Bak(bf=None, input(2)=1, input(3)=2, state=A) % Bak(bf=None, input(2)=2, input(3)=3, state=A) % Bak(bf=None, input(2)=3, input(3)=4, state=A) % Bak(bf=None, input(2)=4, input(3)=1, state=A);
% input(61) = Bax(bf=1, input(2)=2, input(3)=3, state=A) % Bax(bf=None, input(2)=3, input(3)=4, state=A) % Bax(bf=None, input(2)=4, input(3)=5, state=A) % Bax(bf=None, input(2)=5, input(3)=2, state=A) % CytoC(bf=1, state=M);
% input(62) = Bax(bf=1, input(2)=2, input(3)=3, state=A) % Bax(bf=None, input(2)=3, input(3)=4, state=A) % Bax(bf=None, input(2)=4, input(3)=5, state=A) % Bax(bf=None, input(2)=5, input(3)=2, state=A) % Smac(bf=1, state=M);
% input(63) = Bak(bf=1, input(2)=2, input(3)=3, state=A) % Bak(bf=None, input(2)=3, input(3)=4, state=A) % Bak(bf=None, input(2)=4, input(3)=5, state=A) % Bak(bf=None, input(2)=5, input(3)=2, state=A) % CytoC(bf=1, state=M);
% input(64) = Bak(bf=1, input(2)=2, input(3)=3, state=A) % Bak(bf=None, input(2)=3, input(3)=4, state=A) % Bak(bf=None, input(2)=4, input(3)=5, state=A) % Bak(bf=None, input(2)=5, input(3)=2, state=A) % Smac(bf=1, state=M);
% input(65) = CytoC(bf=None, state=C);
% input(66) = Smac(bf=None, state=C);
% input(67) = Smac(bf=None, state=A);
% input(68) = CytoC(bf=None, state=A);
% input(69) = Apaf(bf=1, state=I) % CytoC(bf=1, state=A);
% input(70) = Smac(bf=1, state=A) % XIAP(bf=1);
% input(71) = Apaf(bf=None, state=A);
% input(72) = Apop(bf=None);
% input(73) = Apop(bf=1) % C3(bf=1, state=pro);
% input(74) = Apop(bf=1) % XIAP(bf=1);

% L(bf=None)
out(1,1) = -param(6)*input(1)*input(2) + param(7)*input(20);
% R(bf=None)
out(2,1) = -param(6)*input(1)*input(2) + param(7)*input(20);
% flip(bf=None)
out(3,1) = -param(15)*input(3)*input(23) + param(16)*input(25);
% C8(bf=None, state=pro)
out(4,1) = -param(53)*input(4)*input(38) + param(54)*input(45) - param(9)*input(23)*input(4) + param(10)*input(24);
% BAR(bf=None)
out(5,1) = -param(17)*input(26)*input(5) + param(18)*input(28);
% Apaf(bf=None, state=I)
out(6,1) = -param(29)*input(6)*input(68) + param(30)*input(69);
% C3(bf=None, state=pro)
out(7,1) = -param(34)*input(7)*input(72) + param(35)*input(73) - param(41)*input(26)*input(7) + param(42)*input(29);
% C6(bf=None, state=pro)
out(8,1) = -param(50)*input(31)*input(8) + param(51)*input(34);
% C9(bf=None)
out(9,1) = -param(32)*input(71)*input(9) + param(33)*input(72);
% PARP(bf=None, state=U)
out(10,1) = -param(47)*input(31)*input(10) + param(48)*input(33);
% XIAP(bf=None)
out(11,1) = -param(37)*input(11)*input(72) + param(38)*input(74) - param(39)*input(11)*input(67) + param(40)*input(70) - param(44)*input(11)*input(31) + param(45)*input(32) + param(46)*input(32);
% Bid(bf=None, state=U)
out(12,1) = -param(12)*input(12)*input(26) + param(13)*input(27);
% Bax(bf=None, input(2)=None, input(3)=None, state=C)
out(13,1) = -param(82)*input(13)*input(46) + param(83)*input(48) - param(66)*input(13) + param(67)*input(21);
% Bak(bf=None, input(2)=None, input(3)=None, state=M)
out(14,1) = -param(85)*input(14)*input(47) + param(86)*input(49) - param(73)*input(14)*input(35) + param(74)*input(40);
% Bcl2(bf=None)
out(15,1) = -param(94)*input(15)*input(46) + param(95)*input(51) - param(76)*input(15)*input(35) + param(77)*input(41);
% BclxL(bf=None, state=C)
out(16,1) = -param(91)*input(16)*input(46) + param(92)*input(50) - param(88)*input(16)*input(35) + param(89)*input(44) - param(68)*input(16) + param(69)*input(22);
% Mcl1(bf=None, state=M)
out(17,1) = -param(100)*input(17)*input(47) + param(101)*input(54) - param(80)*input(17)*input(35) + param(81)*input(43);
% CytoC(bf=None, state=M)
out(18,1) = -param(120)*input(18)*input(60) + param(121)*input(63) - param(114)*input(18)*input(59) + param(115)*input(61);
% Smac(bf=None, state=M)
out(19,1) = -param(123)*input(19)*input(60) + param(124)*input(64) - param(117)*input(19)*input(59) + param(118)*input(62);
% L(bf=1) % R(bf=1)
out(20,1) = param(6)*input(1)*input(2) - param(7)*input(20) - param(8)*input(20);
% Bax(bf=None, input(2)=None, input(3)=None, state=M)
out(21,1) = -param(70)*input(21)*input(35) + param(71)*input(39) + param(66)*input(13) - param(67)*input(21);
% BclxL(bf=None, state=M)
out(22,1) = -param(98)*input(22)*input(47) + param(99)*input(53) - param(96)*input(22)*input(46) + param(97)*input(52) - param(78)*input(22)*input(35) + param(79)*input(42) + param(93)*input(50) + param(90)*input(44) + param(68)*input(16) - param(69)*input(22);
% DISC(bf=None)
out(23,1) = -param(9)*input(23)*input(4) + param(10)*input(24) - param(15)*input(3)*input(23) + param(16)*input(25) + param(11)*input(24) + param(8)*input(20);
% C8(bf=1, state=pro) % DISC(bf=1)
out(24,1) = param(9)*input(23)*input(4) - param(10)*input(24) - param(11)*input(24);
% DISC(bf=1) % flip(bf=1)
out(25,1) = param(15)*input(3)*input(23) - param(16)*input(25);
% C8(bf=None, state=A)
out(26,1) = -param(17)*input(26)*input(5) + param(18)*input(28) - param(12)*input(12)*input(26) + param(13)*input(27) - param(41)*input(26)*input(7) + param(42)*input(29) + param(55)*input(45) + param(14)*input(27) + param(43)*input(29) + param(11)*input(24);
% Bid(bf=1, state=U) % C8(bf=1, state=A)
out(27,1) = param(12)*input(12)*input(26) - param(13)*input(27) - param(14)*input(27);
% BAR(bf=1) % C8(bf=1, state=A)
out(28,1) = param(17)*input(26)*input(5) - param(18)*input(28);
% C3(bf=1, state=pro) % C8(bf=1, state=A)
out(29,1) = param(41)*input(26)*input(7) - param(42)*input(29) - param(43)*input(29);
% Bid(bf=None, state=T)
out(30,1) = param(14)*input(27) - param(64)*input(30) + param(65)*input(35);
% C3(bf=None, state=A)
out(31,1) = -param(50)*input(31)*input(8) + param(51)*input(34) - param(47)*input(31)*input(10) + param(48)*input(33) - param(44)*input(11)*input(31) + param(45)*input(32) + param(36)*input(73) + param(52)*input(34) + param(49)*input(33) + param(43)*input(29);
% C3(bf=1, state=A) % XIAP(bf=1)
out(32,1) = param(44)*input(11)*input(31) - param(45)*input(32) - param(46)*input(32);
% C3(bf=1, state=A) % PARP(bf=1, state=U)
out(33,1) = param(47)*input(31)*input(10) - param(48)*input(33) - param(49)*input(33);
% C3(bf=1, state=A) % C6(bf=1, state=pro)
out(34,1) = param(50)*input(31)*input(8) - param(51)*input(34) - param(52)*input(34);
% Bid(bf=None, state=M)
out(35,1) = -param(73)*input(14)*input(35) + param(74)*input(40) - param(70)*input(21)*input(35) + param(71)*input(39) - param(76)*input(15)*input(35) + param(77)*input(41) - param(88)*input(16)*input(35) + param(89)*input(44) - param(78)*input(22)*input(35) + param(79)*input(42) - param(80)*input(17)*input(35) + param(81)*input(43) + param(75)*input(40) + param(72)*input(39) + param(90)*input(44) + param(64)*input(30) - param(65)*input(35);
% C3(bf=None, state=ub)
out(36,1) = param(46)*input(32);
% PARP(bf=None, state=C)
out(37,1) = param(49)*input(33);
% C6(bf=None, state=A)
out(38,1) = -param(53)*input(4)*input(38) + param(54)*input(45) + param(52)*input(34) + param(55)*input(45);
% Bax(bf=1, input(2)=None, input(3)=None, state=M) % Bid(bf=1, state=M)
out(39,1) = param(70)*input(21)*input(35) - param(71)*input(39) - param(72)*input(39);
% Bak(bf=1, input(2)=None, input(3)=None, state=M) % Bid(bf=1, state=M)
out(40,1) = param(73)*input(14)*input(35) - param(74)*input(40) - param(75)*input(40);
% Bcl2(bf=1) % Bid(bf=1, state=M)
out(41,1) = param(76)*input(15)*input(35) - param(77)*input(41);
% BclxL(bf=1, state=M) % Bid(bf=1, state=M)
out(42,1) = param(78)*input(22)*input(35) - param(79)*input(42);
% Bid(bf=1, state=M) % Mcl1(bf=1, state=M)
out(43,1) = param(80)*input(17)*input(35) - param(81)*input(43);
% BclxL(bf=1, state=C) % Bid(bf=1, state=M)
out(44,1) = param(88)*input(16)*input(35) - param(89)*input(44) - param(90)*input(44);
% C6(bf=1, state=A) % C8(bf=1, state=pro)
out(45,1) = param(53)*input(4)*input(38) - param(54)*input(45) - param(55)*input(45);
% Bax(bf=None, input(2)=None, input(3)=None, state=A)
out(46,1) = -2*param(102)*power(input(46), 2) + 2*param(103)*input(55) - param(104)*input(46)*input(55) + param(105)*input(57) - param(106)*input(46)*input(57) + param(107)*input(59) - param(82)*input(13)*input(46) + param(83)*input(48) - param(94)*input(15)*input(46) + param(95)*input(51) - param(91)*input(16)*input(46) + param(92)*input(50) - param(96)*input(22)*input(46) + param(97)*input(52) + 2*param(84)*input(48) + param(93)*input(50) + param(72)*input(39);
% Bak(bf=None, input(2)=None, input(3)=None, state=A)
out(47,1) = -2*param(108)*power(input(47), 2) + 2*param(109)*input(56) - param(110)*input(47)*input(56) + param(111)*input(58) - param(112)*input(47)*input(58) + param(113)*input(60) - param(85)*input(14)*input(47) + param(86)*input(49) - param(98)*input(22)*input(47) + param(99)*input(53) - param(100)*input(17)*input(47) + param(101)*input(54) + 2*param(87)*input(49) + param(75)*input(40);
% Bax(bf=1, input(2)=None, input(3)=None, state=A) % Bax(bf=1, input(2)=None, input(3)=None, state=C)
out(48,1) = param(82)*input(13)*input(46) - param(83)*input(48) - param(84)*input(48);
% Bak(bf=1, input(2)=None, input(3)=None, state=A) % Bak(bf=1, input(2)=None, input(3)=None, state=M)
out(49,1) = param(85)*input(14)*input(47) - param(86)*input(49) - param(87)*input(49);
% Bax(bf=1, input(2)=None, input(3)=None, state=A) % BclxL(bf=1, state=C)
out(50,1) = param(91)*input(16)*input(46) - param(92)*input(50) - param(93)*input(50);
% Bax(bf=1, input(2)=None, input(3)=None, state=A) % Bcl2(bf=1)
out(51,1) = param(94)*input(15)*input(46) - param(95)*input(51);
% Bax(bf=1, input(2)=None, input(3)=None, state=A) % BclxL(bf=1, state=M)
out(52,1) = param(96)*input(22)*input(46) - param(97)*input(52);
% Bak(bf=1, input(2)=None, input(3)=None, state=A) % BclxL(bf=1, state=M)
out(53,1) = param(98)*input(22)*input(47) - param(99)*input(53);
% Bak(bf=1, input(2)=None, input(3)=None, state=A) % Mcl1(bf=1, state=M)
out(54,1) = param(100)*input(17)*input(47) - param(101)*input(54);
% Bax(bf=None, input(2)=1, input(3)=None, state=A) % Bax(bf=None, input(2)=None, input(3)=1, state=A)
out(55,1) = param(102)*power(input(46), 2) - param(103)*input(55) - param(104)*input(46)*input(55) + param(105)*input(57);
% Bak(bf=None, input(2)=1, input(3)=None, state=A) % Bak(bf=None, input(2)=None, input(3)=1, state=A)
out(56,1) = param(108)*power(input(47), 2) - param(109)*input(56) - param(110)*input(47)*input(56) + param(111)*input(58);
% Bax(bf=None, input(2)=1, input(3)=2, state=A) % Bax(bf=None, input(2)=2, input(3)=3, state=A) % Bax(bf=None, input(2)=3, input(3)=1, state=A)
out(57,1) = param(104)*input(46)*input(55) - param(105)*input(57) - param(106)*input(46)*input(57) + param(107)*input(59);
% Bak(bf=None, input(2)=1, input(3)=2, state=A) % Bak(bf=None, input(2)=2, input(3)=3, state=A) % Bak(bf=None, input(2)=3, input(3)=1, state=A)
out(58,1) = param(110)*input(47)*input(56) - param(111)*input(58) - param(112)*input(47)*input(58) + param(113)*input(60);
% Bax(bf=None, input(2)=1, input(3)=2, state=A) % Bax(bf=None, input(2)=2, input(3)=3, state=A) % Bax(bf=None, input(2)=3, input(3)=4, state=A) % Bax(bf=None, input(2)=4, input(3)=1, state=A)
out(59,1) = param(106)*input(46)*input(57) - param(107)*input(59) - param(114)*input(18)*input(59) + param(115)*input(61) - param(117)*input(19)*input(59) + param(118)*input(62) + param(116)*input(61) + param(119)*input(62);
% Bak(bf=None, input(2)=1, input(3)=2, state=A) % Bak(bf=None, input(2)=2, input(3)=3, state=A) % Bak(bf=None, input(2)=3, input(3)=4, state=A) % Bak(bf=None, input(2)=4, input(3)=1, state=A)
out(60,1) = param(112)*input(47)*input(58) - param(113)*input(60) - param(120)*input(18)*input(60) + param(121)*input(63) - param(123)*input(19)*input(60) + param(124)*input(64) + param(122)*input(63) + param(125)*input(64);
% Bax(bf=1, input(2)=2, input(3)=3, state=A) % Bax(bf=None, input(2)=3, input(3)=4, state=A) % Bax(bf=None, input(2)=4, input(3)=5, state=A) % Bax(bf=None, input(2)=5, input(3)=2, state=A) % CytoC(bf=1, state=M)
out(61,1) = param(114)*input(18)*input(59) - param(115)*input(61) - param(116)*input(61);
% Bax(bf=1, input(2)=2, input(3)=3, state=A) % Bax(bf=None, input(2)=3, input(3)=4, state=A) % Bax(bf=None, input(2)=4, input(3)=5, state=A) % Bax(bf=None, input(2)=5, input(3)=2, state=A) % Smac(bf=1, state=M)
out(62,1) = param(117)*input(19)*input(59) - param(118)*input(62) - param(119)*input(62);
% Bak(bf=1, input(2)=2, input(3)=3, state=A) % Bak(bf=None, input(2)=3, input(3)=4, state=A) % Bak(bf=None, input(2)=4, input(3)=5, state=A) % Bak(bf=None, input(2)=5, input(3)=2, state=A) % CytoC(bf=1, state=M)
out(63,1) = param(120)*input(18)*input(60) - param(121)*input(63) - param(122)*input(63);
% Bak(bf=1, input(2)=2, input(3)=3, state=A) % Bak(bf=None, input(2)=3, input(3)=4, state=A) % Bak(bf=None, input(2)=4, input(3)=5, state=A) % Bak(bf=None, input(2)=5, input(3)=2, state=A) % Smac(bf=1, state=M)
out(64,1) = param(123)*input(19)*input(60) - param(124)*input(64) - param(125)*input(64);
% CytoC(bf=None, state=C)
out(65,1) = -param(27)*input(65) + param(28)*input(68) + param(122)*input(63) + param(116)*input(61);
% Smac(bf=None, state=C)
out(66,1) = -param(25)*input(66) + param(26)*input(67) + param(125)*input(64) + param(119)*input(62);
% Smac(bf=None, state=A)
out(67,1) = -param(39)*input(11)*input(67) + param(40)*input(70) + param(25)*input(66) - param(26)*input(67);
% CytoC(bf=None, state=A)
out(68,1) = -param(29)*input(6)*input(68) + param(30)*input(69) + param(31)*input(69) + param(27)*input(65) - param(28)*input(68);
% Apaf(bf=1, state=I) % CytoC(bf=1, state=A)
out(69,1) = param(29)*input(6)*input(68) - param(30)*input(69) - param(31)*input(69);
% Smac(bf=1, state=A) % XIAP(bf=1)
out(70,1) = param(39)*input(11)*input(67) - param(40)*input(70);
% Apaf(bf=None, state=A)
out(71,1) = param(31)*input(69) - param(32)*input(71)*input(9) + param(33)*input(72);
% Apop(bf=None)
out(72,1) = -param(34)*input(7)*input(72) + param(35)*input(73) - param(37)*input(11)*input(72) + param(38)*input(74) + param(36)*input(73) + param(32)*input(71)*input(9) - param(33)*input(72);
% Apop(bf=1) % C3(bf=1, state=pro)
out(73,1) = param(34)*input(7)*input(72) - param(35)*input(73) - param(36)*input(73);
% Apop(bf=1) % XIAP(bf=1)
out(74,1) = param(37)*input(11)*input(72) - param(38)*input(74);

end
