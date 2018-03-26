q = 3;
U = eye(2^q);
line = 0;

fid = fopen('QSD_qasm.txt');
tline = fgetl(fid);
while ischar(tline)
%     disp(tline)
	ss = strsplit(tline,{' q',',q',','});
	if (size(ss{1},2) == 2) 
		if (ss{1} == 'ry')
			g = Ry(-str2num(ss{3}));
		elseif (ss{1} == 'rz')
			g = Rz(-str2num(ss{3}));
		else
% 			g = AP(-str2num(ss{3}));	% unmask for AP
			g = eye(2);
		end
		U = kron(eye(2^(q-1-str2num(ss{2}))),kron(g,eye(2^(str2num(ss{2})))))*U;
	else
		U = U_CX(str2num(ss{2}),str2num(ss{3}),q)*U;
	end
    tline = fgetl(fid);
end
fclose(fid);
U

%%

% 0         0    0.3536    0.3536    0.3536    0.3536    0.3536    0.6124
% 
% 
% 0.6857 - 0.2840i   0.3959 - 0.1640i   0.0693 - 0.0287i  -0.0981 + 0.0406i   0.0693 - 0.0287i  -0.0981 + 0.0406i  -0.0981 + 0.0406i  -0.4338 + 0.1797i
% -0.5627 + 0.2331i  -0.2729 + 0.1130i   0.0538 - 0.0223i   0.2211 - 0.0916i   0.0538 - 0.0223i   0.2211 - 0.0916i   0.2211 - 0.0916i   0.5568 - 0.2306i
% 
% 
% -0.5627 + 0.2331i  -0.2729 + 0.1130i  -0.0538 + 0.0223i  -0.2211 + 0.0916i  -0.0538 + 0.0223i  -0.2211 + 0.0916i  -0.2211 + 0.0916i  -0.5568 + 0.2306i
% 0.0219 - 0.0091i  -0.2679 + 0.1110i  -0.4870 + 0.2017i  -0.3197 + 0.1324i  -0.4870 + 0.2017i  -0.3197 + 0.1324i  -0.3197 + 0.1324i   0.0160 - 0.0066i
% 
% 
% -0.6701 + 0.6701i  -0.1722 + 0.1722i  -0.0045 + 0.0045i   0.0269 - 0.0269i  -0.0045 + 0.0045i   0.0269 - 0.0269i   0.0269 - 0.0269i  -0.1380 + 0.1380i
% 0.4430 - 0.4430i  -0.0550 + 0.0550i  -0.2227 + 0.2227i  -0.2540 + 0.2540i  -0.2227 + 0.2227i  -0.2540 + 0.2540i  -0.2540 + 0.2540i  -0.0892 + 0.0892i
% 
% -0.2932 + 0.7079i  -0.0011 + 0.0027i   0.0896 - 0.2164i   0.1197 - 0.2889i   0.0896 - 0.2164i   0.1197 - 0.2889i   0.1197 - 0.2889i   0.0380 - 0.0916i
% 0.3637 - 0.8780i   0.0716 - 0.1728i  -0.0192 + 0.0463i  -0.0492 + 0.1187i  -0.0192 + 0.0463i  -0.0492 + 0.1187i  -0.0492 + 0.1187i   0.0325 - 0.0785i
% 
% -0.0000 + 0.2008i  -0.0000 - 0.2457i   0.0000 - 0.4828i  -0.0000 - 0.3784i  -0.0000 - 0.4828i  -0.0000 - 0.3784i  -0.0000 - 0.3784i  -0.0000 - 0.0593i
% 0.0000 - 0.7520i   0.0000 - 0.3056i  -0.0000 - 0.0684i   0.0000 - 0.1729i   0.0000 - 0.0684i   0.0000 - 0.1729i  -0.0000 - 0.1729i   0.0000 - 0.4920i
% 
% -0.1239 - 0.2991i  -0.1208 - 0.2915i  -0.2115 - 0.5106i  -0.0711 - 0.1716i  -0.2115 - 0.5106i  -0.0711 - 0.1716i  -0.0711 - 0.1716i   0.1090 + 0.2633i
% -0.0691 - 0.1667i  -0.0722 - 0.1743i   0.0185 + 0.0448i  -0.1219 - 0.2943i   0.0185 + 0.0448i  -0.1219 - 0.2943i  -0.1219 - 0.2943i  -0.3020 - 0.7291i
% 
% -0.2923 - 0.2923i  -0.1090 - 0.1090i  -0.2767 - 0.2767i   0.0852 + 0.0852i  -0.2767 - 0.2767i   0.0852 + 0.0852i   0.0852 + 0.0852i   0.4772 + 0.4772i
% 0.2369 + 0.2369i   0.0536 + 0.0536i   0.2213 + 0.2213i  -0.1407 - 0.1407i   0.2213 + 0.2213i  -0.1407 - 0.1407i  -0.1407 - 0.1407i  -0.5327 - 0.5327i
% 
% -0.0007 - 0.0003i   0.1083 + 0.0448i  -0.1108 - 0.0459i   0.2867 + 0.1188i  -0.1108 - 0.0459i   0.2867 + 0.1188i   0.2867 + 0.1188i   0.7554 + 0.3129i
% 0.3761 + 0.1558i   0.2671 + 0.1106i   0.4862 + 0.2014i   0.0886 + 0.0367i   0.4862 + 0.2014i   0.0886 + 0.0367i   0.0886 + 0.0367i  -0.3800 - 0.1574i
   
% rz q0,-3.141593
% ry q0,-1.749068
% rz q0,-1.570796
% ap q0,-0.785398
% rz q1,1.570796
% cnot q0,q1
% rz q1,-1.570796
% cnot q0,q1
% rz q0,-3.141593
% ry q0,-0.504834
% rz q0,3.141593
% ap q0,-0.000000
% ry q1,1.733418
% cnot q0,q1
% ry q1,0.762379
% cnot q0,q1
% rz q0,-0.535601
% ry q0,-1.437175
% rz q0,-1.157558
% ap q0,-1.178097
% rz q1,-0.785398
% cnot q0,q1
% rz q1,0.960879
% cnot q0,q1
% rz q0,0.626540
% ry q0,-1.040477
% rz q0,-0.626540
% ap q0,0.000000
% rz q1,-1.570796
% cnot q0,q2
% rz q1,-0.021343
% cnot q1,q2
% rz q1,0.021343
% cnot q0,q2
% rz q1,-1.570796
% cnot q1,q2
% rz q0,2.285636
% ry q0,-2.217616
% rz q0,0.654520
% ap q0,-1.178097
% rz q1,-0.785398
% cnot q0,q1
% rz q1,-1.088849
% cnot q0,q1
% rz q0,2.140665
% ry q0,-0.797715
% rz q0,-2.140665
% ap q0,0.000000
% ry q1,0.536833
% cnot q0,q1
% ry q1,0.536833
% cnot q0,q1
% rz q0,4.712389
% ry q0,-1.570796
% rz q0,0.815214
% ap q0,-0.785398
% rz q1,-0.000000
% cnot q0,q1
% rz q1,-0.815214
% cnot q0,q1
% rz q0,-0.000000
% ry q0,-1.570796
% rz q0,-1.570796
% ap q0,0.785398
% ry q1,0.523599
% cnot q0,q2
% ry q1,0.523599
% cnot q1,q2
% ry q1,0.523599
% cnot q0,q2
% ry q1,0.523599
% cnot q1,q2
% rz q0,-3.141593
% ry q0,-2.326379
% rz q0,-1.570796
% ap q0,0.785398
% rz q1,1.570796
% cnot q0,q1
% rz q1,1.570796
% cnot q0,q1
% rz q0,0.000000
% ry q0,-0.815214
% rz q0,-0.000000
% ap q0,0.000000
% ry q1,0.536833
% cnot q0,q1
% ry q1,0.536833
% cnot q0,q1
% rz q0,-2.780856
% ry q0,-1.231483
% rz q0,-2.779386
% ap q0,-0.785398
% rz q1,-0.000000
% cnot q0,q1
% rz q1,-2.293740
% cnot q0,q1
% rz q0,1.931533
% ry q0,-1.231483
% rz q0,-1.931533
% ap q0,-0.000000
% rz q1,-0.785398
% cnot q0,q2
% rz q1,-0.806741
% cnot q1,q2
% rz q1,-0.764055
% cnot q0,q2
% rz q1,-0.785398
% cnot q1,q2
% rz q0,1.234435
% ry q0,-0.814217
% rz q0,3.055660
% ap q0,-0.392699
% rz q1,0.785398
% cnot q0,q1
% rz q1,1.806353
% cnot q0,q1
% rz q0,-2.208546
% ry q0,-1.024349
% rz q0,2.208546
% ap q0,-0.000000
% ry q1,1.733418
% cnot q0,q1
% ry q1,0.762379
% cnot q0,q1
% rz q0,3.141593
% ry q0,-2.636759
% rz q0,1.570796
% ap q0,0.785398
% rz q1,-1.570796
% cnot q0,q1
% rz q1,1.570796
% cnot q0,q1
% rz q0,3.141593
% ry q0,-1.392524
% rz q0,-3.141593
% ap q0,-0.000000
% 
% BO =
%     0.1563   -0.4871   -0.4871   -0.2813   -0.4871   -0.2813   -0.2813   -0.1624
%    -0.4871    0.7188   -0.2813   -0.1624   -0.2813   -0.1624   -0.1624   -0.0938
%    -0.4871   -0.2813    0.7188   -0.1624   -0.2813   -0.1624   -0.1624   -0.0938
%    -0.2813   -0.1624   -0.1624    0.9063   -0.1624   -0.0937   -0.0937   -0.0541
%    -0.4871   -0.2813   -0.2813   -0.1624    0.7188   -0.1624   -0.1624   -0.0938
%    -0.2813   -0.1624   -0.1624   -0.0937   -0.1624    0.9063   -0.0937   -0.0541
%    -0.2813   -0.1624   -0.1624   -0.0937   -0.1624   -0.0937    0.9063   -0.0541
%    -0.1624   -0.0938   -0.0938   -0.0541   -0.0938   -0.0541   -0.0541    0.9688
% 
% BOD =														
%   -0.1444 - 0.0598i   0.4501 + 0.1864i   0.4501 + 0.1864i   0.2598 + 0.1076i   0.4501 + 0.1864i   0.2598 + 0.1076i   0.2598 + 0.1076i   0.1500 + 0.0621i
%    0.4501 + 0.1864i  -0.6640 - 0.2751i   0.2598 + 0.1076i   0.1500 + 0.0621i   0.2598 + 0.1076i   0.1500 + 0.0621i   0.1500 + 0.0621i   0.0866 + 0.0359i
%    0.4501 + 0.1864i   0.2598 + 0.1076i  -0.6640 - 0.2751i   0.1500 + 0.0621i   0.2598 + 0.1076i   0.1500 + 0.0621i   0.1500 + 0.0621i   0.0866 + 0.0359i
%    0.2598 + 0.1076i   0.1500 + 0.0621i   0.1500 + 0.0621i  -0.8373 - 0.3468i   0.1500 + 0.0621i   0.0866 + 0.0359i   0.0866 + 0.0359i   0.0500 + 0.0207i
%    0.4501 + 0.1864i   0.2598 + 0.1076i   0.2598 + 0.1076i   0.1500 + 0.0621i  -0.6640 - 0.2751i   0.1500 + 0.0621i   0.1500 + 0.0621i   0.0866 + 0.0359i
%    0.2598 + 0.1076i   0.1500 + 0.0621i   0.1500 + 0.0621i   0.0866 + 0.0359i   0.1500 + 0.0621i  -0.8373 - 0.3468i   0.0866 + 0.0359i   0.0500 + 0.0207i
%    0.2598 + 0.1076i   0.1500 + 0.0621i   0.1500 + 0.0621i   0.0866 + 0.0359i   0.1500 + 0.0621i   0.0866 + 0.0359i  -0.8373 - 0.3468i   0.0500 + 0.0207i
%    0.1500 + 0.0621i   0.0866 + 0.0359i   0.0866 + 0.0359i   0.0500 + 0.0207i   0.0866 + 0.0359i   0.0500 + 0.0207i   0.0500 + 0.0207i   -0.8950 - 0.3707i

%% QuInE qasm
% 
% ry q0,-0.927295
% rz q0,-6.2832
% rz q1,-3.141593
% 
% ry q1,0.927295
% cnot q0,q1
% ry q1,0.927295
% cnot q0,q1
% 
% rz q0,3.141593
% rz q0,1.570796
% rz q1,-1.570796
% cnot q0,q1
% rz q1,1.570796
% cnot q0,q1
% ry q0,-0.927295
% rz q0,-3.141593
% 
% U =
%   -0.3795 - 0.3795i  -0.3162 - 0.3162i  -0.5060 - 0.5060i   0.0000 + 0.0000i
%   -0.1897 - 0.1897i   0.6325 + 0.6325i  -0.2530 - 0.2530i   0.0000 + 0.0000i
%   -0.5060 - 0.5060i   0.0000 + 0.0000i   0.3795 + 0.3795i  -0.3162 - 0.3162i
%   -0.2530 - 0.2530i   0.0000 + 0.0000i   0.1897 + 0.1897i   0.6325 + 0.6325i
%   
% -0.3795 + 0.3795i   0.3163 - 0.3163i   0.0000 - 0.0000i   0.5060 - 0.5060i
% -0.0000 + 0.6761i   0.0000 + 0.2201i  -0.0000 + 0.7000i  -0.0000 - 0.0679i
% -0.1890 - 0.1890i   0.3252 + 0.3252i   0.4364 + 0.4364i   0.4101 + 0.4101i
% 
% 0.2880    0.2000    0.0000    0.5120
% 0.4571    0.0484    0.4900    0.0046
% 0.0715    0.2115    0.3809    0.3364

%% Now QuAMdq
% 
% rz q0,-0.000000
% ry q0,-0.927295
% rz q0,-1.570796
% ap q0,-0.785398
% rz q1,-1.570796
% cnot q0,q1
% rz q1,-1.570796
% cnot q0,q1
% ry q1,0.927295
% cnot q0,q1
% ry q1,0.927295
% cnot q0,q1
% rz q0,3.141593
% ry q0,-0.927295
% rz q0,-3.141593
% ap q0,-0.000000
%
% BO =
%    -0.2800   -0.6400   -0.6400   -0.3200
%    -0.6400    0.6800   -0.3200   -0.1600
%    -0.6400   -0.3200    0.6800   -0.1600
%    -0.3200   -0.1600   -0.1600    0.9200
%    
% BOD =
%   -0.2800 + 0.0000i  -0.6400 + 0.0000i  -0.6400 - 0.0000i  -0.3200 - 0.0000i
%   -0.6400 + 0.0000i   0.6800 - 0.0000i  -0.3200 - 0.0000i  -0.1600 - 0.0000i
%   -0.6400 + 0.0000i  -0.3200 + 0.0000i   0.6800 + 0.0000i  -0.1600 + 0.0000i
%   -0.3200 + 0.0000i  -0.1600 + 0.0000i  -0.1600 + 0.0000i   0.9200 - 0.0000i
% 
% BOD =
%   -0.1980 + 0.1980i  -0.4525 + 0.4525i  -0.4525 + 0.4525i  -0.2263 + 0.2263i
%   -0.4525 + 0.4525i   0.4808 - 0.4808i  -0.2263 + 0.2263i  -0.1131 + 0.1131i
%   -0.4525 + 0.4525i  -0.2263 + 0.2263i   0.4808 - 0.4808i  -0.1131 + 0.1131i
%   -0.2263 + 0.2263i  -0.1131 + 0.1131i  -0.1131 + 0.1131i   0.6505 - 0.6505i
% 	
% U recomposition matches with both AP and non-AP runs of BOD
% 
% state with AP
% -0.3200 + 0.0000i   0.4000 + 0.0000i   0.4000 - 0.0000i   0.7600 + 0.0000i
% -0.8968 + 0.0000i  -0.0040 + 0.0000i  -0.0040 - 0.0000i   0.4424 + 0.0000i
% -0.7920 - 0.0000i  -0.4050 + 0.0000i  -0.4050 - 0.0000i  -0.2114 - 0.0000i
% 
% state without AP
% -0.2263 - 0.2263i   0.2828 + 0.2828i   0.2828 + 0.2828i   0.5374 + 0.5374i
% -0.0000 - 0.8968i  -0.0000 - 0.0040i   0.0000 - 0.0040i   0.0000 + 0.4424i
% 0.5601 - 0.5601i   0.2863 - 0.2863i   0.2863 - 0.2863i   0.1495 - 0.1495i
%
% state probabilities for both BOD
% 0.1024    0.1600    0.1600    0.5776
% 0.8043    0.0000    0.0000    0.1957
% 0.6273    0.1640    0.1640    0.0447