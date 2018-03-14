function Mk = genMk(k)
	Mk = zeros(2^k,2^k);
	for i = 1:2^k
		for j = 1:2^k
			p = bitand(i-1,bin2gray(j-1));
			Mk(i,j) = (-1)^length(strfind(dec2bin(p),'1'));
		end
	end	
end
