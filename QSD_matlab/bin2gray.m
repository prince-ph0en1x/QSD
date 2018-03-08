function num = bin2gray(num)
	num = bitxor(num,bitshift(num,-1));
end
