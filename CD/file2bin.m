% function [binData]=file2bin(fname)
% Input(s)
%    fname: is the name of the file like 'otello.txt' or '1.bmp'
% Output(s)
%    binData: is the column vector bit-stream of data in binary
%             Sequence of the bytes is preserved
% Description
%    The function returns a bit-stream of data in the filename "fname"
%    The length of the file is prepended to the bit stream so the first
%    four bytes are the length of the original file
% Considerations
%    The file specified by "fname" must be available in the path
%    The length of the file should be less than (2^32-1) bytes

function [binData]=file2bin(fname)
fid = fopen(fname,'r');
%read the raw file
strData = fread(fid,inf);
%32-bit length of the file
file_len = double(typecast(uint32(length(strData)),'uint8')).';
%matrix of the file in binary (file_len prepended)
mat = dec2bin([file_len ; strData]) - '0'; 
matSize = size(mat);
mat = [zeros(matSize(1),8 - matSize(2)) mat]; %make sure of the 8-bit bytes
binData = reshape(mat.',matSize(1) * 8,1); %bit-stream
fclose(fid);
return