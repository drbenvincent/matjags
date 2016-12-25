function A = structsToArrays(S)
% Suppose S is this struct array
%
% S(c).X1(s)
% S(c).X2(s,i)
% S(c).X3(s,i,j)
%
% where s=1:N in all cases
%
% Then we return
% A.X1(c,s)
% A.X2(c,s,i)
% A.X3(c,s,i,j)

C = length(S);
fld = fieldnames(S);
A = [];
for fi=1:length(fld)
	fname = fld{fi};
	tmp = S(1).(fname);
	sz = size(tmp);
	psz = prod(sz);
	data = zeros(C, psz);
	for c=1:C
		tmp = S(c).(fname);
		data(c,:) = tmp(:)';
	end
	if sz(2) > 1 % vector or matrix variable
		data = reshape(data, [C sz]);
	end
	A.(fname) = data;
end
end