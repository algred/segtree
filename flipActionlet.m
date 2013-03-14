function [A T] = flipActionlet(A, T)

cols = size(imgs, 2);

for i = 1:length(A)
	A(i).bbox(:, 1) = cols - (A(i).bbox(:, 1) + A(i).bbox(:, 3) - 1) + 1;	
end

for i = 1:length(T)
	T{i}(:, 2) = cols - (T{i}(:, 2) + T{i}(:, 4) - 1) + 1; 
end
	
A = actionletPos(A, T);

end