% A supplement to red_tide, this function makes the H matrix from a given F
% vector without doing all the other stuff that red_tide does.
% 
% IN:   T = vector of times
% IN:   F = vector of frequencies
% 
% OUT:  H = Regressor matrix built from F such that:
%           H = [sin(2*pi*F(1)*T) cos(2*pi*F(1)*T) ... ]

function H = H_make(T,F)

if isrow(T)
    T = T';
else
end

H = zeros(length(T),2*length(F));

for i = 1:length(F)
    H(:,2*i - 1) = sin(2*pi*F(i)*T);
    H(:,2*i    ) = cos(2*pi*F(i)*T);
end

end