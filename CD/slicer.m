% function [q] = slicer(a , M)
% Input(s)
%    a: the complex input to be sliced/quantized
%    M: constellation size
%             
% Output(s)
%    q: detected symbol
% Description
%   Given the rectangular constellation size, the nearest point
%   in the constellation is returned
% Considerations
%   M must be power of 4 and greater than or equal to 4
%   The constellation coordinates are odd numbers

function [q] = slicer(a , M)

if M == 4
        real_a = real(a);
        imag_a = imag(a);
        if real_a > 0
            if imag_a > 0
                q = 1 + i;
            else
                q = 1 - i;
            end
        else
            if imag_a > 0
                q = -1 + i;
            else
                q = -1 - i;
            end
        end
else
        real_a = real(a);
        imag_a = imag(a);
        if real_a > 0
            if imag_a > 0
                center = (1 + i) * log2(M / 4);
            else
                center = (1 - i) * log2(M / 4);
            end
        else
            if imag_a > 0
                center = (-1 + i) * log2(M / 4);
            else
                center = (-1 - i) * log2(M / 4);
            end
        end
        q = center + quan(a - center, M / 4);
end

return
