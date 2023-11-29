%Convert QPSK Symbols to bits
function info_bits=QPSK2bits(info_symbols)
    L = length(info_symbols);
    info_bits=zeros(2*L,1);
    bit_index = 1;
    for i=1:L
        detected_symbol = slicer(info_symbols(i), 4);
        switch detected_symbol
            case -1 - j
                info_bits(bit_index) = 0;
                bit_index = bit_index + 1;
                info_bits(bit_index) = 0;
                bit_index = bit_index + 1;
            case -1 + j
                info_bits(bit_index) = 0;
                bit_index = bit_index + 1;
                info_bits(bit_index) = 1;
                bit_index = bit_index + 1;
            case  1 + j
                info_bits(bit_index) = 1;
                bit_index = bit_index + 1;
                info_bits(bit_index) = 1;
                bit_index = bit_index + 1;
            case  1 - j
                info_bits(bit_index) = 1;
                bit_index = bit_index + 1;
                info_bits(bit_index) = 0;
                bit_index = bit_index + 1;
        end
    end
end
