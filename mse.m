function mean_squared_error = mse(symbols)
    mean_squared_error = 0;
    for k=1:length(symbols)
        sliced_symbol = slicer(symbols(k), 4);
        mean_squared_error = mean_squared_error + (sliced_symbol-symbols(k))^2;
    end
    mean_squared_error = mean_squared_error/length(symbols);
end