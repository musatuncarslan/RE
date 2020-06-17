function wrapped = some_func(lon)
    wrapped = lon;
    idx_pos = (lon > 0);
    wrapped_pos = mod(lon(idx_pos), 360); 
    wrapped_pos(wrapped_pos == 0) = 360;
    wrapped(idx_pos) = wrapped_pos;
    wrapped(~(idx_pos)) = mod(wrapped(~idx_pos), 360); 
end