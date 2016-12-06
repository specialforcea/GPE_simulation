function status = isposintscalar(x)
status = (~isempty(x)) && isnumeric(x) && isreal(x) && numel(x)==1 && isfinite(x) && (x>0) && (floor(x) == x);