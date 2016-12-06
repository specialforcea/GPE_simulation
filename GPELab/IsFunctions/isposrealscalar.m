function status = isposrealscalar(x)
status = (~isempty(x)) && isnumeric(x) && isreal(x) && numel(x)==1 && isfinite(x) && (min(x,0)==0);

