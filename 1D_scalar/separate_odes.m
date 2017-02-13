function fv = separate_odes(cv,A,B)
c = cv(1,:) + 1i*cv(2,:);
cp = ODEs_define(c,A,B);
fv = [real(cp);img(cp)];
end