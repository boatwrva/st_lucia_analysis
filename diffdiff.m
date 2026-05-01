function d = diffdiff(mat,dim)
n = ndims(mat);
inds_1 = '(';
for n = 1:n
    if n==dim
        inds_1 = [inds_1 '1:end-1,'];
    else
        inds_1 = [inds_1 ':,'];
    end
end
inds_1 = [inds_1(1:end-1) ')'];

inds_2 = '(';
for n = 1:n
    if n==dim
        inds_2 = [inds_2 '2:end,'];
    else
        inds_2 = [inds_2 ':,'];
    end
end
inds_2 = [inds_2(1:end-1) ')'];

eval_str = ['diff(mat' inds_1 ',1,' num2str(dim) ')-flip(diff(flip(mat' inds_2 ',' num2str(dim) '),1,' num2str(dim) ')' ',' num2str(dim) ');'];

d = eval(eval_str);

end
% d = (diff(mat(1:end-1))-fliplr(diff(fliplr(heave_HDSS(2:end)))));
