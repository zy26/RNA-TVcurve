function savedouble(filename, header1, header2, number)
fid = fopen(filename,'w');
fprintf(fid, 'The distance score between "%s" and "%s": \n', header1, header2);
fprintf(fid, '%f\n', number);
fclose(fid);