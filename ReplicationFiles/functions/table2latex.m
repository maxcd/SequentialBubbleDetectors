function table2latex(tt_list, flformat, cformat) 

numel(tt_list)
nrows = []; ncols = [];
cnames ={};
for it = 1:numel(tt_list)
    tt = tt_list{it};
    [a, b] = size(tt);
    nrows(it) = a; 
    ncols(it) = b;
    cnames = [cnames tt.Properties.VariableNames];
end

if ~all( nrows-nrows(1) == 0)
    disp('tables must have identical number of rows.')
else
    nrows = nrows(1);
end

if nargin < 4
    cformat = 'p{2cm}';
end

tdata = {} ;str_head = []; table_names = {}; str_col = [];
for it = 1:numel(tt_list)
    if it>1
        str_head = [str_head, ' & '];
        str_col = [str_col, ' & '];
    end
    tt = tt_list{it};
    if it>1
        tdata{it} = [tdata, tt{:,:}];
    elseif it==1
        tdata = tt{:,:};
    end
    str_head = [str_head ,'\\multicolumn{',num2str(size(tt, 2)) , '}{c}{\\textbf{%s}}'];
    
    [nr, nc] = size(tt);

    str_tmp =  [flformat{it}, repmat( ['& ', flformat{it}, ' '], 1, nc-1)];
    str_col = [str_col, str_tmp];
    table_names= [table_names, ['Table', num2str(it)]];
end
str_head = [str_head, '\\\\ \n'];
str_col = [str_col, '\\\\ \n'];


ncols = sum(ncols);

if nargin < 2
    flformat = '%5.4f'; 
end
%ncols = length(cnames);
%nrows = length(lam0seq)
%strdgp = ['&\\multicolumn{3}{c}{(i) constant mean ($\\mu_0=y_0=1$), $\\beta=0$}&\\multicolumn{3}{c|}{(ii) trend ($\\mu_0=y_0=0$, $\\beta=1$)}  \\\\ \n '];


cc = [' \\multicolumn{1}{ ', cformat, '}{\\centering $%s$} '];
str_names =[ cc, repmat(['&', cc], 1, ncols-1), '\\\\ \\hline \n'];
str1 = ['%5.1f', repmat( ['& ', flformat, ' '], 1, ncols-1) ,   ' \\\\ \n'];


% print the latex table code
fprintf(str_head, table_names{:})
fprintf(str_names, cnames{:})
for irow=1:nrows
    fprintf(str_col, tdata(irow,:))
end

end

function out = getVarName(var)
    out = inputname(1);
end