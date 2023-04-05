function blu_stl = cora_blu_stl_convert(cora_stl)

blu_stl = evalc('disp(cora_stl)');

if contains(blu_stl,'U')
    error("bluSTL does not currently (correctly) support until operator")
end

blu_stl = replace(blu_stl,'G','alw_');
blu_stl = replace(blu_stl,'F','ev_');
blu_stl = replace(blu_stl,'&','and');
blu_stl = replace(blu_stl,'|','or');
blu_stl = replace(blu_stl,'~','not ');
blu_stl = regexprep(blu_stl,'x(\w*)','x$1(t)');

for k=1:size(cora_stl.variables)
   old = strcat('(\w*)',cora_stl.variables{k},'\s');
   new = strcat('x',string(k),'(t)');
   blu_stl = regexprep(blu_stl,old,new);
end

blu_stl = strtrim(blu_stl);

end
