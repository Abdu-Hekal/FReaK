function blu_stl = coraBreachConvert(cora_stl)

blu_stl = str(cora_stl);

if contains(blu_stl,'U')
    error("bluSTL does not currently (correctly) support until operator")
end

if contains(blu_stl,'X')
    error("bluSTL does not currently support next operator")
end

if contains(blu_stl,'R')
    error("bluSTL does not currently support release operator")
end


blu_stl = replace(blu_stl,'G','alw_');
blu_stl = replace(blu_stl,'F','ev_');
blu_stl = replace(blu_stl,'&','and');
blu_stl = replace(blu_stl,'|','or');
blu_stl = replace(blu_stl,'~','not ');

for k=1:size(cora_stl.variables)
   old = strcat('(\w*)',cora_stl.variables{k},'\s');
   new = strcat(cora_stl.variables{k},'[t]');
   blu_stl = regexprep(blu_stl,old,new);
end

blu_stl = strtrim(blu_stl);

end
