function stl = coraBreachConvert(cora_stl)

stl = str(cora_stl);

if contains(stl,'U')
    error("bluSTL does not currently (correctly) support until operator")
end

if contains(stl,'X')
    error("bluSTL does not currently support next operator")
end

if contains(stl,'R')
    error("bluSTL does not currently support release operator")
end


stl = replace(stl,'G','alw_');
stl = replace(stl,'F','ev_');
stl = replace(stl,'&','and');
stl = replace(stl,'|','or');
stl = replace(stl,'~','not ');

for k=1:size(cora_stl.variables)
   old = strcat('(\w*)',cora_stl.variables{k},'(\w*)');
   new = strcat(cora_stl.variables{k},'[t]');
   stl = regexprep(stl,old,new);
end

stl = strtrim(stl);

end
