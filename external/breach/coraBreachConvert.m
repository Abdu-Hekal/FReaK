function stl = coraBreachConvert(cora_stl)

stl = str(cora_stl);

if contains(stl,'U')
    error("KF does not currently (correctly) support until operator")
end

if contains(stl,'X')
    error("KF does not currently support next operator")
end

if contains(stl,'R')
    error("KF does not currently support release operator")
end


stl = replace(stl,'G','alw_');
stl = replace(stl,'F','ev_');
stl = replace(stl,'&','and');
stl = replace(stl,'|','or');
stl = replace(stl,'~','not ');

for k=1:size(cora_stl.variables)
   old = strcat(cora_stl.variables{k},'(?![0-9])');
   new = strcat(cora_stl.variables{k},'[t]');
   stl = regexprep(stl,old,new);
end

stl = strtrim(stl);
end
