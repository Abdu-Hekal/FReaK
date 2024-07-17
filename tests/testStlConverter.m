x = stl('my_var',2);
cora_stl = globally(x(1)<2 & x(2)>3, interval(0,1));
str = coraBreachConvert(cora_stl)


