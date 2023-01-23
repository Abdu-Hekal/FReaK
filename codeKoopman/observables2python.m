function text = observables2python(f,n)

    x = sym('x',[n,1]);
    g = f(x);

    tmp = string(vpa(g));
    text = '';
    
    for i = 1:length(tmp)
        text = [text,newline,tmp{i},','];
    end
    text = text(1:end-1);
    
    text = strrep(text,'cos','np.cos');
    text = strrep(text,'^','**');
    
    for i = 1:n
       text = strrep(text,['x',num2str(i)],['x','[',num2str(i-1),']']); 
    end
    
    text = ['return np.array([',newline,text,'])'];
end