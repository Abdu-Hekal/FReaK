function pgon = convertToPolygon(R,dim)
% convert a reachable set object to a polygon
        
    pgon = []; list = cell(length(R.timePoint.set),1);

    for i = 1:length(R.timePoint.set)-1
       
        int1 = interval(R.timePoint.set{i}); int1 = int1(dim);
        int2 = interval(R.timePoint.set{i+1}); int2 = int2(dim);
        t1 = R.timePoint.time{i};
        t2 = R.timePoint.time{i+1};
        
        V = [t1 t1 t2 t2; ...
             infimum(int1) supremum(int1) supremum(int2) infimum(int2)];
        
        list{i} = polygon(V(1,:),V(2,:));
    end
    
    while length(list) > 1
       cnt1 = 1; cnt2 = 1;
       while cnt2 < length(list)
          list{cnt1} = list{cnt2} | list{cnt2+1};
          cnt1 = cnt1 + 1; cnt2 = cnt2 + 2;
       end
       if cnt2 - 1 < length(list)
          list{cnt1} = list{cnt2};
          cnt1 = cnt1 + 1;
       end
       list = list(1:cnt1-1);
    end
    
    pgon = list{1};
end