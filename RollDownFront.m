function loc = RollDownFront(trace,loc,meanval)

if (loc<length(trace)-1)
    if (trace(loc)>trace(loc+1))||(trace(loc)>(meanval))
        loc = RollDownFront(trace,loc+1,meanval);
    end
else
    return;
end

end
