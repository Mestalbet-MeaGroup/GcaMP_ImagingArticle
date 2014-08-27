%---Recursive Roll down hill Functions---%
function loc = RollDownBack(trace,loc,meanval)

if (loc>1)
    if (trace(loc)>trace(loc-1))||(trace(loc)>(meanval))
        loc = RollDownBack(trace,loc-1,meanval);
    end
else
    return;
end

end
