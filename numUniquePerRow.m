function y = numUniquePerRow(x)

%Sort each row
srtX = sortrows(x.');

y = sum(diff(srtX, 1,1)>=1,1);
end