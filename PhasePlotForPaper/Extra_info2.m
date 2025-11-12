function [rot,digg] = Extra_info2(tdiff,numits)

    smallest = eps;
    weights = weighting(numits/2-1);
    
    diffval1 = tdiff(2:numits/2)-tdiff(1:numits/2-1);
    diffval2 = tdiff(numits/2+2:end) -tdiff(numits/2+1:numits-1);
    diffvalnegint1 = find(diffval1<0);
    diffvalnegint2 = find(diffval2<0); 
    diffval1(diffvalnegint1) = diffval1(diffvalnegint1) + 1; 
    diffval2(diffvalnegint2) = diffval2(diffvalnegint2) + 1; 
    
    dig1 = sum(weights.*diffval1);
    dig2 = sum(weights.*diffval2);
    rot1 = sum([diffval1,diffval2]);

    digg = -log10(max(abs(dig2-dig1),smallest));
    rot  = mod(rot1/length(tdiff),1);

end

function G1 = weighting(iter)
    t1 = linspace(1,iter,iter)./iter;
    t1 = (-t1.*(1-t1)).^(-1);
    g1 = exp(t1);
    g1(1,iter) = 0;
    G1 = g1/sum(g1);
end
