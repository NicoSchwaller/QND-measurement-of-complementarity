function v=visibility(qubit,rho)
    rhoA=PartialTrace(rho,[2],[2,2]);
    rhoB=PartialTrace(rho,[1],[2,2]);
    if qubit=='A'
        v=2*abs(rhoA(1,2));
    end
    if qubit=='B'
        v=2*abs(rhoB(1,2));
    end

end