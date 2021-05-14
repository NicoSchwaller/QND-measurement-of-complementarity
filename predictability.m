function p=predictability(qubit,rho)
    rhoA=PartialTrace(rho,[2],[2,2]);
    rhoB=PartialTrace(rho,[1],[2,2]);
    if qubit=='A'
        p=abs(rhoA(1,1)-rhoA(2,2));
    end
    if qubit=='B'
        p=abs(rhoB(1,1)-rhoB(2,2));
    end
end