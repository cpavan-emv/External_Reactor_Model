function Tad(gas, T0, P0, comp_string)
    gas.TPX=(T0,P0, comp_string)
    gas.equilibrate("HP")
    return gas.T
end

function SL(gas, T0, P0, comp_string)
    gas.TPX=(T0,P0, comp_string)
    flame=ct.FreeFlame(gas)
    flame.solve(loglevel=0)
    return flame.velocity[1]
end