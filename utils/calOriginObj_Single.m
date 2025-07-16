function Obj2 = calOriginObj_Single(Alg, Data, calperobj)

    dec = Alg.result{end}.decs;
    decmin = dec(1,:);
    Obj2 = calperobj(decmin, Data);

end

