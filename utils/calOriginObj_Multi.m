function [HV2,obj_f1] = calOriginObj_Multi(Alg, Data, calperobj1, points)

    decp_f1 = findf1(Alg.result{end});
    obj_f1 = calperobj1(decp_f1,Data);
    HV2 = HVtest(obj_f1, points);

end

