#function per media ed errore

function mean_func(d_slope,v_slope)

    mean_slope_v=mean(v_slope)
    err_slope_v=std(v_slope)
    mean_slope_d=mean(d_slope)
    err_slope_d=std(d_slope)

    return mean_slope_d,err_slope_d,mean_slope_v,err_slope_v
end
