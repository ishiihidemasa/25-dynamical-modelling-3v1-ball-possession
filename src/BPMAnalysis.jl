function get_area_timeseries(res::Result)
    x1 = @view res.positionbyid[1:2, :]
    x2 = @view res.positionbyid[3:4, :]
    x3 = @view res.positionbyid[5:6, :]
    return @. abs(
        (x1[1, :] - x3[1, :]) * (x2[2, :] - x3[2, :]) - (x2[1, :] - x3[1, :]) * (x1[2, :] - x3[2, :])
    ) / 2
end

function get_mean_area(res::Result)
    return mean(get_area_timeseries(res))
end
