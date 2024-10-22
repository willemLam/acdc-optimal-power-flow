function setpt_limited = setpoint_limit(setpt,prev_value,limit_max)
    
    limit_low = prev_value - limit_max;
    limit_high = prev_value + limit_max;
    setpt_limited = min(max(setpt,limit_low),limit_high);
end