function source = Seismic_Source(f0,t0,t)
    source = 1.0 * ( 1 - 2 * (pi*f0*(t-t0))^2 ) * exp ( - (pi*f0*(t-t0))^2 );
end 