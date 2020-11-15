function [chi_d, yp_int_dot] = guidance(eta, startWp, endWp, pi_p, yp_int)
    kappa = 0.3;
    crossTrackErr = crosstrackWpt(endWp(1), endWp(2), startWp(1), startWp(2), eta(1), eta(2));
    dist2Wp = norm(endWp - eta(1:2));
    alongTrackDist = sqrt(dist2Wp^2 - crossTrackErr^2); % Pythagoran thm
    Delta = alongTrackDist/5 + 500;
    K_p = 1 / Delta; % Kp = 1 / delta (avoid singularity)
    K_i = kappa * K_p;
    chi_d = pi_p - atan(K_p * crossTrackErr + K_i * yp_int);
    yp_int_dot = Delta*crossTrackErr / (Delta^2 + (crossTrackErr + kappa*yp_int)^2); 
end

