function [tf,amplitude] = time_failure(response,velocity,bar,omega_eq_2,time_out,time)
    
   
    [ndof,nt,ns]=size(response);
    
    for j=1:ndof
        we2(j,:) = interp1(time,omega_eq_2(j,:),time_out,'pchip')';
    end

    amplitude = zeros(size(response));
    tf = zeros(ns,ndof);
    for i=1:ns
        
        for j=1:ndof
            barrier = bar(j);
            
            x = response(j,:,i);
            dx = velocity(j,:,i);
            A = sqrt(x.^2 + dx.^2./we2(j,:));
            amplitude(j,:,i)=A;
            time_aux = time_out(abs(A) > barrier);

            if numel(time_aux)==0
                tf(i,j) = NaN;
                %warning('Check the barrier');
            else
                tf(i,j) = time_aux(1);
            end
        end

    end