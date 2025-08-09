function q=qf(p,u) % PC for translation, fixed u0x 
q=p.u0x(1:p.nu)'*u(1:p.nu); 