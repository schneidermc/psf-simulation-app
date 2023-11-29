function error = eval_error(stack_simu, stack_exp, gamma)

EXP = stack_exp.^gamma;
SIMU = stack_simu.^gamma;
error = EXP/sum((EXP(:))) - SIMU/sum(SIMU(:)); %for use with lsqnonlin

end