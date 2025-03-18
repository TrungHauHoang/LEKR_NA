function [F,MVMs] = vectorF_1D(Jacobian,U)

F = Jacobian*U;
MVMs = 1;
