clc
clear all
syms Gast I  Krr Krl Klr Kll  GastT


A =  [Gast;I]
At = [GastT I];

K = [Krr Krl; Klr Kll]

KllGEN = At*K*A

simple(KllGEN)
