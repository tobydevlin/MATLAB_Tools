

in.h = 3;
in.Vx = 1;
in.Hrms = 1.18;
in.Hs = in.Hrms*sqrt(2);
in.T = 8;
in.Tp = in.T;
in.d50 = 0.0002;
in.d90 = 0.0003;
in.ks = 0.06;
in.ws = 0.025;

[out1 in1] = vanrijn(in);
[out2 in2] = bijker(in);