A=[
-3	1	0.0497870683678639
-2.90000000000000	1	0.0550232200564072
-2.80000000000000	1	0.0608100626252180
-2.70000000000000	1	0.0672055127397498
-2.60000000000000	1	0.0742735782143339
-2.50000000000000	1	0.0820849986238988
-2.40000000000000	1	0.0907179532894125
-2.30000000000000	1	0.100258843722804
-2.20000000000000	2	0.0554015791811669
-2.10000000000000	2	0.0612282141264910
-2	2	0.0676676416183064
-1.90000000000000	2	0.0747843096113175
-1.80000000000000	2	0.0826494441107933
-1.70000000000000	2	0.0913417620263673
-1.60000000000000	2	0.100948258997328
-1.50000000000000	3	0.0743767200494766
-1.40000000000000	3	0.0821989879805355
-1.30000000000000	3	0.0908439310113375
-1.20000000000000	3	0.100398070637401
-1.10000000000000	3	0.110957027899360
-1	4	0.0919698602928606
-0.900000000000000	4	0.101642414935150
-0.800000000000000	4	0.112332241029305
-0.700000000000000	5	0.0993170607582819
-0.600000000000000	5	0.109762327218805
-0.500000000000000	5	0.121306131942527
-0.400000000000000	6	0.111720007672607
-0.300000000000000	6	0.123469703446953
-0.200000000000000	6	0.136455125512997
-0.0999999999999996	7	0.129262488290851
0	7	0.142857142857143
0.100000000000000	8	0.138146364759456
0.200000000000000	9	0.135711417573352
0.300000000000000	9	0.149984311952889
0.400000000000000	10	0.149182469764127
0.500000000000000	11	0.149883751881830
0.600000000000000	11	0.165647163671864
0.700000000000000	12	0.167812725622540
0.800000000000000	13	0.171195456037882
0.900000000000000	14	0.175685936511211
1	15	0.181218788563936
1.10000000000000	16	0.187760376496652
1.20000000000000	18	0.184450940152030
1.30000000000000	19	0.193120877243118
1.40000000000000	20	0.202759998342234
1.50000000000000	22	0.203713139560821
1.60000000000000	23	0.215349235843266
1.70000000000000	25	0.218957895669088
1.80000000000000	27	0.224061017200480
1.90000000000000	29	0.230548084216527
2	31	0.238356648352602
2.10000000000000	33	0.247459694320232
2.20000000000000	36	0.250694819428726
2.30000000000000	39	0.255748268072173
2.40000000000000	41	0.268857960503454
2.50000000000000	44	0.276874862743261
2.60000000000000	48	0.280494542395869
2.70000000000000	51	0.291759445585742
2.80000000000000	55	0.298993577656310
2.90000000000000	59	0.308036362193950
3	64	0.313836514424807
];

for i=1:size(A,1)
    foldername=strcat(' RMI_lnKA_',num2str(A(i,1),'%g'));
    s=strcat('mkdir',foldername);
    system(s);
	foldername=strcat('RMI_lnKA_',num2str(A(i,1),'%g'));
    filename=strcat(foldername,'/value_start.m');
if A(i,1)<=2
copyfile('value_start.m_bak',filename);
else
copyfile('value_start.m_bak2',filename);   
end
    fid = fopen(filename,'rt+');
    frewind(fid);
    fprintf(fid,'A=%g;\n',A(i,3));
	fprintf(fid,'K=%g;\n',A(i,2));
    fclose(fid);
end