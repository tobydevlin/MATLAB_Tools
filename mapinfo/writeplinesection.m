function writeplinesection(themiffile,themidfile,thex,they,thestart,theend)

fprintf(themiffile,'Pline %d\r\n',length(thex));
fprintf(themiffile,'%11.2f %11.2f\r\n',[thex,they]');
fprintf(themiffile,'      Pen (1,2,0)\r\n');

fprintf(themidfile,'"ARC",%d,%d\r\n',thestart,theend);

end 



