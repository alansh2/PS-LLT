function exportWing(filename,vert,sect)
%EXPORTWING  Export wing configuration to text file
fid = fopen(filename,'w');
fprintf(fid,'%13s %13s %13s %13s %13s %13s %13s\r\n','xle','xte','y','z','twist(deg)','a0','alfZL(deg)');
fprintf(fid,'%13.8f %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f\r\n',[vert [sect;[0 0 0]]].');
fclose(fid);