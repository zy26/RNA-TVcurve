function RunRNAplotCmd(command, input)
	global globalpath
	if (size(globalpath) == 0)
		path = '';
	else
		path = strcat(globalpath, filesep);
	end
	cmd = strcat(path, command, ' < ', input);
	system(cmd)
end