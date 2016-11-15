function RunCmd(command, input, output)
	global globalpath
	if (size(globalpath) == 0)
		path = '';
	else
		path = strcat(globalpath, filesep);
	end
	cmd = strcat(path, command, ' < ', input, ' > ', output);
	system(cmd)
end