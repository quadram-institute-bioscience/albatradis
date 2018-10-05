def write_plot_file(filename, foward, reverse, p_len):
	with open(filename, 'w') as fh:
		for i in range(0, p_len):
			fh.write(str(foward[i])+' '+str(reverse[i])+"\n")
	return
	