#Given a list of avalanches of various sizes, writes them to a .txt file (which some characterizing text)
def avalanche_writing(param_string, characterization, avalanche_sizes):
    avalanche_file = open(param_string + '/avalanche_distrb_' + characterization + '.txt', 'a')
    for i in range(len(avalanche_sizes)):
        avalanche_file.write(str(avalanche_sizes[i]) + '\n')
    avalanche_file.close()
