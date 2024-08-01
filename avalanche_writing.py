#Given a list of avalanches of various sizes, writes them to a .txt file (which some characterizing text that specific time_window, T_min, and T_max)
def avalanche_writing(param_string, characterizing_text, avalanche_sizes, data_type):
    avalanche_file = open(param_string + '/avalanche_distrb_' + characterizing_text + '.txt', 'a')
    for i in range(len(avalanche_sizes)):
        avalanche_file.write(str(avalanche_sizes[i]) + '\n')
    avalanche_file.close()
