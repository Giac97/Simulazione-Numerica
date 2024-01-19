"""
The following function is used to read the output file from the genetic algorithm simulation in order to obtain the data from it such as the distance and the coordinates (in order beginning to end) of the best path
"""
def read_tsp_output(filename):
    data = []
    with open(filename, 'r') as file:
        block = {}
        for line in file:
            line = line.strip()
            if line.startswith('Generation'):
                if block:
                    data.append(block)
                gen_num = int(line.split()[-1])
                block = {'gen_num': gen_num}
            elif line.startswith('Best distance:'):
                gen_dist = float(line.split()[-1])
                block['gen_dist'] = gen_dist
            elif line.startswith('Best path:'):
                path = []
                for path_line in file:
                    path_line = path_line.strip()
                    if not path_line:
                        break
                    x, y = map(float, path_line.split())
                    path.append((x, y))
                block['path'] = path
    if block:
        data.append(block)
    return data
