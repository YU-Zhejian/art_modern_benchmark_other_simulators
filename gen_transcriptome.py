import numpy as np
import sys


if __name__ == '__main__':
    bit_gen = np.random.PCG64()
    rng = np.random.Generator(bit_gen)
    for i in range(int(sys.argv[1])):
        print(f'>transcript{i}')
        print(''.join(rng.choice(['A', 'C', 'T', 'G'], int(sys.argv[2]))))
