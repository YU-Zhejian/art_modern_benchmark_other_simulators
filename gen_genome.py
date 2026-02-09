import numpy as np
import sys


if __name__ == '__main__':
    bit_gen = np.random.PCG64()
    rng = np.random.Generator(bit_gen)
    print('>genome')
    for i in range(int(sys.argv[1]) // 1024):
        print(''.join(rng.choice(['A', 'C', 'T', 'G'], size=1024)))

