import random
import sys


if __name__ == '__main__':
    rng = random.SystemRandom()
    with open('genome.fa', 'w') as f:
        f.write('>genome\n')
        for i in range(int(sys.argv[1])):
            f.write(rng.choice('ACTG') + '\n')
