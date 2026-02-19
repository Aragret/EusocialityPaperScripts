#!/usr/bin/python3

import random

def random_replace(text, token, replace, num_replacements):
    num_tokens = text.count(token)
    points = [0] + sorted(random.sample(range(1,num_tokens+1),num_replacements)) + [num_tokens+1]
    return replace.join(token.join(text.split(token)[i:j]) for i,j in zip(points,points[1:]))

with open('/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/relax/bootstrap_tree', 'r') as f:
    tree = f.read()

    for i in range(1, 11):
        replaced_tree = random_replace(tree, 'Background', 'Foreground', 3)

        print(replaced_tree.count('Foreground'))

        with open ('/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/relax/bootstrap_{0}/bootstrap.tree'.format(i), 'w') as out:
            out.write(replaced_tree)
