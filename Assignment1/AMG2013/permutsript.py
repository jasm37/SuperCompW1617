from itertools import combinations
file = open('permutw','w')
a = ('-march=native','-fomit-frame-pointer','floop-block','-floop-interchange','-floop-strip-mine','-funroll-loops','-flto');
file.write('\n'.join('\n'.join(' '.join(str(i) for i in c) for c in combinations(a, i)) for i in range(1, len(a)+1)));
