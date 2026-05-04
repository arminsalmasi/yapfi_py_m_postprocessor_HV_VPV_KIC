import time
import numpy as np

elnames = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'] * 10
phnames = ['P1', 'P2', 'P3', 'P4'] * 10

mf = [np.random.rand(100, 100) for _ in elnames]
mur = [np.random.rand(100, 100) for _ in elnames]
phf = [np.random.rand(100, 100) for _ in phnames]

np.set_printoptions(threshold=np.inf)

def old_way():
    h_str = 'POINT_DATA 10000\n'
    k = 0
    for i in elnames:
        h_str = h_str + 'SCALARS ' + 'mole-fraction(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default' + '\n'  \
                        + np.array_str(mf[k][:]).replace('[', '').replace(']', '')+'\n' \
                        + 'SCALARS ' + 'chemical-potential(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default'+'\n' \
                        + np.array_str(mur[k][:]).replace('[', '').replace(']', '')+'\n'
        k += 1
    k = 0
    for i in phnames:
        h_str = h_str + 'SCALARS ' + 'phase-fraction(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default' + '\n' \
                        + np.array_str(phf[k][:]).replace('[', '').replace(']', '')+'\n'
        k += 1
    return h_str

def new_way():
    parts = ['POINT_DATA 10000\n']
    k = 0
    for i in elnames:
        parts.append('SCALARS ' + 'mole-fraction(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default' + '\n' \
                        + np.array_str(mf[k][:]).replace('[', '').replace(']', '')+'\n' \
                        + 'SCALARS ' + 'chemical-potential(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default'+'\n' \
                        + np.array_str(mur[k][:]).replace('[', '').replace(']', '')+'\n')
        k += 1
    k = 0
    for i in phnames:
        parts.append('SCALARS ' + 'phase-fraction(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default' + '\n' \
                        + np.array_str(phf[k][:]).replace('[', '').replace(']', '')+'\n')
        k += 1
    return "".join(parts)

start = time.time()
r1 = old_way()
end1 = time.time() - start

start = time.time()
r2 = new_way()
end2 = time.time() - start

print(f"Old way: {end1:.4f}s")
print(f"New way: {end2:.4f}s")
assert r1 == r2
